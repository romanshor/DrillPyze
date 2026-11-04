"""
Numba time-stepping kernel for the coupled axial-torsional drillstring model.

See ``docs/derivation/README.md`` for the equations.  The kernel is written over plain
arrays so it compiles to tight machine code; :class:`DrillPyze.Simulation`
owns the arrays and provides the user-facing API.
"""

import numpy as np
from numba import njit

from .geometry import normal_force

# --------------------------------------------------------------------------
# indices into the scalar state vector ``sc``
S_T = 0           # simulation time (s)
S_WTD = 1         # top drive speed (rad/s)
S_TAUM = 2        # motor torque (N m)
S_TAU0 = 3        # pipe torque at top (N m)
S_INT = 4         # PI integrator state
S_LPTAU = 5       # low-pass state of tau_0 (ZTorque)
S_HPLP = 6        # internal low-pass state of the high-pass filter
S_WFILT = 7       # filtered top-drive speed
S_WSP = 8         # rate-limited speed set point (rad/s)
S_BITD = 9        # bit depth (m MD)
S_HOLED = 10      # hole depth (m MD)
S_WBIT = 11       # bit speed (rad/s)
S_TOB = 12        # torque on bit (N m)
S_VBIT = 13       # bit axial velocity (m/s)
S_WOB = 14        # weight on bit (N)
S_VBLK = 15       # block velocity (m/s, + down)
S_HKLD = 16       # hookload (N, effective tension at surface)
S_AD_INT = 17     # autodriller integrator
S_BLKPOS = 18     # block travel since start (m, + down)
S_BITSTUCK = 19   # 1 if bit rotationally stuck
N_SCALARS = 20

# indices into the float parameter vector ``pf``
P_DT = 0
P_KT = 1
P_KA = 2
P_WC = 3
P_VC = 4
P_ITD = 5
P_KP = 6
P_KI = 7
P_ZC = 8
P_THP = 9
P_TLP = 10
P_TEN = 11
P_TAUMAX = 12
P_RAMP = 13
P_RBIT = 14
P_MUBS = 15
P_MUBK = 16
P_TSTICK = 17     # bit sticking torque at zero WOB (N m)
P_KROP = 18
P_W0 = 19
P_QS_TOP = 20      # static compression at the top face (N)
P_AD_KP = 21       # autodriller gains (m/s per N)
P_AD_KI = 22
P_AD_VMAX = 23     # autodriller max block speed (m/s)
P_HKLD_REF = 24    # off-bottom reference hookload for autodriller (N)
N_PF = 25

# indices into the integer parameter vector ``pi``
I_CTRL = 0         # 0 internal PI/ZTorque, 1 external motor torque, 2 prescribed speed
I_AXIAL = 1        # 0 prescribed block velocity, 1 autodriller (WOB)
I_COUPLE = 2       # 1: 2-D Coulomb friction on (v, omega r) with live-tension F_N (default)
                   # 0: independent axial/torsional friction with static F_N
N_PI = 3

# recorded channels
REC_CHANNELS = ('time', 'omega_td', 'motor_torque', 'pipe_torque', 'omega_bit',
                'tob', 'v_block', 'hookload', 'v_bit', 'wob', 'bit_depth',
                'hole_depth', 'omega_sp')
N_REC = len(REC_CHANNELS)


@njit(cache=True, fastmath=False)
def _face(aL, bR, ZL, ZR):
    """Godunov face solution: returns (velocity, force) from incoming invariants."""
    f = (aL - bR) / (1.0 / ZL + 1.0 / ZR)
    return aL - f / ZL, f


@njit(cache=True)
def _sign(x):
    if x > 0.0:
        return 1.0
    elif x < 0.0:
        return -1.0
    return 0.0


@njit(cache=True)
def _bit(alpha_t, alpha_a, Zt, Za, on_bottom, w_prev, stuck_prev, pf):
    """
    Bit boundary: the axial-torsional coupling point.

    Coulomb-type bit-rock interaction:
      * stuck (|omega_b| <= omega_c):  |TOB| <= T_s0 + mu_bs * r_b * WOB
      * sliding:                       TOB = mu_bk * r_b * WOB * sign(omega_b)
    Penetration: ROP = k_r * max(WOB - W_0, 0) * |omega_b|, unilateral contact.

    Returns (omega_b, tob, v_b, wob, stuck).
    """
    if (not on_bottom) or alpha_a <= 0.0:
        # off bottom, or the string is lifting the bit
        return alpha_t, 0.0, alpha_a, 0.0, 0.0

    r_b = pf[P_RBIT]
    mus = pf[P_MUBS]
    muk = pf[P_MUBK]
    ts0 = pf[P_TSTICK]
    kr = pf[P_KROP]
    w0 = pf[P_W0]
    wc = pf[P_WC]

    w = alpha_t if alpha_t > 0.0 else 0.0
    wob = 0.0
    v = 0.0
    tob = 0.0
    stuck = 0.0
    for _ in range(8):
        # axial: penetration law and characteristic relation v_b = alpha_a - WOB/Za
        q_lin = Za * alpha_a                     # WOB if the bit does not advance
        aw = abs(w)
        if q_lin <= w0 or kr * aw == 0.0:
            wob = q_lin
            v = 0.0
        else:
            wob = (alpha_a + kr * w0 * aw) / (1.0 / Za + kr * aw)
            v = kr * (wob - w0) * aw
        # torsional: Coulomb bit friction; omega_b = alpha_t - TOB/Zt
        tau_stop = Zt * alpha_t                  # TOB that holds the bit at rest
        if (abs(w_prev) <= wc or stuck_prev > 0.5) and abs(tau_stop) <= ts0 + mus * r_b * wob:
            tob = tau_stop
            w_new = 0.0
            stuck = 1.0
        else:
            tob = _sign(alpha_t) * min(muk * r_b * wob, abs(tau_stop))
            w_new = alpha_t - tob / Zt
            stuck = 0.0
        if abs(w_new - w) < 1e-12:
            w = w_new
            break
        w = w_new
    return w, tob, v, wob, stuck


@njit(cache=True)
def _coulomb_2d(va, vt, m_a, m_t, c):
    """
    Implicit (backward-Euler) 2-D Coulomb friction on one cell.

    Given the predicted sliding velocity (va, vt) = (v, omega r), the masses
    per length in each direction (m_a = rho A, m_t = rho J / r^2) and the
    friction impulse per length c = dt * mu * F_N, return the new velocity

        V = V* - (dt/M) * mu F_N * V/|V|,

    i.e. friction opposes the *new* sliding direction exactly.  If the
    momentum |M V*| is within c the cell sticks (V = 0).  Otherwise
    V_i = m_i V*_i s / (m_i s + c), where the new speed s = |V| solves

        g(s) = sum_i (m_i V*_i / (m_i s + c))^2 - 1 = 0,

    a convex decreasing function solved by Newton's method from s = |V*|.
    """
    pa = m_a * va
    pt = m_t * vt
    if pa * pa + pt * pt <= c * c:
        return 0.0, 0.0
    if c <= 0.0:
        return va, vt
    s = np.sqrt(va * va + vt * vt)
    for _ in range(50):
        da = m_a * s + c
        dt_ = m_t * s + c
        ra = pa / da
        rt = pt / dt_
        g = ra * ra + rt * rt - 1.0
        dg = -2.0 * (ra * ra * m_a / da + rt * rt * m_t / dt_)
        s_new = s - g / dg
        if s_new < 0.0:
            s_new = 0.0
        if abs(s_new - s) <= 1e-14 * (s + 1e-30):
            s = s_new
            break
        s = s_new
    return va * m_a * s / (m_a * s + c), vt * m_t * s / (m_t * s + c)


@njit(cache=True, fastmath=True)
def advance(nsteps, sc, pf, pi,
            # grid
            h, Kt, mt, Zt, Ka, ma, Za, rc, wb_len, kx, ky, kz, zx, zy, zz,
            mus, muk, qs_c, fn0,
            # state
            w, tau, v, qd,
            # inputs (held constant over the call)
            w_sp_target, v_blk_in, taum_ext, wtd_ext, wob_sp,
            # recording
            rec, rec_idx, rec_every, rec_phase,
            # scratch
            wf, tf, vf, qf):
    """
    Advance the model ``nsteps`` time steps.  State arrays are updated in
    place.  Every ``rec_every`` steps a row is written into ``rec`` starting at
    ``rec_idx``; returns (new rec_idx, new rec_phase).
    """
    n = h.shape[0]
    dt = pf[P_DT]
    kt = pf[P_KT]
    ka = pf[P_KA]
    wc = pf[P_WC]
    vc = pf[P_VC]
    I = pf[P_ITD]
    ctrl = pi[I_CTRL]
    axial_mode = pi[I_AXIAL]
    coupled = pi[I_COUPLE] == 1

    # filter coefficients (exact discretisation of first-order lags)
    a_en = np.exp(-dt / pf[P_TEN]) if pf[P_TEN] > 0.0 else 0.0
    a_lp = np.exp(-dt / pf[P_TLP]) if pf[P_TLP] > 0.0 else 0.0
    a_hp = np.exp(-dt / pf[P_THP]) if pf[P_THP] > 0.0 else 0.0

    nrec = rec.shape[0]

    for step in range(nsteps):
        # ---------------------------------------------------- top drive input
        wtd = sc[S_WTD]
        if ctrl == 0:
            # rate-limited set point
            dw = w_sp_target - sc[S_WSP]
            lim = pf[P_RAMP] * dt
            if pf[P_RAMP] > 0.0 and abs(dw) > lim:
                sc[S_WSP] += lim * _sign(dw)
            else:
                sc[S_WSP] = w_sp_target
            # measurements (speed filter, torque band-pass for ZTorque)
            sc[S_WFILT] = a_en * sc[S_WFILT] + (1.0 - a_en) * wtd
            sc[S_LPTAU] = a_lp * sc[S_LPTAU] + (1.0 - a_lp) * sc[S_TAU0]
            sc[S_HPLP] = a_hp * sc[S_HPLP] + (1.0 - a_hp) * sc[S_LPTAU]
            tau_bp = sc[S_LPTAU] - sc[S_HPLP]
            e = sc[S_WSP] - sc[S_WFILT] - pf[P_ZC] * tau_bp
            u = pf[P_KP] * e + pf[P_KI] * (sc[S_INT] + e * dt)
            tmax = pf[P_TAUMAX]
            if tmax > 0.0 and abs(u) > tmax:
                u = tmax * _sign(u)          # saturate, freeze integrator
            else:
                sc[S_INT] += e * dt
            sc[S_TAUM] = u
        elif ctrl == 1:
            sc[S_TAUM] = taum_ext
            sc[S_WSP] = w_sp_target

        # ------------------------------------------------ top drive dynamics
        beta_t1 = w[0] - tau[0] / Zt[0]
        if ctrl == 2:
            w_new_td = wtd_ext
        else:
            w_new_td = (I * wtd / dt + sc[S_TAUM] + Zt[0] * beta_t1) / (I / dt + Zt[0])
        tau0 = Zt[0] * (w_new_td - beta_t1)
        if ctrl == 2:
            sc[S_TAUM] = tau0 + I * (w_new_td - wtd) / dt
            sc[S_WSP] = wtd_ext
        sc[S_WTD] = w_new_td
        sc[S_TAU0] = tau0
        wf[0] = w_new_td
        tf[0] = tau0

        # ---------------------------------------------------- drawworks input
        if axial_mode == 1:
            wob_surf = pf[P_HKLD_REF] - sc[S_HKLD]
            err = wob_sp - wob_surf
            vb = pf[P_AD_KP] * err + pf[P_AD_KI] * (sc[S_AD_INT] + err * dt)
            vmax = pf[P_AD_VMAX]
            if vb > vmax:
                vb = vmax
            elif vb < -vmax:
                vb = -vmax
            else:
                sc[S_AD_INT] += err * dt
            sc[S_VBLK] = vb
        else:
            sc[S_VBLK] = v_blk_in
        beta_a1 = v[0] - qd[0] / Za[0]
        q0 = Za[0] * (sc[S_VBLK] - beta_a1)
        vf[0] = sc[S_VBLK]
        qf[0] = q0
        sc[S_HKLD] = -(q0 + pf[P_QS_TOP])

        # --------------------------------------------------- interior faces
        for j in range(n - 1):
            a_t = w[j] + tau[j] / Zt[j]
            b_t = w[j + 1] - tau[j + 1] / Zt[j + 1]
            wf[j + 1], tf[j + 1] = _face(a_t, b_t, Zt[j], Zt[j + 1])
            a_a = v[j] + qd[j] / Za[j]
            b_a = v[j + 1] - qd[j + 1] / Za[j + 1]
            vf[j + 1], qf[j + 1] = _face(a_a, b_a, Za[j], Za[j + 1])

        # ---------------------------------------------------------- bit face
        alpha_t = w[n - 1] + tau[n - 1] / Zt[n - 1]
        alpha_a = v[n - 1] + qd[n - 1] / Za[n - 1]
        on_bottom = sc[S_BITD] >= sc[S_HOLED] - 1e-9
        wb_, tob, vb_, wob, stuck = _bit(alpha_t, alpha_a, Zt[n - 1], Za[n - 1],
                                         on_bottom, sc[S_WBIT], sc[S_BITSTUCK], pf)
        wf[n] = wb_
        tf[n] = tob
        vf[n] = vb_
        qf[n] = wob
        sc[S_WBIT] = wb_
        sc[S_TOB] = tob
        sc[S_VBIT] = vb_
        sc[S_WOB] = wob
        sc[S_BITSTUCK] = stuck

        # -------------------------------------------- cell update + sources
        for j in range(n):
            w_old = w[j]
            v_old = v[j]
            # flux update
            tau[j] -= dt * Kt[j] / h[j] * (wf[j + 1] - wf[j])
            ws = w[j] - dt / (mt[j] * h[j]) * (tf[j + 1] - tf[j])
            qd[j] -= dt * Ka[j] / h[j] * (vf[j + 1] - vf[j])
            vs = v[j] - dt / (ma[j] * h[j]) * (qf[j + 1] - qf[j])
            # implicit viscous damping
            ws /= (1.0 + dt * kt)
            vs /= (1.0 + dt * ka)
            r = rc[j]
            if coupled:
                # ---- 2-D Coulomb inclusion on the sliding velocity (v, omega r)
                #      with the normal force from the live tension: axial load
                #      changes torque and rotation changes drag
                T = -(qs_c[j] + qd[j])
                fn = normal_force(T, wb_len[j], kx[j], ky[j], kz[j], zx[j], zy[j], zz[j])
                if (v_old / vc) ** 2 + (w_old / wc) ** 2 <= 1.0:
                    mu = mus[j]
                else:
                    mu = muk[j]
                va, vt = _coulomb_2d(vs, ws * r, ma[j], mt[j] / (r * r), dt * mu * fn)
                v[j] = va
                w[j] = vt / r
            else:
                # ---- independent 1-D Coulomb inclusions (Aarsnes & Shor eqs. 4, 20-22)
                fn = fn0[j]
                # torsion: torque per length needed to stop the cell this step
                t_stop = mt[j] * ws / dt
                if abs(w_old) <= wc:
                    cap = mus[j] * fn * r
                else:
                    cap = muk[j] * fn * r
                if t_stop > cap:
                    t_stop = cap
                elif t_stop < -cap:
                    t_stop = -cap
                w[j] = ws - dt * t_stop / mt[j]
                # axial: force per length needed to stop the cell this step
                f_stop = ma[j] * vs / dt
                if abs(v_old) <= vc:
                    cap = mus[j] * fn
                else:
                    cap = muk[j] * fn
                if f_stop > cap:
                    f_stop = cap
                elif f_stop < -cap:
                    f_stop = -cap
                v[j] = vs - dt * f_stop / ma[j]

        # -------------------------------------------------- bookkeeping
        sc[S_T] += dt
        sc[S_BITD] += sc[S_VBIT] * dt
        sc[S_BLKPOS] += sc[S_VBLK] * dt
        if sc[S_BITD] > sc[S_HOLED]:
            sc[S_HOLED] = sc[S_BITD]

        rec_phase += 1
        if rec_phase >= rec_every:
            rec_phase = 0
            if rec_idx < nrec:
                rec[rec_idx, 0] = sc[S_T]
                rec[rec_idx, 1] = sc[S_WTD]
                rec[rec_idx, 2] = sc[S_TAUM]
                rec[rec_idx, 3] = sc[S_TAU0]
                rec[rec_idx, 4] = sc[S_WBIT]
                rec[rec_idx, 5] = sc[S_TOB]
                rec[rec_idx, 6] = sc[S_VBLK]
                rec[rec_idx, 7] = sc[S_HKLD]
                rec[rec_idx, 8] = sc[S_VBIT]
                rec[rec_idx, 9] = sc[S_WOB]
                rec[rec_idx, 10] = sc[S_BITD]
                rec[rec_idx, 11] = sc[S_HOLED]
                rec[rec_idx, 12] = sc[S_WSP]
                rec_idx += 1
    return rec_idx, rec_phase
