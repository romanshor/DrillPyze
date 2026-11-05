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
S_BITSTUCK = 19   # 1 if bit rotationally sticking
S_PHIB = 20       # bit rotation angle (rad), regenerative model
S_WINT = 21       # integral of |omega_bit| dt (rad), for averages over slow updates
N_SCALARS = 22

# indices into the float parameter vector ``pf``
P_DT = 0
P_KT = 1           # unused (damping is per cell, see advance())
P_KA = 2           # unused
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
P_LEVS = 14       # bit torque per unit WOB while sticking (m): TOB <= P_LEVS * WOB
P_LEVK = 15       # bit torque per unit WOB while sliding (m):  TOB  = P_LEVK * WOB
P_KCUT = 16       # regenerative: axial cutting force per unit summed blade depth, S*a (N/m)
P_BLADES = 17     # regenerative: number of blades
P_KROP = 18       # quasi-static: ROP = P_KROP * max(WOB - P_W0, 0) * |omega|  (m/(N rad))
P_W0 = 19         # quasi-static: contact (threshold) WOB (N)
P_QS_TOP = 20      # static compression at the top face (N)
P_AD_KP = 21       # autodriller gains (m/s per N)
P_AD_KI = 22
P_AD_VMAX = 23     # autodriller max block speed (m/s)
P_HKLD_REF = 24    # off-bottom reference hookload for autodriller (N)
P_AXC = 25         # axial static contrast c: mu_s,axial = mu_k + c (mu_s - mu_k)
P_WEAK = 26        # 0: step static->kinetic, 1: exponential weakening mu_k + (mu_s - mu_k) exp(-s)
P_BLK_ACC = 27     # block acceleration limit (m/s^2), 0 = unlimited
P_HP_T = 28        # high-pass filter coefficient for the torsional Stokes damping (per step)
P_HP_A = 29        # same, axial
P_KV = 30          # 1 if Kelvin-Voigt material damping is active
N_PF = 31

# indices into the integer parameter vector ``pi``
I_CTRL = 0         # 0 internal PI/ZTorque, 1 external motor torque, 2 prescribed speed
I_AXIAL = 1        # 0 prescribed block velocity, 1 autodriller (WOB)
I_COUPLE = 2       # 1: 2-D Coulomb friction on (v, omega r) with live-tension F_N (default)
                   # 0: independent axial/torsional friction with static F_N
I_BIT = 3          # 0 Coulomb, 1 cutting (quasi-static), 2 regenerative cutting
N_PI = 4
BIT_MODELS = ('coulomb', 'cutting', 'regenerative')

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
    Quasi-static bit boundary (Coulomb and cutting models).

    Torque:  sticking (|omega_b| <= omega_c):  |TOB| <= LEVS * WOB
             sliding:                          TOB  =  LEVK * WOB * sign(omega_b)
    Penetration (depth of cut per revolution = 2 pi v / omega):
             ROP = KROP * max(WOB - W0, 0) * |omega_b|
    with the characteristic relations omega_b = alpha_t - TOB/Zt and
    v_b = alpha_a - WOB/Za, and unilateral contact (WOB >= 0).

    Returns (omega_b, tob, v_b, wob, sticking).
    """
    if (not on_bottom) or alpha_a <= 0.0:
        # off bottom, or the string is lifting the bit
        return alpha_t, 0.0, alpha_a, 0.0, 0.0

    lev_s = pf[P_LEVS]
    lev_k = pf[P_LEVK]
    kr = pf[P_KROP]
    w0 = pf[P_W0]
    wc = pf[P_WC]

    w = alpha_t if alpha_t > 0.0 else 0.0
    wob = 0.0
    v = 0.0
    tob = 0.0
    stuck = 0.0
    for _ in range(8):
        # axial: penetration law and v_b = alpha_a - WOB/Za
        q_lin = Za * alpha_a                     # WOB if the bit does not advance
        aw = abs(w)
        if q_lin <= w0 or kr * aw == 0.0:
            wob = q_lin
            v = 0.0
        else:
            wob = (alpha_a + kr * w0 * aw) / (1.0 / Za + kr * aw)
            v = kr * (wob - w0) * aw
        # torsional: TOB proportional to WOB, with sticking
        tau_stop = Zt * alpha_t                  # TOB that holds the bit at rest
        if (abs(w_prev) <= wc or stuck_prev > 0.5) and abs(tau_stop) <= lev_s * wob:
            tob = tau_stop
            w_new = 0.0
            stuck = 1.0
        else:
            tob = _sign(alpha_t) * min(lev_k * wob, abs(tau_stop))
            w_new = alpha_t - tob / Zt
            stuck = 0.0
        if abs(w_new - w) < 1e-12:
            w = w_new
            break
        w = w_new
    return w, tob, v, wob, stuck


@njit(cache=True)
def _bit_regen(alpha_t, alpha_a, Zt, Za, u, dt, w_prev, stuck_prev, phi, surf, pf):
    """
    Regenerative cutting bit (Richard, Germay & Detournay 2007 type).

    The hole bottom under the bit is stored as a profile ``surf`` over M
    angular bins.  ``n`` blades sit at phi + 2 pi j / n.  The depth of cut of
    blade j is d_j = max(u - surf[bin_j], 0): the material left by the
    previous blade to pass that angle, so the depth of cut carries the delay
    of one blade passage (the regenerative effect).

        WOB = S a sum_j d_j        (axial cutting force; zero when lifted off)
        TOB = WOB a / (2 zeta)     (cutting torque; sticking when stalled)

    The axial contact is solved implicitly for the new bit depth.
    Returns (omega_b, tob, v_b, wob, sticking).
    """
    kc = pf[P_KCUT]
    n = int(pf[P_BLADES])
    lev = pf[P_LEVK]
    wc = pf[P_WC]
    M = surf.shape[0]
    dth = 2.0 * np.pi / M

    # surface under each blade
    s_bl = np.empty(n)
    for j in range(n):
        ang = phi + 2.0 * np.pi * j / n
        b = int(np.floor(ang / dth)) % M
        s_bl[j] = surf[b]

    # implicit axial contact: u1 = u + dt (alpha_a - WOB/Za), WOB = kc sum (u1 - s_j)+
    u_free = u + dt * alpha_a
    u1 = u_free
    for _ in range(n + 2):
        m = 0
        ssum = 0.0
        for j in range(n):
            if s_bl[j] < u1:
                m += 1
                ssum += s_bl[j]
        if m == 0:
            u_new = u_free
        else:
            g = dt * kc / Za
            u_new = (u_free + g * ssum) / (1.0 + g * m)
        if abs(u_new - u1) < 1e-15:
            u1 = u_new
            break
        u1 = u_new
    wob = 0.0
    for j in range(n):
        if u1 > s_bl[j]:
            wob += kc * (u1 - s_bl[j])
    v = (u1 - u) / dt

    # torque: cutting torque proportional to WOB, sticking when stalled
    tau_stop = Zt * alpha_t
    if (abs(w_prev) <= wc or stuck_prev > 0.5) and abs(tau_stop) <= lev * wob:
        return 0.0, tau_stop, v, wob, 1.0
    tob = _sign(alpha_t) * min(lev * wob, abs(tau_stop))
    return alpha_t - tob / Zt, tob, v, wob, 0.0


@njit(cache=True)
def _sweep(phi_old, phi_new, u, surf, n):
    """Blades moving forward from phi_old to phi_new leave the surface at depth u."""
    if phi_new <= phi_old:
        return
    M = surf.shape[0]
    dth = 2.0 * np.pi / M
    for j in range(n):
        off = 2.0 * np.pi * j / n
        b0 = int(np.floor((phi_old + off) / dth))
        b1 = int(np.floor((phi_new + off) / dth))
        # bins the blade has left behind this step (it cuts the bin it occupies when it leaves)
        last = min(b1, b0 + M)
        for b in range(b0, last):
            k = b % M
            if u > surf[k]:
                surf[k] = u


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


@njit(cache=True)
def _kv_solve(x, m, h, k, dt, cp, dp):
    """
    Backward-Euler Kelvin-Voigt step, in place (Thomas algorithm):

        m_j h_j (x_j' - x_j) / dt = k_{j+1} (x_{j+1}' - x_j') - k_j (x_j' - x_{j-1}')

    with face conductances k (length n+1, k_0 = k_n = 0: no viscous flux
    through the top drive / bit faces).  Unconditionally stable and
    momentum-conserving; dissipates sum k (dx)^2 >= 0.
    """
    n = x.shape[0]
    for j in range(n):
        mh = m[j] * h[j] / dt
        a = -k[j]                    # sub-diagonal
        b = mh + k[j] + k[j + 1]     # diagonal
        c = -k[j + 1]                # super-diagonal
        d = mh * x[j]
        if j == 0:
            cp[j] = c / b
            dp[j] = d / b
        else:
            den = b - a * cp[j - 1]
            cp[j] = c / den
            dp[j] = (d - a * dp[j - 1]) / den
    x[n - 1] = dp[n - 1]
    for j in range(n - 2, -1, -1):
        x[j] = dp[j] - cp[j] * x[j + 1]


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
            wf, tf, vf, qf,
            # regenerative bit: hole-bottom profile over angular bins
            surf,
            # damping (per cell / per face), see docs/derivation sec. 5.6
            k0t, k0a, k1t, k1a, wbar, vbar, kvt, kva, fy):
    """
    Advance the model ``nsteps`` time steps.  State arrays are updated in
    place.

    Damping: k0* (1/s) act on the full velocity (Couette + calibration),
    k1* (1/s) on the oscillatory part w - wbar (Stokes layer; wbar is a
    per-cell low-pass state), kv* (N s/m resp. N m s) are Kelvin-Voigt face
    conductances (n+1, zero at the ends), fy (N/m) the yield-stress drag.  Every ``rec_every`` steps a row is written into ``rec`` starting at
    ``rec_idx``; returns (new rec_idx, new rec_phase).
    """
    n = h.shape[0]
    dt = pf[P_DT]
    a_hp_t = pf[P_HP_T]
    a_hp_a = pf[P_HP_A]
    use_kv = pf[P_KV] > 0.5
    w_old_a = np.empty(n)
    v_old_a = np.empty(n)
    wst = np.empty(n)
    vst = np.empty(n)
    cpr = np.empty(n)
    dpr = np.empty(n)
    wc = pf[P_WC]
    vc = pf[P_VC]
    I = pf[P_ITD]
    ctrl = pi[I_CTRL]
    axial_mode = pi[I_AXIAL]
    coupled = pi[I_COUPLE] == 1
    bit_model = pi[I_BIT]
    n_blades = int(pf[P_BLADES])
    axc = pf[P_AXC]
    weak_exp = pf[P_WEAK] > 0.5

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
            acc = pf[P_BLK_ACC]
            dv = v_blk_in - sc[S_VBLK]
            if acc > 0.0 and abs(dv) > acc * dt:
                sc[S_VBLK] += acc * dt * _sign(dv)
            else:
                sc[S_VBLK] = v_blk_in
        beta_a1 = v[0] - qd[0] / Za[0]
        q0 = Za[0] * (sc[S_VBLK] - beta_a1)
        v_top = sc[S_VBLK]
        # the hook can only hold the string up: hookload >= 0 (the block stacks out)
        if -(q0 + pf[P_QS_TOP]) < 0.0:
            q0 = -pf[P_QS_TOP]
            v_top = beta_a1 + q0 / Za[0]
        vf[0] = v_top
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
        if bit_model == 2:
            wb_, tob, vb_, wob, stuck = _bit_regen(alpha_t, alpha_a, Zt[n - 1], Za[n - 1], sc[S_BITD], dt,
                                                   sc[S_WBIT], sc[S_BITSTUCK], sc[S_PHIB], surf, pf)
        else:
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
            w_old_a[j] = w[j]
            v_old_a[j] = v[j]
            # flux update
            tau[j] -= dt * Kt[j] / h[j] * (wf[j + 1] - wf[j])
            ws = w[j] - dt / (mt[j] * h[j]) * (tf[j + 1] - tf[j])
            qd[j] -= dt * Ka[j] / h[j] * (vf[j + 1] - vf[j])
            vs = v[j] - dt / (ma[j] * h[j]) * (qf[j + 1] - qf[j])
            # implicit linear fluid damping: Couette (+ calibration) on the full
            # velocity, Stokes layer on the part above the high-pass corner
            wst[j] = (ws + dt * k1t[j] * wbar[j]) / (1.0 + dt * (k0t[j] + k1t[j]))
            vst[j] = (vs + dt * k1a[j] * vbar[j]) / (1.0 + dt * (k0a[j] + k1a[j]))

        # implicit Kelvin-Voigt material damping: one tridiagonal solve each
        if use_kv:
            _kv_solve(wst, mt, h, kvt, dt, cpr, dpr)
            _kv_solve(vst, ma, h, kva, dt, cpr, dpr)

        for j in range(n):
            w_old = w_old_a[j]
            v_old = v_old_a[j]
            ws = wst[j]
            vs = vst[j]
            fyj = fy[j]
            r = rc[j]
            if coupled:
                # ---- 2-D Coulomb inclusion on the sliding velocity (v, omega r)
                #      with the normal force from the live tension: axial load
                #      changes torque and rotation changes drag
                T = -(qs_c[j] + qd[j])
                fn = normal_force(T, wb_len[j], kx[j], ky[j], kz[j], zx[j], zy[j], zz[j])
                # static limit: an ellipse with axial and torsional static factors
                mu_sa = muk[j] + axc * (mus[j] - muk[j])
                mu_st = mus[j]
                vt_s = ws * r
                spd = np.sqrt(vs * vs + vt_s * vt_s)
                ea = vs / spd if spd > 0.0 else 0.0
                et = vt_s / spd if spd > 0.0 else 1.0
                s_old = np.sqrt((v_old / vc) ** 2 + (w_old / wc) ** 2)
                # resisting force capacities (N/m): wall friction + mud yield stress
                fa_s = mu_sa * fn + fyj
                ft_s = mu_st * fn + fyj
                if spd > 0.0 and fa_s > 0.0 and ft_s > 0.0:
                    f_sd = 1.0 / np.sqrt((ea / fa_s) ** 2 + (et / ft_s) ** 2)
                else:
                    f_sd = ft_s
                if s_old <= 1.0:
                    # near rest: stick if the static ellipse can hold the cell
                    pa = ma[j] * vs / (dt * fa_s) if fa_s > 0.0 else 1e30
                    pt = mt[j] * ws / (r * dt * ft_s) if ft_s > 0.0 else 1e30
                    if pa * pa + pt * pt <= 1.0:
                        v[j] = 0.0
                        w[j] = 0.0
                        continue
                    f_res = f_sd
                elif weak_exp:
                    f_res = muk[j] * fn + fyj + (f_sd - muk[j] * fn - fyj) * np.exp(-s_old)
                else:
                    f_res = muk[j] * fn + fyj
                va, vt = _coulomb_2d(vs, vt_s, ma[j], mt[j] / (r * r), dt * f_res)
                v[j] = va
                w[j] = vt / r
            else:
                # ---- independent 1-D Coulomb inclusions (Aarsnes & Shor eqs. 4, 20-22)
                fn = fn0[j]
                # torsion: torque per length needed to stop the cell this step
                t_stop = mt[j] * ws / dt
                if abs(w_old) <= wc:
                    cap = (mus[j] * fn + fyj) * r
                elif weak_exp:
                    cap = ((muk[j] + (mus[j] - muk[j]) * np.exp(-abs(w_old) / wc)) * fn + fyj) * r
                else:
                    cap = (muk[j] * fn + fyj) * r
                if t_stop > cap:
                    t_stop = cap
                elif t_stop < -cap:
                    t_stop = -cap
                w[j] = ws - dt * t_stop / mt[j]
                # axial: force per length needed to stop the cell this step
                f_stop = ma[j] * vs / dt
                mu_sa = muk[j] + axc * (mus[j] - muk[j])
                if abs(v_old) <= vc:
                    cap = mu_sa * fn + fyj
                elif weak_exp:
                    cap = (muk[j] + (mu_sa - muk[j]) * np.exp(-abs(v_old) / vc)) * fn + fyj
                else:
                    cap = muk[j] * fn + fyj
                if f_stop > cap:
                    f_stop = cap
                elif f_stop < -cap:
                    f_stop = -cap
                v[j] = vs - dt * f_stop / ma[j]

        # mean (low-passed) velocities seen by the Stokes-layer damping
        for j in range(n):
            wbar[j] += a_hp_t * (w[j] - wbar[j])
            vbar[j] += a_hp_a * (v[j] - vbar[j])

        # -------------------------------------------------- bookkeeping
        sc[S_T] += dt
        sc[S_BITD] += sc[S_VBIT] * dt
        sc[S_BLKPOS] += sc[S_VBLK] * dt
        sc[S_WINT] += abs(sc[S_WBIT]) * dt
        if bit_model == 2:
            phi_new = sc[S_PHIB] + sc[S_WBIT] * dt
            _sweep(sc[S_PHIB], phi_new, sc[S_BITD], surf, n_blades)
            # keep the angle bounded without disturbing bin alignment
            sc[S_PHIB] = phi_new - 2.0 * np.pi * np.floor(phi_new / (2.0 * np.pi)) if phi_new > 1e6 else phi_new
            if sc[S_BITD] > sc[S_HOLED]:
                sc[S_HOLED] = sc[S_BITD]
        elif sc[S_BITD] > sc[S_HOLED]:
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
