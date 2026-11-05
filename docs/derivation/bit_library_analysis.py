"""
Re-express each Bit Model Library model as forces vs depth of cut at the
sample-data parameters, and compute the sensitivity of each library ROP output
to each input.  Reproduces the numbers in bit-models.md.

Usage:
    git clone https://github.com/Kartikay100/Bit-Model-Library
    python bit_library_analysis.py path/to/Bit-Model-Library
"""
import math, os, sys
import numpy as np, pandas as pd
LIB = sys.argv[1] if len(sys.argv) > 1 else 'Bit-Model-Library'
d0 = pd.read_excel(os.path.join(LIB, 'SampleData', 'Input_File.xlsx'))
d0.columns = [str(c).strip() for c in d0.columns]
P = d0.iloc[1].to_dict()                                   # middle sample row
a = P['bit_radius (inches)']; eps = P['rock_intrinsic_energy (lbf/in²)']; zeta = P['zeta (unitless)']
sig = P['max_contact_pressure (lbf/in²)']; l = P['length_of_wear_flat (inches)']; mu = P['mu (unitless)']
alpha = math.radians(P['alpha (degrees)']); beta = math.radians(P['beta (degrees)']); dc = P['dc (inches)']
nc = P['nc (unitless)']; Rc = P['Rc (lbf/in²)']; Rp = P['Rp (lbf/in²)']; sc_ = P['sigma_c (lbf/in²)']; st = P['sigma_t (lbf/in²)']
w = P['w (inches)']; thf = math.radians(P['theta_f (degrees)']); s0 = P['sigma_0 (lbf/in²)']
print('sample row:', {k: round(v, 4) if isinstance(v, float) else v for k, v in P.items() if k in ('rotational_speed (rpm)','wob (lbf)','tob (lbf·ft)','rop (ft/hr)','d (inches)')})
W, T = P['wob (lbf)'], P['tob (lbf·ft)']
print(f"implied bit friction 2T/(aW): TOB as lbf·ft -> {2*T*12/(a*W):.2f};  TOB as lbf·in -> {2*T/(a*W):.2f}")

d = np.linspace(0.0005, 0.02, 200)       # depth of cut per revolution (in)

# ---- Detournay & Defourny (1992), textbook form
Wf_dd = a * l * sig;  Wc = zeta * eps * a * d
Tc = eps * a**2 * d / 2;  Tf = mu * 1.0 * a * Wf_dd / 2
print(f"\nD&D: W_f = {Wf_dd:.0f} lbf ({100*Wf_dd/W:.0f}% of sample WOB); T_f = {Tf/12:.0f} lbf·ft; "
      f"cutting slope dW/dd = {zeta*eps*a:.0f} lbf/in, dT/dW (cutting) = a/2ζ = {a/(2*zeta):.2f} in")

# ---- Ziaja (1985): per-cutter chord area Sc(d) and wear area Sw
def Sc(dd):
    ca = math.cos(alpha)
    return math.cos(beta) * (0.24*math.pi*dc**2 - (dc/(6*ca))*(4*np.sqrt(2*dc*ca*(dc*ca-dd)) - np.sqrt(16*dc**2*ca**2-dd**2))
                              - (dd/(2*ca**2))*np.sqrt(16*dc**2*ca**2-dd**2))
x = l*math.cos(alpha)*math.sin(alpha)
Sw = (4/3)*x*(math.tan(alpha)+1/math.tan(alpha))*math.sqrt(dc*x*math.cos(alpha)-x**2)/math.cos(alpha)
sc = Sc(d)
lo = d < 0.01
n_exp = np.polyfit(np.log(d[lo]), np.log(np.abs(sc[lo])), 1)[0]
print(f"\nZiaja: Sc(d) over d = 0.0005-0.01 in fits Sc ∝ d^{n_exp:.2f};  Sc(0) = {Sc(1e-9):.2e} in² (should be 0)")
print(f"       normal/tangential cutting ratio = tan(alpha) = {math.tan(alpha):.3f};  wear area Sw = {Sw:.2e} in²")
print(f"       wear-flat force per cutter Rp·Sw = {Rp*Sw:.2f} lbf  (x{nc:.0f} cutters = {nc*Rp*Sw:.0f} lbf; library multiplies by nc twice = {nc*nc*Rp*Sw:.0f} lbf)")

# ---- Che et al. (2016): per-cutter normal force linear in d
t1 = -math.cos(alpha)/(4*math.sin(math.pi/4-alpha/2)**2 + 4*st/sc_)
t2 = math.sqrt(math.cos(alpha)**2 + 4*(st/sc_)*math.sin(math.pi/4-alpha/2)**2 + 4*math.sin(math.pi/4-alpha/2)**4)
t3 = math.sqrt(st/sc_+1)/(2*math.sin(math.pi/4-alpha/2))
t4 = 0.5*math.cos(math.pi/4-alpha/2)/math.sin(math.pi/4-alpha/2)
g = math.cos(alpha)*(t1 + t2/(4*math.sin(math.pi/4-alpha/2)**2 + 4*st/sc_) + t3 + t4)
print(f"\nChe: per-cutter F_n = σc·w·g(α, σt/σc)·d with g = {g:.3f}  ->  normal stress on cut = {g*sc_:.0f} psi "
      f"(vs D&D ζε = {zeta*eps:.0f} psi); linear in d, no torque output")

# ---- Gerbaud et al. (2006): F_n/A = tan(θf+α)·ε  -> ζ_eff
print(f"\nGerbaud: normal/tangential = tan(θf+α) = {math.tan(thf+alpha):.2f} (an effective ζ); "
      f"wear-flat stress σ0 = {s0:.0f} psi; chamfer term zero in sample (A_ch = 0)")

# ---- sensitivity of library ROP outputs (as coded) to each input, sample row
import importlib.util, sys
spec = importlib.util.spec_from_file_location('m', os.path.join(LIB, 'BitModel', 'models.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
base = {k: P[k] for k in P}
def rops(p):
    M = m.PDC_BIT_MODELS(p)
    return dict(DD=M.detournay_rop()[0], Ziaja=M.ziaja_rop()[2], Che=M.che_model()[2], Gerbaud=M.gerbaud_normal_force_rop()[2])
r0 = rops(base)
print('\nlibrary ROP (ft/h) at sample row:', {k: round(v, 1) for k, v in r0.items()}, ' measured:', P['rop (ft/hr)'])
rows = []
for k in ['rock_intrinsic_energy (lbf/in²)', 'zeta (unitless)', 'max_contact_pressure (lbf/in²)', 'length_of_wear_flat (inches)',
          'mu (unitless)', 'alpha (degrees)', 'beta (degrees)', 'dc (inches)', 'nc (unitless)', 'Rc (lbf/in²)', 'Rp (lbf/in²)',
          'd (inches)', 'sigma_c (lbf/in²)', 'sigma_t (lbf/in²)', 'w (inches)', 'A (inches²)', 'sigma_0 (lbf/in²)', 'theta_f (degrees)',
          'wob (lbf)', 'rotational_speed (rpm)']:
    p = dict(base); p[k] = base[k] * 1.1
    r = rops(p)
    rows.append([k] + [ (r[mm]/r0[mm]-1)/0.1 if r0[mm] else np.nan for mm in r0])
S = pd.DataFrame(rows, columns=['parameter (+10%)'] + list(r0)).set_index('parameter (+10%)').round(2)
print('\nelasticity of each library ROP to each input (dROP/ROP per dX/X):')
print(S.replace(0.0, '·').to_string())
