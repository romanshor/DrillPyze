# Rock strength, pore pressure and temperature at the bit

This note derives how DrillPyze sets the rock strength seen by the bit
(`DrillPyze.formation`) and how it models the wellbore temperature that feeds
it (`DrillPyze.thermal`). The bit models of [`bit-models.md`](bit-models.md) all
scale with the drilling strength $`S = \zeta\varepsilon`$. When
`BitParameters.strength_source = 'formation'`, the intrinsic specific energy is
taken from the rock's confined compressive strength:

```math
\varepsilon = \frac{\text{CCS}_\text{eff}}{\eta},\qquad S = \frac{\zeta\,\text{CCS}_\text{eff}}{\eta}.
```

$`\eta`$ is a cutting efficiency, equal to 1 for a sharp cutter. Two lines of
evidence support linking $`\varepsilon`$ to compressive strength:
* scratch tests show that $`\varepsilon`$ is well correlated with UCS for a sharp
  cutter in ductile mode [Mitaim et al. 2004];
* field MSE methods equate bit specific energy with CCS divided by a mechanical
  efficiency [Teale 1965; Pessier & Fear 1992; Caicedo et al. 2005].

$`\text{CCS}_\text{eff}`$ combines three effects, described in §1–3. All depths
are in m, pressures in Pa and temperatures in °C.

---

## 1. Confinement by the overbalance (P1)

Mohr–Coulomb strength under a confining stress equal to the differential
pressure $`\Delta p`$ across the rock at the bit gives
[Caicedo, Calhoun & Ewy 2005, eq. 2]

```math
\text{CCS} = \text{UCS} + \Delta p\,\frac{1+\sin\varphi}{1-\sin\varphi}
= \text{UCS} + \Delta p + \frac{2\,\Delta p\,\sin\varphi}{1-\sin\varphi}.
```

### Permeable rock (drained)

Pore pressure has time to equalise at the cutting face:

```math
\Delta p_\text{HP} = p_\text{bh} - p_p .
```

### Low-permeability, low-porosity rock (undrained)

When the bit removes the overburden from the rock beneath it, the rock is
unloaded from $`\sigma_v`$ to $`p_\text{bh}`$ faster than pore fluid can flow. The
pore pressure then drops by Skempton's response to the change in mean stress
[Skempton 1954]. Following Caicedo et al. (eqs. 11–13), the effective pore
pressure is

```math
p_p^{*} = p_p - B\,\frac{\sigma_v - p_\text{bh}}{3},
\qquad \Delta p_\text{LP} = p_\text{bh} - p_p^{*},
```

with $`B = 1`$ in their formulation. That is the default here, and it can be
adjusted when lab $`B`$ values are available. For FORGE, OU measured Skempton
and Biot coefficients on 58-32 and 78-32 cores [GDR submission 1464].

### Blending by porosity

Effective porosity $`\phi_e`$ is the proxy for drainage, and the strength is
blended linearly [Caicedo et al. 2005]:

| $`\phi_e`$ | CCS used |
|---|---|
| $`\le 0.05`$ | $`\text{CCS}_\text{LP}`$ |
| $`\ge 0.20`$ | $`\text{CCS}_\text{HP}`$ |
| in between | linear blend of the two |

Granite has very low permeability (FORGE matrix permeability below 30 µD) and
very low porosity (matrix porosity below 0.5%) [Utah FORGE review 2026]. It is
therefore treated as fully undrained. For FORGE-like conditions this raises CCS
well above the drained value. For example, at 2,500 m TVD with hydrostatic pore
pressure and a 25 kPa/m overburden gradient, $`\Delta p_\text{LP}`$ is about
12 MPa larger than $`\Delta p_\text{HP}`$.

**Context.** Detournay & Atkinson (2000) showed that in shear-dilatant,
low-permeability rock, cutting at drilling speeds is undrained. Detournay & Tan
(2002) found that shale specific energy then depends on $`p_\text{bh}`$ and not on
$`p_\text{bh} - p_p`$: $`\varepsilon = \varepsilon_0 + m\,p_\text{bh}`$. The
Skempton construction above is the engineering counterpart. Both say the far
field pore pressure matters little in tight rock. Whether FORGE granite follows
the drained or the undrained limit can be tested directly against the field
ROP data.

---

## 2. Thermal weakening (T1)

```math
\text{UCS}(T) = \text{UCS}_0\,\bigl[1 - \beta\,\max(T - T^{*}, 0)\bigr],
```

where $`T`$ is the undisturbed formation temperature at the bit, from the
geothermal gradient. The rock ahead of the bit has not yet been cooled
[Wang et al. 2023, assumption 6].

$`T^{*} = 193`$ °C and $`\beta = 1.41\times10^{-3}`$ /K are a least-squares fit
(RMS error 3.8%) to the residual UCS of granite heated to 250–650 °C and air
cooled [Xu et al. 2019]:

| T (°C) | 250 | 350 | 450 | 550 | 650 |
|---|---|---|---|---|---|
| measured loss | 12% | 14% | 38% | 53% | 63% |
| fitted loss | 8% | 22% | 36% | 50% | 64% |

At FORGE temperatures (about 200–230 °C) the effect is small, 1–5%.

Caveat: these are residual strengths measured after cooling. Strength measured
at temperature can differ, and so can the effect of repeated quenching.

---

## 3. Cooling-induced thermal stress at the hole bottom (T2)

Circulating mud is cooler than the rock at the bit by
$`\Delta T = T_\text{rock} - T_\text{mud,bit}`$. A cooled surface layer that is
constrained laterally carries an in-plane thermoelastic tensile stress
[e.g. Timoshenko & Goodier 1970]:

```math
\sigma_\text{th} = \frac{\alpha E\,\Delta T}{1-\nu}.
```

For FORGE granite ($`\alpha = 6\times10^{-6}`$ /K, $`E = 55`$ GPa, $`\nu = 0.26`$) this
is about 0.45 MPa/K, or 45 MPa per 100 K. That is comparable to the tensile
strength (about 11 MPa) and to the overbalance term.

The cooled layer is thin. Between blade passes ($`t_p = 2\pi/(n\omega)`$, with $`n`$
blades) heat diffuses to a depth

```math
\delta = \sqrt{\pi\,D_\text{th}\,t_p},
\qquad D_\text{th} = \frac{k}{\rho c}\approx 0.9\ \text{mm}^2/\text{s},
```

while each blade cuts to depth $`d_b = \text{ROP}\cdot t_p`$. Only the cooled
fraction of the cut is weakened:

```math
\text{CCS}_\text{eff} = \max\!\Bigl(\text{CCS} - \kappa\,\sigma_\text{th}\,\min(1, \delta/d_b),\ \ f_\text{floor}\,\text{UCS}(T)\Bigr).
```

* $`\kappa \in [0, 1]`$ is the efficiency of the thermal stress. It is the one
  parameter that has to be calibrated; the default is 0.5.
* $`f_\text{floor} = 0.1`$ keeps the strength physical.
* The term makes the thermal benefit depend on ROP, RPM and blade count. It is
  largest at low ROP, high RPM and with many blades.

Status: this formulation is DrillPyze's own, built from standard
thermoelasticity and conduction. SPE-217718 studies the effect of rapid cooling
on cutting hot dry rock experimentally and is the natural source for
calibrating $`\kappa`$.

---

## 4. Wellbore temperature: the growing-wellbore drift-flux model

The mud temperature at the bit and the bottomhole pressure come from the model
of **Wang, Khaled, Luu, Ashok & van Oort (2023), SPE-214836-MS**. That model
couples thermohydraulic flow to a wellbore that grows as it is drilled, and it
was validated against Utah FORGE 16(78)-32. It builds on Evje & Fjelde (2002),
Fallah et al. (2019) and Khaled et al. (2023c).

DrillPyze implements it in the **single-phase liquid** form (gas volume fraction
zero), which is the case for water-based mud.

### Conservation equations

For each flow path (inside the drillstring, flowing down; in the annulus,
flowing up), the mixture energy equation (paper eq. 4) becomes, per cell of
length $`\Delta x`$:

```math
\rho c A\,\Delta x\,\frac{dT}{dt}
= \rho c Q\,(T_\text{upstream} - T) + \dot Q_\text{wall}\,\Delta x
+ \dot m_\text{cut} c_r\,(T_\text{rock} - T)\,\delta_{k,\text{bottom}} .
```

* $`Q`$ is the pump rate.
* The mass and momentum equations (eqs. 1, 3) reduce to $`v = Q/A`$. Assumption 1
  of the paper, ROP ≪ mud velocity, justifies this.
* Mud density follows paper eq. 2:

  ```math
  \rho(T,p) = \bigl(1 - (T - T_\text{ref})\,c_\text{exp}\bigr)\rho_\text{ref} + \frac{p - p_\text{ref}}{v_\text{sound}^2}.
  ```

  It is integrated up the annulus, together with the frictional loss (Fanning,
  Newtonian), to give $`p_\text{bh}`$ and the ECD used in §1.

### Radial heat transfer (paper eqs. 6–9)

```math
\dot q_\text{ds} = \frac{T_\text{ann} - T_\text{ds}}{R_\text{ds}},\qquad
\dot q_\text{ann} = \frac{T_\text{for} - T_\text{ann}}{R_\text{for}} + \frac{T_\text{ds} - T_\text{ann}}{R_\text{ds}} .
```

* $`R_\text{ds}`$ is inside convection, plus conduction through the pipe wall,
  plus outside convection.
* $`R_\text{for}`$ is wall convection, plus conduction through the casing where
  cased, plus conduction to the first formation node.
* Each convective resistance is $`1/(h\pi D)`$, with $`h = \text{Nu}\,k/D_h`$.
* The paper does not state its Nusselt correlation. DrillPyze uses
  Dittus–Boelter, $`\text{Nu} = 0.023\,\text{Re}^{0.8}\text{Pr}^{0.4}`$ [Dittus &
  Boelter 1930], for turbulent flow and $`\text{Nu} = 3.66`$ for laminar flow.

### Formation

Each axial cell has 12 radial conduction cells, log-spaced out to 5 m beyond
the wall, solved implicitly. The outer boundary is held at the undisturbed
geothermal temperature (paper boundary conditions; Khaled et al. 2023c).

### Growing wellbore (paper eqs. 10, 14–20; Figs. 1–3)

* **Bottom cell.** The bottom cell lengthens at the ROP. When it reaches 1.5 cell
  lengths it is cut in two, and both halves keep the state of the original cell
  ("grow-cut-grow").
* **Cuttings.** Cuttings enter the bottom annulus cell at the undisturbed
  formation temperature. DrillPyze uses the rock's density and heat capacity for
  them; the paper assumes mud properties.
* **Bottom formation cells.** As the hole deepens, they relax toward the
  undisturbed temperature:

```math
T^\text{bottom}_{t+\Delta t} = T^\text{bottom}_t + \frac{\text{ROP}\,\Delta t}{L_\text{bottom}}\bigl(T^\text{bottom,formation}_t - T^\text{bottom}_t\bigr).
```

* **Moving string.** The string geometry is volume-averaged over each cell at
  every update (Fig. 1). The ROP-scaled area-change terms of eqs. 11–13 are of
  order ROP/$`v`$ and are omitted.

### Numerics

Each step uses backward Euler for the fluid and the formation, coupled by three
Gauss–Seidel iterations: drillstring downstream, annulus upstream from the bit,
then a tridiagonal solve for each formation layer. The time step is at most
2 s, and the model is updated from the drilling simulation every 2 s of
simulated time.

### Verification against SPE-214836-MS

`examples/forge_thermal_spe214836.py` and `tests/test_thermal_rock.py` use the
paper's FORGE case: about 5,900 ft MD, 550 gpm, 40 °C surface temperature and
69.1 °C/km, with Table 2 properties. The inlet temperature is not given in the
paper; 35 °C is assumed.

| quantity | SPE-214836 | DrillPyze |
|---|---|---|
| downhole temperature at pump start → after 30 min | ≈ 310 → 210 °F (Fig. 6) | 325 → 209 °F |
| downhole rise from drilling after 30 min, 40 / 60 / 72 ft/h | ≈ 1.3 / 2.2 / 2.7 °F (Fig. 7) | 1.25 / 1.87 / 2.24 °F |
| rise linear in ROP | yes | yes |
| ROP → 0 converges to the static wellbore | yes | yes (test) |
| static well stays at the geothermal profile | – | exact (test) |
| energy balance while circulating | – | within 3% (test) |

**Open item.** The time history of the drilling-induced rise differs:
* In the paper it jumps to most of its final value within a minute and then
  plateaus.
* In DrillPyze it builds steadily. About 60% of it comes from the bit advancing
  into hotter rock along the 69 °C/km gradient; cuttings and the freshly
  exposed wall give a prompt rise of about 0.4 °F at 60 ft/h.
* Possible causes include the unstated inlet temperature, the paper's mesh, and
  details of its eqs. 10 and 19 (for example, whether the cuttings enthalpy is
  taken relative to the mud temperature or to a reference temperature).
* The result is insensitive to the DrillPyze cell size (2–30 m).

---

## 5. Parameters and FORGE defaults

| parameter | default | source |
|---|---|---|
| UCS | 110 MPa | FORGE core 91–126 MPa [Ye et al. 2022, via Utah FORGE review 2026] |
| friction angle φ | 36° | FORGE granitoid model value [Finnila & Podgorney 2020]; core 32–48° |
| porosity | 0.005 | matrix < 0.5% [Utah FORGE review 2026] |
| pore pressure gradient | 9.79 kPa/m | hydrostatic |
| overburden gradient | 25 kPa/m | assumed |
| Skempton B | 1.0 | Caicedo et al. 2005 |
| surface T, gradient | 40 °C, 69.1 °C/km | SPE-214836 Table 2 |
| rock ρ, c, k | 2800 kg/m³, 930 J/kg K, 2.31 W/m K | SPE-214836 Table 2 |
| α, E, ν | 6×10⁻⁶ /K, 55 GPa, 0.26 | Podgorney et al. 2021; Utah FORGE review 2026 |
| T*, β | 193 °C, 1.41×10⁻³ /K | fit to Xu et al. 2019 |
| κ | 0.5 | to be calibrated (§3) |
| mud c, k, expansion, viscosity | 3750 J/kg K, 0.75 W/m K, 10⁻⁴ /K, 20 mPa·s | SPE-214836 Table 2 |

---

## References

* **Wang, N., Khaled, M.S., Luu, A., Ashok, P. and van Oort, E., 2023.** Downhole Temperature Estimation of a Growing High-Temperature Wellbore Using a Modified Drift Flux Modeling Approach. SPE Annual Technical Conference and Exhibition, San Antonio, Texas. SPE-214836-MS. https://doi.org/10.2118/214836-MS
* Khaled, M.S., Wang, N., Ashok, P. and van Oort, E., 2023c. Downhole heat management for drilling shallow and ultra-deep high enthalpy geothermal wells. *Geothermics* 107, 102604. https://doi.org/10.1016/j.geothermics.2022.102604
* Evje, S. and Fjelde, K.K., 2002. Relaxation schemes for the calculation of two-phase flow in pipes. *Math. Comput. Model.* 36(4–5), 535–567.
* Fallah, A. et al., 2019. An integrated thermal and multi-phase flow model for estimating transient temperature dynamics during drilling operations. SPE/IADC-194083-MS.
* Caicedo, H.U., Calhoun, W.M. and Ewy, R.T., 2005. Unique ROP predictor using bit-specific coefficient of sliding friction and mechanical efficiency as a function of confined compressive strength impacts drilling performance. SPE/IADC-92576-MS. CCS equations as given in US Patent 7,412,331 (https://patents.google.com/patent/US7412331) and the Skempton low-permeability form in US Patent 7,555,414.
* Skempton, A.W., 1954. The pore-pressure coefficients A and B. *Géotechnique* 4(4), 143–147.
* Detournay, E. and Atkinson, C., 2000. Influence of pore pressure on the drilling response in low-permeability shear-dilatant rocks. *Int. J. Rock Mech. Min. Sci.* 37(7), 1091–1101. https://doi.org/10.1016/S1365-1609(00)00050-2
* Detournay, E. and Tan, C.P., 2002. Dependence of drilling specific energy on bottom-hole pressure in shales. SPE/ISRM Rock Mechanics Conference, Irving, Texas, 526–535.
* Mitaim, S., Dagrain, F., Richard, T., Detournay, E. and Drescher, A., 2004. A novel apparatus to determine the rock strength parameters. 9th National Convention on Civil Engineering, Thailand. https://orbi.umons.ac.be/handle/20.500.12907/1650
* Teale, R., 1965. The concept of specific energy in rock drilling. *Int. J. Rock Mech. Min. Sci.* 2, 57–73.
* Pessier, R.C. and Fear, M.J., 1992. Quantifying common drilling problems with mechanical specific energy and a bit-specific coefficient of sliding friction. SPE-24584-MS.
* Xu, C., Sun, Q., Pan, X., Zhang, W. and Wang, Y., 2019. Variation on thermal damage rate of granite specimen with thermal cycle treatment. *High Temp. Mater. Process.* 38, 849–855. https://doi.org/10.1515/htmp-2019-0046
* SPE-217718, 2025. Study examines thermal effects on drilling performance in hot dry rock (*JPT* synopsis). https://jpt.spe.org/study-examines-thermal-effects-on-drilling-performance-in-hot-dry-rock-restricted
* Timoshenko, S.P. and Goodier, J.N., 1970. *Theory of Elasticity*, 3rd ed. McGraw-Hill (thermal stresses).
* Dittus, F.W. and Boelter, L.M.K., 1930. Heat transfer in automobile radiators of the tubular type. *Univ. Calif. Publ. Eng.* 2, 443–461.
* Utah FORGE: A Decade of Innovation — Comprehensive Review of Field-Scale Advances (Part 1), 2026. *Processes* 14(3), 512. https://doi.org/10.3390/pr14030512 (compiles Ye et al. 2022, ARMA-2022-0272; Finnila & Podgorney 2020, SGP-TR-216; Podgorney et al. 2021, *Energies* 14, 4758).
* Utah FORGE 5-2615: Well 58-32 and 78-32 Poroelastic Properties (dataset), GDR submission 1464. https://gdr.openei.org/submissions/1464
