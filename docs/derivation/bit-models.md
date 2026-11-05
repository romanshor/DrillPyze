# Bit–rock interaction: from the Bit Model Library to three reduced models

This note analyses the four PDC bit models in the
[Bit Model Library](https://github.com/Kartikay100/Bit-Model-Library)
(SPE-223715-MS) and reduces them to three models, each with two or three
parameters, for use as the bit boundary condition in DrillPyze.

The models in the library are:
* Ziaja (1985);
* Detournay & Defourny (1992);
* Gerbaud et al. (2006);
* Che et al. (2016).

The script `bit_library_analysis.py` in this folder reproduces every number
below from the library's code and sample data.

Notation:
* $a = D/2$ is the bit radius;
* $\omega$ is the bit speed;
* $v$ is the axial bit velocity;
* $d = 2\pi v/\omega$ is the depth of cut per revolution;
* $W$ is weight on bit (WOB) and $T$ is torque on bit (TOB).

---

## 1. What each library model is

### Detournay & Defourny (1992)

The bit force is split into a cutting part and a frictional-contact part (wear
flats):

$$
W = \underbrace{\zeta\varepsilon a d}_{W_c} + \underbrace{a\ell\sigma}_{W_f},
\qquad
T = \underbrace{\tfrac12\varepsilon a^2 d}_{T_c} + \underbrace{\tfrac12\mu\gamma a W_f}_{T_f}.
$$

* **Inputs:** $\varepsilon$ (intrinsic specific energy), $\zeta$, $\ell$
  (wear-flat length), $\sigma$ (contact pressure), $\mu$, $\gamma$ (bit
  constant) and $a$. That is six rock and bit parameters plus geometry.
* **What the output depends on:** in the library's ROP calculation, the result
  depends only on the product $\zeta\varepsilon$ and on $W_f$. The
  elasticities are −0.91 to $\varepsilon$, −0.91 to $\zeta$ and −0.14 to each of
  $\sigma$ and $\ell$; $\mu$ has no effect.

### Ziaja (1985)

This is a cutter-level model.
* **Forces per cutter:** the cutting forces on each cutter are $R_c S_c$, where
  $S_c(d)$ is the chord area of a round cutter. They are split by the back rake
  angle into a normal part, $R_c S_c\sin\alpha$, and a tangential part,
  $R_c S_c\cos\alpha$. A wear-flat force $R_p S_w$ is added.
* **Inputs:** $R_c$, $R_p$, $\alpha$, $\beta$, $d_c$, $n_c$, wear length and
  $\mu$.
* **Normal-to-tangential ratio:** it is fixed at $\tan\alpha$, which is 0.21 in
  the sample.
* **ROP:** the library computes ROP as $\text{rpm}\cdot d\cdot n_c$, with $d$ as
  an *input*. So ROP has no dependence on WOB (elasticity 0) and is
  proportional to $d$ and $n_c$.

### Che et al. (2016)

This is a single-cutter chipping and crushing model.
* **Normal force:** $F_n = \sigma_c\,w\,g(\alpha, \sigma_t/\sigma_c)\,d$, which
  is **linear in depth of cut**. In the sample, $g = 1.57$.
* **Use in the library:** it is scaled to the bit with $n_c$, after subtracting
  the Detournay–Defourny wear-flat force.
* **What ROP depends on:** WOB (elasticity 1.14), $\sigma_c$ (−0.9), cutter
  width (−0.91) and $W_f$ (−0.14). $\sigma_t$ and $\alpha$ hardly matter
  (≤ 0.08).
* **Torque:** the model gives none.

### Gerbaud et al. (2006)

This is a cutter model with cutting, chamfer and wear-flat terms.
* **Cutting force:** the normal cutting force is
  $\tan(\theta_f + \alpha)\,A\,\varepsilon$. So Gerbaud's effective
  normal-to-tangential ratio is $\tan(\theta_f+\alpha)$, which is 1.80 in the
  sample.
* **ROP:** in the library ROP follows the cut area $A$ (elasticity 0.93) and
  $\theta_f$ (2.2), and has **no dependence on WOB**, because $A$ is supplied as
  an input.

---

## 2. What they have in common

With the bookkeeping removed, all four models share one structure:

$$
\boxed{\;W = W_f + S\,a\,d, \qquad T = \frac{a}{2\zeta}\,W_c + \mu_f\,\frac{a}{2}\,W_f\;}
$$

1. **Cutting forces are proportional to the depth of cut.** This is exact in
   Detournay & Defourny, Che and Gerbaud. Ziaja's chord area is nominally
   nonlinear, but see §5.
   * Every rock and cutter input ($\sigma_c$, $\sigma_t$, $R_c$, $\alpha$,
     $\beta$, $\theta_f$, $d_c$, $w$, $n_c$, $A$) enters only through two
     numbers. The first is the **drilling strength** $S = \zeta\varepsilon$, the
     axial cutting force per unit cut area. The second is the **force ratio**
     $\zeta$, normal over tangential cutting force.
   * Che and Gerbaud are best read as ways to *estimate* $S$ (and $\zeta$) from
     rock strength and cutter geometry, not as different laws.
2. **Wear-flat or chamfer contact adds a force that does not depend on depth of
   cut.** This is $W_f$, with its frictional torque. It appears as $a\ell\sigma$
   in Detournay & Defourny, $R_p S_w$ in Ziaja and $\sigma_0 A_f$ in Gerbaud,
   and Che subtracts it.
3. **All four are rate independent.** None of them has velocity weakening.
   Torque depends on $d$ and $W$, not on $\omega$. For drillstring dynamics this
   matters: at the bit, stick-slip then has to come from somewhere other than
   these laws. The options are a static/kinetic friction difference, the
   regenerative (delayed) depth of cut, or along-string friction.

Three things identify the bit from field data:
* the slope of TOB against WOB, which gives $a/(2\zeta)$;
* the slope of WOB against $d = 2\pi\,\text{ROP}/\omega$, which gives $Sa$;
* that line's intercept, which gives $W_f$.

These are the three numbers the reduced models keep.

---

## 3. Simplifications and why

| Simplification | Justification | Effect in the sample data |
|---|---|---|
| Collapse every rock and cutter input into $S$ and $\zeta$ | Every model is linear in $d$; the inputs only scale the slope (§2) | Exact |
| Contact friction $\mu_f = \mu\gamma \approx 1/\zeta$, so $T = aW/(2\zeta)$ over the whole range | Removes a parameter that cannot be told apart from $\zeta$ without data at very low WOB | Sample has $\mu\gamma = 0.9$ and $1/\zeta = 1.0$, so the torque error is $\tfrac a2 W_f(1/\zeta-\mu\gamma) \approx$ 1.7% of TOB |
| Linear rather than chord-area dependence on $d$ | Over 0.0005–0.01 in, the Ziaja area as implemented scales as $d^{-0.02}$, so it is not usable as a nonlinearity. For a correct chord area the nonlinearity matters only at depths of cut well below field values | – |
| Depth of cut per bit revolution, $d = 2\pi v/\omega$, not per cutter × $n_c$ | Avoids double counting cutters (§5) and needs no cutter layout | – |
| Keep a static/kinetic difference only in the friction-dominated model | Field evidence of bit-induced velocity weakening is for worn bits in the contact regime | – |

---

## 4. The three reduced models

All three are implemented in the kernel (`BitParameters.model`), selectable in
the UI, and share the drilling strength $S$.

### Model 1: Coulomb (friction-dominated)

Parameters: $\mu_s$, $\mu_k$, $S$.

$$
|T| \le \mu_s\,W\,\tfrac{D}{3}\ \ (\text{sticking}),\qquad
T = \mu_k\,W\,\tfrac{D}{3}\ \ (\text{sliding}),\qquad
\text{ROP} = \frac{\omega\,W}{2\pi a S}.
$$

* **Physics:** this is the Detournay & Defourny model in the contact limit
  ($W \approx W_f$, a worn or dull bit), where torque is friction,
  $T = \mu\gamma\tfrac a2 W$. To that it adds the static/kinetic difference
  that gives velocity weakening and bit-induced stick-slip.
* **Fitting:**
  * $\mu_k$ from TOB/WOB while drilling;
  * $\mu_s$ from the breakout torque on bottom;
  * $S$ from the slope of ROP/RPM against WOB.

### Model 2: Cutting (Detournay & Defourny, quasi-static)

Parameters: $S$, $\zeta$, $W_f$.

$$
W = W_f + S\,a\,d,\qquad T = \frac{a}{2\zeta}\,W,\qquad
\text{ROP} = \frac{\omega\,(W-W_f)^+}{2\pi a S}.
$$

* **Physics:** a sharp to moderately worn PDC bit with a cutting response.
  It has a WOB threshold before penetration and **no velocity weakening**. The
  bit sticks only when the string cannot supply $aW/(2\zeta)$.
* **Fitting:**
  * $\zeta$ from the slope of TOB against WOB;
  * $S$ and $W_f$ from a straight-line fit of WOB against
    $d = 2\pi\,\text{ROP}/\omega$.

### Model 3: Regenerative cutting (Richard, Germay & Detournay 2007)

Parameters: $S$, $\zeta$, plus the blade count $n$ (bit geometry).

$$
W = S\,a\sum_{j=1}^{n} d_j,\qquad T = \frac{a}{2\zeta}\,W,\qquad
d_j = \bigl(u_b - s(\theta_j)\bigr)^+ .
$$

* **Physics:** the same cutting law, but the depth of cut of each blade is
  measured against the **hole-bottom profile** $s(\theta)$ left by the previous
  blade. The profile is stored over 720 angular bins and updated as the blades
  sweep. This brings in the one-blade-passage delay (the regenerative effect),
  which couples axial and torsional motion at the bit.
* **Behaviour in the model:**
  * Lift-off and stalling come out naturally: if the bit rises above the
    profile, $W = 0$; if it stops, the profile stops being cut.
  * With enough damping it settles to exactly the Model 2 state with $W_f = 0$
    (test `test_regenerative_bit_converges_to_cutting_when_stable`).
  * With light damping it self-excites into coupled bit bounce and stick-slip
    (test `test_regenerative_bit_self_excites_when_lightly_damped`).
* **Fitting:** the same $S$ and $\zeta$ as Model 2. The steady states agree, so
  **no extra fitting** is needed. $n$ comes from the bit specification.

### Choosing between them

| Situation | Model |
|---|---|
| Worn or dull bit; torque dominated by friction; on-bottom stick-slip seen in the data | Coulomb |
| PDC bit at steady drilling; ROP and torque fitting; controller design for torque on bottom | Cutting |
| Studying bit bounce, coupled axial–torsional instability, or how blade count, ROP and RPM affect stability | Regenerative |

---

## 5. Issues found in the library code and sample data

These came up while reducing the models. They don't affect the reduced models
above, which follow the published equations, but they explain some of the
sensitivities in §1 and are worth checking in the library.

1. **Detournay & Defourny torque:**
   * `detournay_wob_tob` omits the factor ½ on both $T_c$ and $T_f$, so the
     torque comes out twice the published value.
   * `detournay_rop` inverts the torque with the ½ included, so the two
     functions are inconsistent.
   * `d_avg` is computed but not used, so ROP comes from $d_{WOB}$ only.
2. **Ziaja:**
   * The chord-area expression for $S_c(d)$ does not go to zero as $d \to 0$
     (0.19 in² at $d = 0$) and is almost independent of $d$. This suggests a
     transcription error.
   * The wear-flat term is multiplied by $n_c$ twice: 4,970 lbf instead of
     108 lbf in the sample.
3. **Depth of cut in the cutter models:** Ziaja, Che and Gerbaud use
   ROP = rpm·$d$·$n_c$. This treats every cutter as removing a full layer each
   revolution. In Che, $n_c$ cancels because force is also divided by $n_c$. In
   Gerbaud it does not, so ROP scales roughly as $n_c^2$ (elasticity 2.1).
4. **Gerbaud:** the cut area $A$ is an input rather than being computed from
   the depth of cut, so the calculation is circular, and ROP does not respond
   to WOB.
5. **Sample data units:** with the `tob` column read as lbf·ft, the implied bit
   friction $2T/(aW)$ is 8.5, which is physically implausible. Read as lbf·in it
   is 0.71, which is typical, so the column is probably lbf·in.

## References

* Detournay, E. and Defourny, P., 1992. A phenomenological model for the drilling action of drag bits. *Int. J. Rock Mech. Min. Sci.*, 29(1), 13–23.
* Richard, T., Germay, C. and Detournay, E., 2007. A simplified model to explore the root cause of stick–slip vibrations in drilling systems with drag bits. *J. Sound Vib.*, 305(3), 432–456.
* Ziaja, M.B., 1985. Mathematical model of the polycrystalline diamond bit drilling process and its practical application. *SPE/IADC Drilling Conference*.
* Gerbaud, L., Menand, S. and Sellami, H., 2006. PDC bits: all comes from the cutter/rock interaction. *IADC/SPE 98988*.
* Che, D., Zhu, W.-L. and Ehmann, K.F., 2016. Chipping and crushing mechanisms in orthogonal rock cutting. *Int. J. Mech. Sci.*, 119, 224–236.
* Bit Model Library: SPE-223715-MS, SPE/IADC Drilling Conference 2025.
