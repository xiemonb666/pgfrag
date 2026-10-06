# Source review: mixed-mode LEFM field

Status: source prepared; numerical validation remains for the main agent to authorize.

## Scope and provenance

The permitted change is limited to the experiment evaluator and its formula note under `pgfrag/experiments/mixed_lefm_benchmark_20261006/`. It is a bounded query wrapper for the leading field of a locally straight, z-extruded crack. It contains no solver/provider calls, sample cases, automatic fracture law, toughness input, or calibration.

Formula anchors reviewed in the source pages:

- Li and Maurini (2019), Eq. (A.1), PDF p. 519: Mode I/II Cartesian displacement field and the plane-strain value `kappa = 3 - 4 nu`; [paper PDF](https://www.lmm.jussieu.fr/~corrado/articles/cm-a19-LiMau.pdf).
- Bower, *Applied Mechanics of Solids*, §9.3.1: full leading rectangular stress and displacement expressions, Mode III `u3`, and explicit plane-strain scope for in-plane displacement; [textbook section](https://solidmechanics.org/Text/Chapter9_3/Chapter9_3.php).
- Moës, Dolbow, Belytschko (1999), Eq. (16), journal pp. 135–136: four crack-tip branch functions; [publisher DOI record](https://doi.org/10.1002/%28SICI%291097-0207%2819990910%2946%3A1%3C131%3A%3AAID-NME726%3E3.0.CO%3B2-J) and [author-uploaded full text at Eq. (16)](https://www.researchgate.net/publication/51992357_A_Finite_Element_Method_for_Crack_Growth_without_Remeshing).
- University of Liège CM3 fracture mechanics course, Eq. (mixedModeG): isotropic mixed-mode energy release identity and plane-strain `E'`; [course equation](https://www.ltas-cm3.ulg.ac.be/FractureMechanics/print.php?p=Lecture3_P4). Bower §9.3.4 independently gives the plane-strain Mode I special case.

## Critical source-level checks

- The `sqrt(r) h(theta)` derivative is implemented by its explicit polar chain rule: `dx = A/sqrt(r) [0.5 cos(theta) h - sin(theta) h']`; `dy = A/sqrt(r) [0.5 sin(theta) h + cos(theta) h']`. Each Mode I/II/III angular derivative is written analytically.
- Normalized branch values apply `q = sqrt(r/R)` to each angular function; the gradients use the derivative of that same `q h(theta)` field, `q/r` times its polar angular/radial chain-rule terms.
- A finite cylinder clips queries to `r <= R`, `|z| <= L/2`. The field remains z-invariant and does not assert a physical finite crack-front solution at the end caps.
- At exact `r = 0`, field values are withheld with `UNKNOWN_FRONT`. At `x < 0, y = 0`, absent an explicit side the field returns `UNKNOWN_BRANCH_FACE`; `upper/lower` force `theta = +pi/-pi`. The code never infers a side from signed zero or shifts the point.
- `evaluate_local_point` accepts exact crack-frame coordinates and returns world-frame tensors. The world wrapper keeps exact comparisons and never adds a tolerance if an arbitrary rotation round-trip perturbs local `y = 0`.
- Local `epsilon = sym(grad u)` is passed through full 3-D isotropic Hooke law. The Mode I/II displacement has `u_z = 0`, hence `epsilon_zz = 0` and `sigma_zz = lambda (epsilon_xx + epsilon_yy)`; Mode III contributes `tau_xz` and `tau_yz`. The physical crack-face traction uses the local-y normal and its three components vanish on both explicit traces.
- For proper `Q`, vectors use `Q u`, gradients use `Q grad(u) Q^T`, and symmetric tensors use `Q T Q^T`. Frame construction rejects nonorthogonal and reflected axes.
- With `K` in Pa sqrt(m), `E'` and `mu` in Pa, each energy-rate term `K^2 / modulus` has units Pa m = J/m^2. Branch functions use `sqrt(r/R)`, so their values are dimensionless and their gradients are 1/m.

This review is algebraic/source-only. No numerical probes, NumPy imports, provider/solver execution, sample execution, or pytest were run. Face traction, rigid-rotation covariance, stress/Voigt consistency, and the energy identity still need the later small finite validation requested by the main agent.
