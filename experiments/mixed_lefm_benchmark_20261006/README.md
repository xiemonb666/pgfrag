# Mixed-mode straight-front LEFM source

`mixed_lefm_field.py` is a pure-standard-library analytic reference for an isotropic, finite-cylinder query region around an extruded straight crack. Local axes are `(x, y, z)`: the crack surface is `x < 0, y = 0`, propagation points toward `+x`, and the straight front follows `+z`. The field is the classical leading Williams field and stays independent of `z`; the finite radius and length only bound queries. It does not represent a finite crack front turning at the cylinder end caps.

For Mode I and II, with `kappa = 3 - 4 nu` (plane strain) and `mu = E / (2 (1 + nu))`, the evaluator uses

```text
u_x = KI  / (2 mu) sqrt(r / 2 pi) cos(theta/2) (kappa - cos(theta))
    + KII / (2 mu) sqrt(r / 2 pi) sin(theta/2) (kappa + 2 + cos(theta))
u_y = KI  / (2 mu) sqrt(r / 2 pi) sin(theta/2) (kappa - cos(theta))
    - KII / (2 mu) sqrt(r / 2 pi) cos(theta/2) (kappa - 2 + cos(theta))
u_z = KIII / mu sqrt(2 r / pi) sin(theta/2)
```

Each term is differentiated analytically as `sqrt(r) h(theta)`, using `dr/dx = cos(theta)`, `dr/dy = sin(theta)`, `dtheta/dx = -sin(theta)/r`, and `dtheta/dy = cos(theta)/r`; no finite differences or near-tip cutoff are used. The local displacement gradient is rotated with `Q grad(u) Q^T`; displacement, strain, and stress tensors are rotated with the validated proper frame `Q`. The 3-D isotropic constitutive law is applied to the symmetric gradient, so plane-strain `epsilon_zz = 0` produces the full `sigma_zz`, and Mode III retains `tau_xz` and `tau_yz`. Engineering strain Voigt order is `[xx, yy, zz, gamma_xy, gamma_yz, gamma_xz]`; stress order is `[xx, yy, zz, xy, yz, xz]`.

The four crack-tip branch basis values are normalized as `sqrt(r/R) * [sin(theta/2), cos(theta/2), sin(theta/2) sin(theta), cos(theta/2) sin(theta)]`. Thus they are dimensionless and their analytic local gradients are in `1/m`; `R` is the caller-supplied positive reference length. Their associated displacement coefficients have units meters.

Exact front queries return `UNKNOWN_FRONT` without field values. An exact negative-x face query without a side returns `UNKNOWN_BRANCH_FACE`; callers can explicitly request `upper` (`theta = +pi`) or `lower` (`theta = -pi`) for a one-sided trace. In a rotated frame, use `evaluate_local_point((x, 0.0, z), case=..., branch_trace="upper")` or `"lower"` to preserve exact local `y = 0`; it returns the field tensors rotated to world coordinates. `evaluate_point` uses exact comparisons after its world-to-local transform: an explicit trace is rejected if that round trip does not land exactly on the face, and an unqualified point is classified as a face only at exact local `y = 0`. No tolerance is used to infer or snap a trace. The implementation never chooses a trace from signed zero and never shifts `r = 0`. On either explicit trace, traction on the physical crack face (normal to local y) is zero analytically; the full stress may still contain tangential components. Out-of-cylinder points return `OUTSIDE_DOMAIN`. Nonfinite/invalid SI inputs, frames, or branch selectors are rejected.

## Formula sources

- The Mode I/II Cartesian displacement equations and plane-strain `kappa` are verified against Li and Maurini, “Crack kinking in a variational phase-field model of brittle fracture with strongly anisotropic surface energy,” *Journal of the Mechanics and Physics of Solids* 125 (2019), Eq. (A.1), PDF p. 519: [paper PDF](https://www.lmm.jussieu.fr/~corrado/articles/cm-a19-LiMau.pdf). The same signs and factors are independently checked against A. F. Bower, *Applied Mechanics of Solids*, §9.3.1, displayed displacement expressions and its note that in-plane displacements apply to plane strain: [§9.3](https://solidmechanics.org/Text/Chapter9_3/Chapter9_3.php).
- The full leading rectangular stress expressions, Mode III displacement, and plane-strain scope are given in Bower, *Applied Mechanics of Solids*, §9.3.1, same [source page](https://solidmechanics.org/Text/Chapter9_3/Chapter9_3.php). This implementation derives all six stress components instead from the analytic full gradient and 3-D isotropic Hooke law, preserving the plane-strain `sigma_zz` term.
- The four branch functions are Moës, Dolbow, and Belytschko, “A finite element method for crack growth without remeshing,” *International Journal for Numerical Methods in Engineering* 46 (1999), Eq. (16), journal pp. 135–136: [publisher DOI record](https://doi.org/10.1002/%28SICI%291097-0207%2819990910%2946%3A1%3C131%3A%3AAID-NME726%3E3.0.CO%3B2-J), [author-uploaded full text with Eq. (16)](https://www.researchgate.net/publication/51992357_A_Finite_Element_Method_for_Crack_Growth_without_Remeshing).
- The isotropic mixed-mode energy release identity `G = (KI^2 + KII^2)/E' + KIII^2/(2 mu)`, with `E' = E/(1 - nu^2)` for plane strain, is stated as Eq. (mixedModeG) in the University of Liège CM3 fracture mechanics course notes: [equation page](https://www.ltas-cm3.ulg.ac.be/FractureMechanics/print.php?p=Lecture3_P4). Bower §9.3.4 independently gives the plane-strain Mode I relation `G = (1 - nu^2) KI^2 / E`: [§9.3](https://solidmechanics.org/Text/Chapter9_3/Chapter9_3.php).

This is an elastic verification primitive only. `G` is not compared with a toughness, and the module contains no material calibration, crack initiation rule, or autonomous fracture behavior. No solver, provider, sample, NumPy, or test run is part of this source preparation.
