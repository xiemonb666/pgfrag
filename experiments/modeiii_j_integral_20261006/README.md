# Straight-front Mode III J-integral diagnostic

This standalone component computes a contour integral from field samples supplied by its caller. It does not solve for a displacement field and does not infer that supplied data came from a solver. Each result records the caller's provenance label and sets `field_provenance_verified` to `false`.

## Formula and convention

For a straight crack front along global `+z`, a `z`-invariant antiplane field has `u_x=u_y=0` and nonzero `u_z(x,y)`. Rice's contour form is

\[
J=\oint_\Gamma [W n_x-t_z u_{z,x}]\,ds,
\qquad W=\frac{\mu}{2}(u_{z,x}^2+u_{z,y}^2),
\qquad t_z=\mu\nabla u_z\cdot\mathbf n.
\]

The supplied unit normal points outward from the region enclosing the front. This sign convention produces positive `J` for the standard Mode III singular field. Rice's Eq. (10) gives `J=(1+nu) KIII^2/E` for Mode III; since `mu=E/[2(1+nu)]`, this becomes

\[
J=G=\frac{K_{III}^2}{2\mu}.
\]

For coordinates in meters, `u_z` in meters, `grad_u_z` in `m/m` (dimensionless), `mu` in Pa, and line weights `ds` in meters, `W` and each integrand term are in Pa (`J/m^3`). The integral is in `Pa m = J/m^2`.

The primary reference is J. R. Rice, “A Path Independent Integral and the Approximate Analysis of Strain Concentration by Notches and Cracks,” *Journal of Applied Mechanics* 35(2), 379–386 (1968), DOI [10.1115/1.3601206](https://doi.org/10.1115/1.3601206). The paper's author-hosted full text is [Rice's 1968 paper (Harvard)](https://esag.harvard.edu/rice/015_Rice_PathIndepInt_JAM68.pdf). Its Eq. (8a) writes the line integral on a circular path as `r ∫[W cos(theta) - T ∂u/∂x] dtheta`; Eq. (10) gives the Mode III factor relation used above.

## Input contract

`integrate_modeiii_j` accepts four external arrays: contour coordinates `(n,3)` in meters ordered along the contour, outward unit normals `(n,3)`, positive physical line weights `(n,)` in meters, and displacement gradients `grad_u_z=(n,3)` in `m/m`. It also requires positive `shear_modulus_Pa`, the exact `INPUT_UNITS` declaration, and a nonempty `field_provenance` string. It accepts at most 256 samples. The output includes every sample's normal, weight, gradient, energy density, traction, both integrand terms, and weighted contribution so a caller can recompute the sum.

The component returns `UNKNOWN` without a numeric `J` for missing or malformed data, nonfinite data, nonpositive weights or modulus, inconsistent sample counts, non-unit/non-outward normals, a non-closed weighted-normal sum, an ordered contour that does not encircle the front, non-planar contours, non-z-invariant gradients, exact front samples, and exact branch-face samples without an explicit one-sided trace. It never nudges coordinates to avoid a singular or branch point. Geometry checks are diagnostics and do not prove the caller's contour topology.

## Analytic validation artifact

Run `run_analytic_validation.py` with the workspace Python runtime. It reads the fixed `E`, `nu`, and `KIII` values from `modeiii_pu_operator_20261005/fixed_case.json`, then samples the read-only exact field reference from `modeiii_field_reference_20261005/field_reference.py`. The generated contours are two radii, two rotations, and two Gauss-Legendre orders per smooth arc, split exactly at the negative-x branch crossing. All query points are interior to smooth arcs; no point is shifted by epsilon. There are 224 total field samples, within the 256-sample cap. The declared maximum relative-error tolerance is `1e-4`; the actual maximum appears in the JSON.

The output JSON stores the actual input sample arrays and per-sample energy/traction terms. The PNG plots relative error against the analytic `KIII^2/(2mu)` value. This verifies integration on a manufactured field only; it is not evidence of a solved field, validated physical `G`, material calibration, or crack growth.

The six small tests in `test_j_integral.py` cover missing provenance/inputs, invalid units, invalid weights, inward normals, and nonfinite samples.
