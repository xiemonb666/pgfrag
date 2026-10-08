# Rest-frame material field prototype

`rest_material_field.py` defines a pure, three-dimensional field evaluated in a named material rest frame. Its random Fourier directions and phases are generated once from an integer seed and frozen in the field instance. A value depends only on the physical rest coordinate in metres and the immutable field parameters; node IDs, mesh order, query order, and mesh density do not participate in field generation.

For `M` equal-weight modes, the returned dimensionless heterogeneity is

\[
h(\mathbf{x}) = \frac{a}{M}\sum_{j=1}^{M}
\cos(\mathbf{k}_j\!\cdot\!\mathbf{x}+\phi_j),
\qquad |h(\mathbf{x})|\le a<1.
\]

The user supplied `amplitude` `a` is the maximum relative deviation, not an RMS amplitude. With independent uniform phases and isotropic unit directions, the ensemble variance is `a²/(2M)`, and the normalized ensemble radial covariance is

\[
C(r)=\frac{\sin(kr)}{kr}.
\]

The spectral shell radius is `k = q_e / ell`, where `ell` is `correlation_length_m` and `q_e = 2.199123071161498` is the first positive root satisfying `sin(q_e)/q_e = exp(-1)`. Thus `ell` is the first e-fold radius of the ideal isotropic ensemble covariance. It is not a fit to experimental measurements.

The realized field contains a finite set of directions and phases. Its spectrum is a finite collection of lines on one wave-number shell; a fixed realization is anisotropic and oscillatory, and its empirical correlation over a small specimen need not cross `1/e` at exactly `ell`. The sum is smooth, continuous, bounded, and generally quasi-periodic over unbounded space. It is not a Gaussian random field, a Weibull strength distribution, or a calibrated ceramic material model. The finite-mode approximation also does not imply independent material samples at separate points.

Positive caller-supplied parameters use the explicit bounded map `p(x) = p0 * (1 + h(x))`, so `p(x)` lies in `[p0(1-a), p0(1+a)]` and stays strictly positive. The module does not supply nominal values or assert any source for `E`, `Gc`, strength, or other material constants. It records `UNCALIBRATED`, synthetic source lineage, no upstream measured-material artifact, coordinate units `m`, and dimensionless field output units.

When a rigid current-frame transform is provided, coordinates are mapped back to the immutable rest frame before evaluation. Rotating or translating the mesh therefore carries its existing field values with the material. Scaling, shear, frame-ID mismatch, non-metre units, non-finite values, and non-positive correlation lengths are rejected.

## Review artifact

`run_review.py` executes the focused invariance/validation checks and saves a 14 by 14 `z=0` sample grid, a PNG heat map, a detached parameter manifest, and a receipt to `pgfrag/tmp/review_20261006/rest_material_field_20261006/`. The corrected pass accounts for 31 check attempts (including the rejected non-finite query) plus 196 saved sample points. The receipt also records 26 point attempts from an earlier failed implementation pass, for 253 cumulative attempts under the 256-point ceiling. Rendering and receipt hashing do not call the field evaluator.
