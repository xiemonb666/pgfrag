# 3-D linear MLS primitive

This folder contains a bounded Python/NumPy mathematical primitive for linear
moving-least-squares interpolation in material coordinates. It accepts plain
NumPy arrays with at most 500 sources, 500 queries, and three field components;
budget violations raise `ValueError` before conversion or output allocation.
It computes field values and analytic spatial gradients. It is not a CPM solver,
a crack or fracture model, or a production mechanics backend.

`mls3d_mls.py` is the implementation. `test_mls3d_mls.py` records source-only
acceptance checks. Tests and solver execution are intentionally deferred until
the root agent reviews this freeze.

Read [CONTRACT.md](../tmp/review_20261005/mls3d_kernel_20261005/CONTRACT.md)
and [SOURCE_REVIEW.md](../tmp/review_20261005/mls3d_kernel_20261005/SOURCE_REVIEW.md)
for behavior, scope, and source qualifications.
