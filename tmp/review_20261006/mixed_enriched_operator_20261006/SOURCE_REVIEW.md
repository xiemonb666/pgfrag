# Source review

The primary physical layout is explicitly `global_front`: one `SOURCE_BASE` coefficient per source and displacement component, followed by four `GLOBAL_FRONT` functions per component. At 27 sources this is 31 scalar functions and 93 vector coefficients. The adapter extracts each `N_s` and its gradient once from the reviewed per-source helper row, then imports the four common branch values/gradients once. It does not construct 27 copies of each global branch and delete columns afterward.

The separate `per_source_front_diagnostic` layout maps all five per-source scalar terms and produces the full 405-column raw matrix. It is for visible diagnostic comparison only. Under linearly reproducing base functions, two branch combinations give exact zero fields per vector component:

```text
x F3 - y F4 + y F1 = 0
x F4 + y F3 - y F2 = 0
```

`branch_dependency_mode_coefficients` encodes the corresponding six candidate coefficient vectors. `rigid_body_mode_coefficients` returns the six physical translation/rotation candidates for either layout. Neither function certifies null-space rank. The 93-DOF spectrum still must be inspected under a separately authorized mechanical validation before physical eligibility can be claimed.

Rows are accepted only when status, canonical source order, scope-tagged semantic IDs, contribution order, finite shapes/gradients, and independently reconstructed engineering-strain B columns agree. Positive volume/surface weights and finite SI material/traction values are required. Isotropic `C` is full 6x6; supplied `C` is checked symmetric positive definite using a small Cholesky validation. `K` and force are serialized as immutable little-endian float64 bytes with units N/m and N, respectively. The implementation returns raw accumulation and includes no K symmetrization, ridge, pseudoinverse, load projection, KKT, solve, eigentruncation, provider, MLS, or fracture model.

The frozen full-3D route encounters two separate explicit 375 guards: `pgfrag/experiments/enriched_dof_contract_20261005/enriched_dof_contract.py::MAX_TOTAL_DOF` through `validate_problem_budget`, and `pgfrag/experiments/efgm_elastic_operator_20261005/efgm_elastic_operator.py::MAX_DOF` inside the reused array assembler before K allocation. The full-3D wrapper invokes both. They are source-level resource bounds, not consequences of B/C/K mathematics; both files and their caps remain untouched.

`operator_arithmetic_check.py` is a planned independent checker source and was not run. Static syntax/JSON parsing is the only allowed validation in this stage. The checker does not compute an eigen-spectrum; the root's six-zero 93-DOF rank gate remains pending.
