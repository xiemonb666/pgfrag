# Enriched DOF registry and provider contract (2026-10-05)

Status: `SOURCE_ONLY / NOT_COUPLED / R2_NONE`. This contract freezes the
identity and shape-row boundary for a first straight-front Mode III prototype.
It does not change the existing MLS kernel or elastic operator, and does not
authorize a solver run, crack-growth run, or promotion.

## Scope and reuse boundary

The first provider is the antiplane scalar field `u_z` for a three-dimensional
cylinder with a straight front. Its exact front function is

```text
F1(x,y) = sqrt(r) sin(theta/2),  r=sqrt(x^2+y^2), theta in (-pi, pi).
grad F1 = (-sin(theta/2), cos(theta/2), 0) / (2 sqrt(r)).
```

The branch cut is the physical negative-x crack face. At `r=0`, the value or
gradient required by a row is unavailable and the provider must return
`UNKNOWN`; it may not sample at `r+epsilon`. On the cut, sided limits must be
associated with the two physical crack-face sides. A local coordinate switch
may not introduce another cut.

The row is one scalar shape contribution and its analytic spatial gradient per
registered coefficient. Those scalars map to physical displacement columns
`(0,0,phi)` and engineering-strain columns
`[0,0,dphi/dz,0,dphi/dy,dphi/dx]`, in the existing order
`[eps_xx,eps_yy,eps_zz,gamma_xy,gamma_yz,gamma_xz]`. Thus the existing
constitutive `D` and `B.T @ D @ B` convention can be reused after building a
matrix with one column per registered coefficient. The old provider's
`3*node_count` layout and its node-dense row contract are not the enriched
coefficient layout.

For base MLS functions and the first front function, the mathematical row is:

```text
base:              phi_i(q), grad(phi_i)(q)
shifted front F1:  phi_i(q) (F1(q)-F1(X_i)),
                   (F1(q)-F1(X_i)) grad(phi_i)(q) + phi_i(q) grad(F1)(q)
unshifted front:   phi_i(q) F1(q),
                   F1(q) grad(phi_i)(q) + phi_i(q) grad(F1)(q)
```

Shifted enrichment is only defined when `F1(X_i)` has a unique physical
branch value. A source on the cut or front cannot be assigned a made-up
`F1(X_i)`; no epsilon offset or side guess is permitted. Such a coefficient
must be unavailable/`UNKNOWN`, or the caller must choose a separately reviewed
unshifted basis. Unshifted use has an explicit rank and approximation gate:
the augmented basis must be checked for additional null modes and linear
dependence, with no pseudoinverse or regularization used to conceal them.

## Registry identity and budgets

The key is exactly
`(source_id, field_component, basis_kind, basis_index)`. In the source module,
`DofKey.canonical_id` encodes this tuple as compact JSON. It is stable and
collision-safe for the same canonical source string, component, basis kind,
and index. `basis_index=0` is reserved for `BASE`; the current non-base kinds
use index 1. `HEAVISIDE`, `FRONT_F1_SHIFTED`, and `FRONT_F1_UNSHIFTED` are
distinct identities. The front kinds register only the first front function
F1. Registration does not establish that an H or front source qualifies
topologically. Equal keys are idempotent. Dense coefficient slots are
assigned only at `freeze()` and sorted by the full key, so registration order
cannot change the slots.

`source_id` identifies the original support source. Registering an enrichment
for an existing source adds a coefficient, not another source node. A source
with BASE and front F1 has one unique source ID and two distinct coefficient
IDs. The registry has hard caps of 375 total coefficients and 125 unique
source IDs. A problem preflight has a 256 quadrature-point cap. Counts must be
validated before copying coordinate arrays or allocating provider rows,
quadrature workspaces, or stiffness matrices; the registry enforces source and
coefficient caps on each addition. `validate_problem_budget` performs the
scalar-count preflight and allocates no arrays.

## Provider protocol and fail-closed rules

`AntiplaneShapeProvider.evaluate(xyz)` returns an `AntiplaneShapeSample` with
one `ScalarShapeContribution(value, gradient)` for every frozen canonical DOF
in exactly the same order. Every registered key must use `field_component=UZ`.
Rows are dense over registry coefficients, not nodes. Each value and each
gradient component is the physical field contribution for a unit coefficient,
including any enrichment product rule; it is not a nodal displacement value.
The gradient is the analytic derivative with respect to reference/material
coordinates.

Only `OK` rows expose numerical contributions. `UNKNOWN` and `BLOCKED` rows
carry no partial IDs or values and stop the whole quadrature/assembly request.
Malformed IDs, dimensions, components, or non-finite values raise
`ProviderContractError`. The caller must not return partial stiffness or
replace unavailable values with zeros. No pseudoinverse, rank truncation, or
regularization is part of this contract.

This protocol validates row structure only. It does not prove that a source
qualifies for H or front enrichment, nor does it choose a Heaviside side
convention or implement its cut quadrature. In particular, the module
intentionally has no `topology_valid: bool` field and does not treat a boolean
assertion as a topology certificate. Before any future coupling, an
independent geometric record must demonstrate that each H support is fully
split by the real finite crack surface and excludes the front, and that each
front support has the required front-topology evidence. The source-only
prototype qualifies no arbitrary source by default. Continuous, uncut MLS
remains mandatory; sharp visibility source dropping and ghost jumps outside
the real crack are disallowed.

## Planned fake-field checks (not run)

The following are acceptance checks for a later small test file. They are
planned only; this source-only task did not run tests or import NumPy.

1. Register one UZ BASE and one UZ shifted-front F1 key for the same source;
   verify two distinct canonical IDs, one source node, stable slots under
   reversed registration order, and idempotent duplicate registration. Also
   verify a HEAVISIDE key has its own ID while registration alone grants no
   source-topology qualification.
2. Register 375 unique coefficients and 125 source IDs successfully; verify
   the next distinct coefficient and the 126th source fail before registry
   mutation. Verify a 376-DOF and 257-quadrature preflight is rejected before
   a fake provider is asked to allocate rows.
3. Use a tiny fake provider to return values and gradients in canonical order;
   verify every enrichment coefficient produces one row column and the scalar
   `u_z` strain mapping matches the six-component engineering-shear order.
4. Return UNKNOWN at a fake front point and BLOCKED at a fake unsupported
   point; verify both rows carry no partial data and the harness stops without
   constructing a partial matrix.
5. Verify wrong coefficient ordering, an `UX` key in the antiplane registry,
   NaN/Inf, wrong tuple sizes, and missing UNKNOWN diagnostics fail closed.
6. Verify no fake topology boolean can qualify a source: topology eligibility
   remains an external, evidence-backed gate and is not set by this module.

## Promotion boundary

Even if the contract and later fake checks pass, the result remains
`NOT_COUPLED`. Independent gates still need to cover support topology, extra
null modes, physical boundary mapping, front/cut quadrature, rank and
conditioning without regularization, quadrature convergence, balance, energy,
and the prescribed Mode III solution. No R2 or ground-truth claim follows from
the registry or its fake-field checks.

## Inputs read

- `pgfrag/tmp/review_20261005/astra_front_shadow_decision_20261005/DECISION.md`
- `pgfrag/experiments/efgm_elastic_operator_20261005/efgm_elastic_operator.py`
- `pgfrag/experiments/mls3d_kernel_20261005/mls3d_shape.py`

