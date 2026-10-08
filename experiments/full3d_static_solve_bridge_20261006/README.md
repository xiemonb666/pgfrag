# Full 3-D static solve bridge

This bridge is the next step after the full-vector elastic operator. It will
consume one frozen operator diagnostic saved result and perform only the
bounded static solve. It will not import or call the MLS provider, field
evaluator, quadrature builder, or stiffness assembler. The saved artifact must
carry the operator, load, semantic registry, six rigid coefficient vectors,
the 384 unique quadrature rows, basis values and gradients, exact field data,
full-vector boundary traction, statuses, units, and source hashes. The exact
field contract is `static_solve_replay.v1`; the operator source pins and result
location are frozen, and the upstream operator diagnostic has run once. The
bridge input is
`tmp/review_20261006/full3d_cracked_elastic_operator_20261006/run01/result.json`.
The case pins its SHA-256. Run01 made one bounded attempt, then blocked during
saved-registry validation because `basis_row_ids` is serialized directly under
`static_solve_replay`, not under its `registry` object. The child imported NumPy
on entering the numeric worker, but stopped before gauge assembly, spectrum
analysis, or a direct solve. The reader now follows the pinned wire layout.
Run02 consumed only the saved operator and field arrays, added no queries, and
completed one direct 168-by-168 KKT solve.

The one fixed case has 162 physical displacement coefficients: 27 sources,
BASE and global `FRONT_F1_UNSHIFTED`, and UX/UY/UZ. The manufactured field is
the Mode III cylinder benchmark with `R=L=0.01 m`, `E=1e6 Pa`, `nu=0.3`, and
`KIII=100 Pa*sqrt(m)`. Its exact coefficient vector is constructed through
the saved, ordered registry: all BASE entries are zero, all UX/UY enrichment
entries are zero, and each UZ enrichment entry is `A*sqrt(R)`. This remains a
synthetic manufactured case, not ceramic data, a local-front method, or a
moving-front model.

The bridge uses the author's exact saved volume quadrature: GL4 x GL16 x GL4,
256 points, with its saved physical weights and saved basis rows. The load
comes from the same diagnostic's 128 side-and-cap boundary rows, including all
traction components. It adds no off-grid points and makes no new provider
queries. The 384 saved positions are the total location budget.

The six gauge rows are physical-volume moments. For the saved volume centroid
`xbar`, define three translations `r=e_i` and three normalized rotations
`r=omega_i x (x-xbar)/R`. With `V=sum_q w_q`, build
`C_lj=(1/V) sum_q w_q r_l(x_q) dot N_j(x_q)` and
`g_l=(1/V) sum_q w_q r_l(x_q) dot u_exact(x_q)`. Thus `C` is dimensionless,
`g` is in metres, and the physical multipliers are in newtons. The nonzero
manufactured target is retained. For this positive Mode III amplitude the
`r_rx` target is positive; its continuum value is
`16*A*sqrt(R)/(21*pi)`, while the recorded target is computed from the saved
256-point rule. This gauge fixes only the six rigid-body representatives; it
does not alter the applied load. The saved rigid-mode rotation center must be
present so the bridge can express the six columns about the actual volume
quadrature centroid before checking the gauge and wrench.

Before the one direct solve, the saved K must resolve exactly six expected
rigid null modes and 156 positive modes, with no significant negative modes.
Any UNKNOWN source row or extra/negative mode blocks the solve. The code may
diagonally scale K and the six constraint rows for numerical stability, while
retaining the raw K, C, f, and g. It may not pin coefficients, project loads,
delete eigenvalues, use a pseudoinverse, or regularize the matrix.

The constrained equations are

```text
[ K  C.T ] [ d      ] = [ f ]
[ C   0   ] [ lambda ]   [ g ]
```

Report the raw weak-equilibrium residual, `C*d-g`, `U=0.5*d.T*K*d`,
`W=d.T*f`, and `2U-W+g.T*lambda`; the expected identity is
`2U-W=-g.T*lambda`. Also report all six rigid-vector dot-load values beside
the integrated side, top, bottom, and total force and moment. Compare the
solved displacement and six engineering strains against the exact field only
at the saved 256 volume rows. Keep that solved field separate from the exact
manufactured coefficient reconstruction.

The planned solve attempt is one run with a 60-second wall limit, an 8 MiB
explicit numeric-array cap, and an 8 MiB result cap. This directory contains
the saved-result reader and the frozen-case contract. The pinned upstream
source manifest is available in
`tmp/review_20261006/full3d_cracked_elastic_operator_20261006/SOURCE_PINS.json`;
the matching operator result is pinned in `static_solve_case.json`. The bridge
run record and run01 receipts are retained under
`pgfrag/tmp/review_20261006/full3d_static_solve_bridge_20261006/`. Run02
returned `SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED` in 0.754 seconds including the
outer process (0.547 seconds in the child). Its result is
`pgfrag/tmp/review_20261006/full3d_static_solve_bridge_20261006/run02/static_solve_bridge_result.json`,
SHA-256 `ef66783e2961dd77f34785a99e1b397210ccc0beb6e8fa2979678aa4992a0ce3`.
The spectrum resolved six rigid modes and 156 positive modes with no
significant negative modes; four strict negative eigenvalues were within the
roundoff tolerance. The nonzero `r_rx` gauge target was `5.031213792347516e-6`
m. Weak-equilibrium and gauge residual norms were `1.70e-16 N` and
`5.35e-21 m`; `U=1.3133527026982052e-6 J`, `W=2.626705405396418e-6 J`, and
`g^T lambda=4.236127247359188e-21 J`, with energy-identity defect
`2U-W+g^T lambda=-3.39e-21 J`. Solved-field relative errors at the saved
volume rows were 1.90% in displacement and 5.31% in strain. Run02 did not
invoke a provider or reassemble the operator, and remains a synthetic
diagnostic rather than a material or fracture qualification.
