# Opening-history independent review design

**Status: DESIGN ONLY.** No opening-control result is pinned here; this note does not authorize or report a numeric run.

## Scope and fixture

Review a saved 12-commit, scalar Mode-I opening history. Treat opening as prescribed input with no applied-force vector. The law's `external_work_j_per_m2` is the signed traction-opening line integral, not a solver's external-force work.

Use a toy fixture, not ceramic calibration: `K0=1e8 Pa/m`, `tpeak=1e3 Pa`, `Gc=1e-2 J/m^2`. This gives `delta0=tpeak/K0=1e-5 m`, `deltaf=2*Gc/tpeak=2e-5 m`, and `deltaf>delta0`. Set one physical interface area to `A=1e-4 m^2`; count it once, so `A*Gc=1e-6 J`. Do not sum both face copies as separate fracture area.

## Independent equations

For each committed state, use finite `delta>=0`, `kappa=max(previous_kappa,delta)`, and the frozen material signature. Recompute

- `d=0` for `kappa<=delta0`; `d=1` for `kappa>=deltaf`; otherwise `d=deltaf*(kappa-delta0)/(kappa*(deltaf-delta0))`;
- `traction=(1-d)*K0*delta` and `psi=0.5*(1-d)*K0*delta^2`;
- `D=0` through `delta0`; `D=Gc*(kappa-delta0)/(deltaf-delta0)` between thresholds; `D=Gc` at and beyond `deltaf`;
- `W=psi+D`, with signed step work `W_i-W_(i-1)`.

For fixed history, independently integrate unloading/reloading work as `0.5*(1-d)*K0*(delta_i^2-delta_(i-1)^2)`. When a new maximum is reached from below the old maximum, include reload-to-old-maximum work and then the exact piecewise-linear envelope integral. Split any comparison path at `delta0`, old maxima, `deltaf`, and turning points. Subdividing a straight segment must leave its summed exact work unchanged within predeclared SI tolerances.

## Twelve-state path

Use an initial virgin state `S0`, then exactly 12 unique committed records in this order:

`[0.5 delta0, delta0, 1.25 delta0, 1.5 delta0, delta0, 0.5 delta0, 1.25 delta0, 1.75 delta0, deltaf, 0.5 deltaf, 0.9 deltaf, 1.1 deltaf]`.

This path covers virgin elasticity, two partial-damage maxima, unload, reload below the prior maximum, new damage, complete failure, unloading, and reopening after failure. The last opening exceeds `deltaf` to expose any healing or residual traction.

## Acceptance gates for a later saved-result audit

1. Bind result, case, receipt, law source, and all state IDs by SHA and exact schema. Require 12 ordered commits after `S0`; reject missing fields, duplicate IDs, nonfinite or negative opening, material-signature changes, or unsupported modes. Do not import the producer law or call a solver.
2. Compare every saved `kappa`, `damage`, traction, `psi`, `D`, and `W` against the equations above. Require `kappa` and damage never decrease; partial states have `0<d<1`.
3. On states 5-6, require damage, `kappa`, and `D` unchanged while `psi` and signed `W` fall. State 7 reloads below the old maximum without adding damage or dissipation. State 8 raises both history and damage.
4. At state 9 require `d=1`, `traction=0`, `psi=0`, and `A*D=1e-6 J`. States 10-12 keep `d=1`, `D=Gc`, `psi=0`, and zero traction; cumulative work stays at the fully dissipated value.
5. Recompute each step's `Delta W=Delta psi+Delta D`; compare exact work for the same path under harmless subdivisions. Use one-side area weights once and keep per-area quantities distinct from joules.

Record tolerances before reading numeric values. A passing result would qualify only this toy opening-history ledger; it would not establish solver coupling, crack growth, mixed-mode response, or ceramic properties.
