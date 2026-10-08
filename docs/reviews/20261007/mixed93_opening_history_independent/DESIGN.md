# Independent review of the global opening-history result

**Status: run01 checker attempt was NOT_ACCEPTED; corrected source is ready for a possible run02 after separate root review and exact GO. No accepted audit result exists.** This is the actual fixed-face 93-DOF control path, not a single-point scalar trajectory. The producer result and receipt are bound by path and SHA-256. Run01 is preserved and failed with `KeyError: analytic_alpha` because the saved interface schema does not contain that field. The checker now derives analytic opening from its reconstructed normal row and keeps the producer rigid-trace array bound exactly before applying the existing physical trace tolerance. The exclusive independent output path is `run02/result.json`; execution remains disabled pending a separate exact root GO.

## Frozen input contract

The producer is the 468-line `run_mixed93_opening_control.py` source with its frozen case and source-pin set. Read its saved result as JSON only and verify it against the exact driver, case, saved elastic input, and receipt hashes. Do not load NPZ, import producer/author modules, run a solver/eigensolver, or issue provider queries.

Reconstruct all 16 one-sided interface points from the saved elastic trace. For each point, verify `J=Hplus-Hminus`, the saved normal row, and `alpha_i=2*sqrt(-x_i/R)`. The prescribed vector is `q(Q)=Q*q_per_Q`, where `q_per_Q` and physical gauge multiplier per Q come from the saved `K_raw` and `C_phys` response. `Qstar=deltaf/min(alpha_i)`.

Use the frozen 12 fractions, in order: `[0, .05, .15, .30, .12, .30, .50, .75, 1.10, .65, 0, 1.10]*Qstar`. Body force, boundary force, and physical gauge target are zero. The index-82 control reaction is nevertheless nonzero and is recomputed as `(K*q + C_phys.T*lambda_phys + f_gamma)[82]`; noncontrol entries must satisfy the frozen residual tolerance. Do not call this an unloaded/free-reaction path.

## Pointwise law and independent work

Use the frozen toy values `K0=1e8 Pa/m`, `tpeak=1e3 Pa`, `Gc=1e-2 J/m^2`, `delta0=1e-5 m`, `deltaf=2e-5 m`. Keep per-area values in J/m² and integrated values in J. Sum the 16 one-side weights once: `A=R*L=1e-4 m^2`, so complete-failure dissipation is `A*Gc=1e-6 J`.

At each point and step, recompute `delta_i=alpha_i*Q_i`, `kappa_i=max(kappa_(i-1),delta_i)`, damage `d=0` through `delta0`, `d=1` at/above `deltaf`, otherwise `d=clamp(deltaf*(kappa-delta0)/(kappa*(deltaf-delta0)),0,1)`. Then `t=(1-d)*K0*delta`, `S=.5*(1-d)*K0*delta^2`, and `Diss=I_envelope(kappa)-S(kappa)`, where `I_envelope` is independently evaluated with a piecewise antiderivative split at `delta0` and `deltaf`.

Recompute each unload/reload interval exactly from fixed committed damage and each new-maximum interval as reload to the old maximum plus the exact envelope integral. Compare per-point work increments, cumulative work, and area-weighted primitive totals without importing the scalar primitive. Require nondecreasing `kappa`, damage, and dissipation at every point.

## Independent gates

1. Validate 12 ordered steps × 16 unique QP IDs, exact `trajectory_binary_le_f64` format and shapes, strict base64 decode, byte lengths, and every little-endian float64 against its JSON trajectory value (`q_physics`, damage, `kappa_m`, opening). Never silently reshape, repair, or project data.
2. Recompute `K*q`, `C_phys.T*lambda`, and `f_gamma=sum(A_i*J_i.T*n*t_i)` from saved matrices and point data. Check the saved bulk/gauge/cohesive vectors, control reaction sign, and noncontrol force residual. Check upper/lower physical traction signs and equal/opposite face forces.
3. Require partial damage at step 3, preserved history through unload step 4 and reload step 5, all points fully damaged at step 8, zero opening with no healing at step 10, and zero traction with retained damage/dissipation on reload step 11.
4. Recompute bulk storage `0.5*q.T*K*q`, exact bulk step work `0.5*(q_per_Q.T*K*q_per_Q)*(Q_i^2-Q_(i-1)^2)`, and the exact area-weighted cohesive interval work. Check each energy increment against change in bulk storage + cohesive storage + cohesive dissipation using relative tolerance `2e-9` and absolute tolerance `1e-22 J`.
5. Independently require final area-weighted dissipation `1e-6 J`; verify exact path work is invariant when the same stored piecewise path is algebraically subdivided at its existing branch/turning points. This subdivision check is bookkeeping only, not solver or time-step convergence.

Any future pass is diagnostic for this toy fixed-face path. It does not qualify ceramic material data, autonomous crack growth, mixed mode, or general physical accuracy.
