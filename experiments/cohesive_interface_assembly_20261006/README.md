# Scalar Mode-I cohesive interface assembly source

This bounded source module connects the existing scalar cohesive law to a
mechanical interface residual and tangent. At each physical quadrature point,
the face interpolation matrices `H+` and `H-` have shape `3 x ndof`; the module
defines `J = H+ - H-`, `jump = J d`, and `delta_n = n . jump`. It adds
`w J.T n t` to the resisting force gradient and
`w (J.T n) k_alg (n.T J)` to the consistent tangent.
Global displacement DOFs are ordered `(ux, uy, uz)` per node.

`one_side_area_weight_m2` counts one physical interface side once. Positive and
negative face nodal forces are returned separately and their resultants are
checked equal and opposite. Explicit positive/negative geometric points and
global nodal positions are required. Before assembly, each face operator must
reproduce rigid translations and its declared point, the face points must
coincide, and `J` must annihilate all three rigid rotations.
This nodal primitive is capped at 12 nodes / 36 DOFs and 16 quadrature points;
both limits reject before allocating the dense `ndof x ndof` tangent. It does
not support the separate 162-enriched-DOF coupling path. Point/translation
checks use `1e-12 m` absolute plus `1e-10` relative tolerance; rigid-rotation
residual uses `1e-15 m` absolute plus `1e-10` relative to the rigid field scale.
Resultant balance uses a separate `1e-12 N` absolute plus `1e-10` relative
force tolerance.

Every call evaluates `update_state` from the supplied committed immutable state
for that point. The trial result includes trial states and energy totals but
does not write them back. Repeating a Newton trial with the same committed
mapping is read-only. The caller controls any future commit boundary.

All input/output quantities use SI units: nodal positions and displacement in
metres, unitless normal and `J`, area in m², traction in Pa, force in N, tangent
in N/m, and integrated energies in J. History keys are unique
`(interface_id, quadrature_id)` pairs, and the history mapping must match the
quadrature point IDs exactly. Compression, tangential slip, and mixed-mode
opening raise `UnsupportedInterfaceModeError` and are UNKNOWN in this module.

The tangent is symmetric for this scalar normal law but may be indefinite in
the softening range. At history equality, `delta0`, and `deltaf`, results mark
the tangent nondifferentiable and identify the selected one-sided branch. This
module contains no solver, Newton loop, crack-front coupling, or calibration;
the shared material fixture remains illustrative and is not a ceramic fit.
No numerical mechanics checks were run in this source-only pass. Before a
mechanics acceptance, add finite-difference tangent checks away from kinks,
repeat-trial idempotence checks against one committed state, force-balance
checks, and one-side `Gc`/energy-ledger checks.
The proposed checker source is `check_interface_numerics.py`; root has not
authorized its execution. It declares one QP for the four-node paired-face
toy, 17 `assemble_trial` calls total (9 directional-FD calls across elastic,
softening, and unloading; 1 repeat; 1 full failure; 2 unsupported modes; 2 ID
failures; 2 pre-allocation cap failures), 11 scalar update/traction calls
inside assembly, plus one direct scalar update to prepare committed unload
history. Cap rejection cases use 17 QPs and 13 nodes. The declared runtime
budget is 10 s maximum. FD uses `h=1e-10 m`, `atol=50 N/m`, `rtol=2e-5`;
force uses `1e-12 N + 1e-10 relative`, moment `1e-12 N·m + 1e-10 relative`,
and energy `1e-10 J + 1e-12 relative`. These are planned checks, not results.

`check_source_contract.py` uses only the Python standard library to parse this
source and inspect the JSON contract schema. It does not import or execute the
assembler, scalar helper, or any mechanics path. The run receipt is saved under
`pgfrag/tmp/review_20261006/cohesive_interface_assembly_20261006/`.

The bilinear opening envelope and area relation follow Alfano et al. (2009),
[Eqs. 5–6](https://paulino.princeton.edu/journal_papers/2009/IJF_09_ModeIFracture.pdf).
Zhang & Paulino (2005), [Eqs. 26–28](https://paulino.princeton.edu/journal_papers/2005/IJP_05_CohesiveZoneModeling.pdf),
document a bilinear cohesive relation and maximum-separation history for
irreversibility. Hillerborg et al. (1976), [original paper](https://www.codil.or.kr/filebank/original/RK/OTMCRK500148/OTMCRK500148.pdf),
is the fictitious-crack/fracture-energy context. These sources inform the law;
they do not calibrate the fixture values.
