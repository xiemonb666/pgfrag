# Local front support topology safe v2

This is a source-only, separately versioned classifier copied from
`local_front_support_contract_20261006/support_topology.py`; provenance and the
specific edits are recorded in the module header and the focused review receipt.
The frozen provider, its manifest, and the original experiment are untouched.

`Omega` is a closed box and a support is `Omega ∩ closed_ball(center, h)`. The
finite crack is `[x_min, front) × {y_plane} × [z_min, z_max]`; its start and z
ends are included and its front edge is excluded. The finite front is the
closed segment `{front} × {y_plane} × [z_min, z_max]`. H requires the whole
box-clipped plane section inside the finite crack, two nonempty y sides, and
positive clearance from the closed front. Exact equality at a closed crack
endpoint that is also clipped by the box is accepted. Nonzero endpoint gaps
within the topology tolerance and every actual overrun return `UNKNOWN`. No
epsilon offset or coordinate snapping is used.

Run `focused_geometry_check.py` to reproduce the one bounded geometry check;
it reads the original fixed case as data and calls only this source-only
classifier.
