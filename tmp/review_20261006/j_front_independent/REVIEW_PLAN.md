# Independent Mode III J and front-step review

This review reads the saved Mode III contour samples and independently sums
the Rice integral. The readback recomputes

    W = mu/2 * (gx^2 + gy^2)
    tz = mu * (gx*nx + gy*ny)
    J = fsum((W*nx - tz*gx)*ds)

from serialized global xyz coordinates, normals, displacement gradients, and
physical ds weights. It does not import the J integrator, field provider, or
solver. It independently reconstructs the rotated ellipse outward normals
and the global Mode III gradient from each point, checks exact front/branch
samples are absent, and checks ordered winding, signed area, normal closure,
branch-face signs, and the analytic KIII^2/(2 mu) value.

The current J artifact identifies itself as an analytic integrator validation
on a manufactured field, not a solved field, measured G, or fracture result.
The criterion here is narrow: all saved per-row terms must recompute from the
external samples; global frame normals and gradients must agree; all eight
contours must be closed/outward/encircle the front; exact front and branch
points must be absent; all independent J values must meet the artifact's
declared 1e-4 quadrature tolerance and higher order must reduce each paired
error.

One source-level contract issue to check and report is whether the 256-row cap
is enforced before any conversion/copy. For the front-step authoring in
progress, acceptance is only that a trial advance may be proposed after
G > Gc. It must not commit propagation solely because a threshold was crossed:
the new candidate geometry/field needs mechanical and energy reanalysis first.
Unknown inputs must remain UNKNOWN, with no fabricated real growth or branch.
