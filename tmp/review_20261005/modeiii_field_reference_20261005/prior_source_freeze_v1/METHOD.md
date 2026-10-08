# Straight-front Mode III field reference

This experiment defines a manufactured analytic field for source-level checks
before any solver coupling. The cylinder has radius R = 0.01 m and length
L = 0.01 m. The crack surface is the negative-x half-plane at y = 0 inside
-L/2 <= z <= L/2. Its straight front is the z axis. The fixed data are
KIII = 100 Pa sqrt(m), E = 1e6 Pa, and nu = 0.3.

Use global xyz coordinates and

    r = sqrt(x^2 + y^2), theta = atan2(y, x), -pi < theta < pi
    F1 = sqrt(r) sin(theta/2)
    grad(F1) = [-sin(theta/2)/(2 sqrt(r)),
                 cos(theta/2)/(2 sqrt(r)), 0]
    mu = E / (2 (1 + nu))
    A = KIII / mu sqrt(2/pi)
    u = (0, 0, A F1)
    tau_xz = mu A dF1/dx
    tau_yz = mu A dF1/dy
    W = mu A^2 / (8 r)
    U = KIII^2 R L / (2 mu)

The stress tensor is symmetric, with only xz and yz shear entries. For any
unit physical normal n, the applied traction is sigma*n. This gives generally
nonzero traction on the cylindrical side wall and both end caps. On either
crack-face trace the traction tends to zero. The field is equilibrated because
F1 is harmonic away from the crack and front; its r^-1/2 gradient has finite
integrated strain energy around the straight front.

The principal atan2 branch cut is the negative x axis, which exactly coincides
with the physical crack surface throughout this bounded cylinder. A query on
the front or exact branch face is UNKNOWN. The implementation never chooses a
side from signed zero and never offsets or regularizes r=0. Points beyond the
closed cylinder are OUTSIDE_DOMAIN. The geometry has no admissible points
beyond the finite crack ends, so this case cannot test continuation into an
exterior medium past an end cap. An off-crack path at x=0, y>0 checks that no
artificial shadow jump is introduced there.

The prepared companion tests check analytic derivatives by central
differences away from the cut, the sqrt(r) face jump, positive-x and off-crack
continuity, zero face-traction limit, side-wall and both end-cap tractions,
and total energy by polar quadrature. Their source is prepared but not
executed pending root review. No proper rotation, curved front, plasticity,
dynamics, material calibration, solver, CAD, or numerical run is included.

The implementation accepts at most 256 finite coordinate queries in one
batch and at most 2 MiB of output. Inputs are plain C-contiguous NumPy arrays.
SI case data are bounded before evaluation; non-representable float64 results
raise an error instead of being clipped.
