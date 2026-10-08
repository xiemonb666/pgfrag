# Mode III field reference source review

Status: SOURCE_ONLY / ROOT_REVIEW_REQUIRED / TESTS_NOT_EXECUTED

The source implements the accepted first-stage straight-front manufactured
case in global xyz. It has no solver, enrichment provider, CAD, or benchmark
coupling. The query and output limits are 256 rows and 2 MiB. Exact r=0 and
exact negative-x branch-face samples return UNKNOWN; no epsilon shift is used.
End-cap and side-wall loads both come from the full tensor traction sigma*n.

The field is internally consistent as a manufactured solution. Away from the
crack and front, F1 is harmonic in x,y and independent of z, so div(sigma)=0.
The energy density integrates to

    integral_0^L integral_0^(2 pi) integral_0^R
      mu A^2/(8 r) * r dr dtheta dz
    = mu pi A^2 R L / 4
    = KIII^2 R L / (2 mu).

The exact outer traction therefore supplies a self-equilibrated pure-traction
problem when the six rigid modes are constrained by a later solver. This
source does not implement or test those constraints.

Known scope uncertainty: the finite crack reaches both cylinder end caps. The
field is a z-invariant infinite straight-front solution restricted to that
finite cylinder, with manufactured exact loads on the side wall and caps. It
does not describe a crack that terminates inside a larger medium, and the
domain contains no exterior points beyond the crack ends on which to test an
end continuation. The off-crack continuity test is limited to a smooth path
away from the cut and front. The source is not a finite-cylinder material
calibration, curved-front solution, or fracture benchmark.

Prepared tests cover central-difference gradients away from the cut, face
traces, positive-x and off-crack continuity, one-sided face traction, all
external boundary traction classes, and an independent polar quadrature
energy reference. They were not run. The next permitted action is one focused
root-authorized test after this source review.
