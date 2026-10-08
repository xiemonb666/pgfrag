# Independent review: straight-front Mode III J

## Finding

The saved contour data pass an independent recomputation of the Mode III
Rice integral. The calculation used only the serialized coordinates,
global normals, global displacement gradients, line weights, and saved
material constants. It did not import or call the J-integral implementation,
analytic field provider, or solver. Per-row terms were recomputed and summed
with Python math.fsum.

The eight cases contain 224 samples total. The reference is
KIII^2/(2*mu) = 0.013 J/m^2. Independent values agree with the serialized
J values to within a maximum relative difference of 1.33440323226758e-16; the maximum absolute
relative error against the analytic value is 8.679744315207977e-5, below
the artifact's stated 1e-4 quadrature tolerance. Raising arc order from 12
to 16 reduces the absolute error for every radius and rotation. Doubling both
ellipse semiaxes leaves the result invariant to roundoff.

The global-frame checks also pass. Reconstructing the exact Mode III gradient
at every saved coordinate gives a maximum absolute difference of
4.84869951806108e-19 m/m; reconstructing each rotated ellipse's outward
normal in global xyz gives a maximum vector difference of
3.1401849173675503e-16. Each contour has both negative-x crack-face traces,
no exact front or branch-face sample, unit in-plane normals, positive
outward signed area, weighted normal closure below 1.8e-16 relative, and
2*pi winding. The independent values and per-case checks are in
independent_readback.json; the bounded execution receipt is in readback_run/.

This establishes only that the integrator reproduces the expected contour
integral on the supplied manufactured analytic field samples. The saved
artifact is explicitly ANALYTIC_INTEGRATOR_VALIDATION_NOT_SOLVED_FIELD;
its caller-provided provenance remains unverified. It is not a solver result,
measured G, material calibration, or crack-growth result.

## Source issue requiring correction

The documented 256-sample maximum is checked after conversion and private
copy. _numeric_array calls astype(np.float64, copy=True) at
j_integral.py:83, while integrate_modeiii_j first rejects
n > MAX_SAMPLES at j_integral.py:137. Thus four oversized inputs can be
copied before the cap is applied. np.asarray at line 71 can also allocate
when a caller supplies a large list-like object. Preflight each input's row
count, byte count, shape, and common sample count before any conversion or
copy; for list-like input, either preflight sequence dimensions safely or
require bounded plain ndarrays. Add a regression that rejects 257 rows before
any allocation/copy hook is invoked.

## Narrow acceptance

Keep the current analytic contour result as a formula/integration check only.
Before treating this module as a bounded external-data interface, move the
sample and byte guards ahead of conversion and copy. The current saved runs do
not sample the exact negative-x branch face, so they do not validate the
one_sided_crack_face input. That argument is one label for the whole batch;
if a future contour supplies exact points from both lips, one scalar cannot
identify the side of every row. A safe narrow rule is to retain UNKNOWN for
exact cut points, or add an explicit side label per sample and test the two
limits independently.

For the front-step authoring in progress, a positive G-Gc may only create a
trial proposal. It must not commit the advance on threshold crossing alone.
Commit requires reanalysis of the proposed geometry, mechanics, and energy
balance. Missing or UNKNOWN inputs cannot be converted into a real advance,
branch, or passed energy check.

## Independent front-step review

The current front-step source and its saved execution receipt are hash-bound:
the source, test, fixture runner, four serialized fixture states, and SVG each
match the SHA-256 recorded in `ACTUAL_TEST_RECEIPT.json`. The receipt reports
15 fixture tests, four fixture states, no solver/native runtime, and no
production writes. This confirms which source the fixture run exercised; it
does not turn the fixtures into mechanics evidence.

The transition rule is appropriately narrow for a source-only state machine.
`G < Gc` arrests, equality does not propose, and `G > Gc` creates a pending
geometry candidate while retaining the old active version. `confirm_candidate`
requires a candidate-bound new surface state, a parent-state and candidate-
geometry binding, a declared current mechanical reanalysis, and energy
receipt/hash/unit fields with finite terms. It checks the serialized state
hash, retained geometry/ownership, candidate strip area, and `Gc * area`.
Missing or mismatched evidence leaves the old version active. The only saved
successful transition is explicitly `COMMITTED_FIXTURE_ONLY` from
`EXPLICIT_FIXTURE_NOT_SOLVER`; that state still says `NOT_COUPLED`. This is
not evidence of physical advance.

The controller records `external_work - delta_potential - Gc_delta_area` as a
residual but applies no balance tolerance or pass criterion. That is an
intentional limit documented by the source, so acceptance is limited to
receipt-gated state-machine wiring. Hashes, solver/status labels, and
`REVIEWED` declarations are caller-supplied and are not verified as authentic
receipts by this module. A production integration must establish that trust
boundary and independently review the energy method before allowing a real
surface version to become mechanically active. Only the declared straight
finite `+x` front is supported; turning and branching remain unsupported.

No new front-step tests were run during this independent review. The checked
receipt already matches the exact current source and fixtures; the source is
still explicitly experimental, uncalibrated, and `NOT_COUPLED`.
