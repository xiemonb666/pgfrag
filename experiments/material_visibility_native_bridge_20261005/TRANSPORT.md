# Native mask transport

`native_mask_transport.py` compiles one strict JSON visibility ledger into a
stdlib TOML transport for a future, versioned native `InitialBondPolicy`
constructor. It is a bounded serialization and consistency layer. It does not
run Julia, construct a bond graph, authorize initialization, or qualify a
material asset.

The JSON file is the input artifact. The compiler records the SHA-256 of its
actual bytes, checks the JSON schema and the ledger's internal point-order,
pair-map, mask, count, and zero-degree relationships, then writes a TOML
projection. External point IDs are preserved and explicitly mapped to native
indices `1..N` in the exact point-row order. Candidate pairs are canonicalized
by external point ID, while `native_i` and `native_j` follow point-row order;
therefore `native_i` may be greater than `native_j`. A native reader must
canonicalize native indices before graph comparison and must retain the
external-ID/order binding. The output is written through a temporary file,
read back through Python's `tomllib`, and compared by value and type with the
projection built from the input. The compiler refuses an output that already
exists at preflight and checks again immediately before replacement. This is
not a cross-process no-clobber guarantee: another process can create the path
between the final check and replacement. The compiler rejects duplicate JSON keys,
non-standard JSON numbers, unsupported fields, and values outside its fixed
limits (64 MiB input, 128 MiB TOML output, 5,000 points, and 250,000 candidate
pairs). It does not create missing output directories.

The result is always `TRANSPORT_ONLY_NATIVE_UNVERIFIED`; the native gate stays
`UNVERIFIED`, and `initialization_authorized` stays false. A positive-degree
mask eligibility value describes only the supplied ledger's internal shape.
Unknown pairs, an empty candidate set, an empty included mask, and zero-degree
points remain blocking states. Planned initial physics-event fields remain
planned and are not observations from a simulation.

## Validation boundary

The compiler does not establish that the classifier source is correct, that
the source CAD or point payload at a referenced path still has the declared
bytes, or that the classifier result and raw-hit interval sidecar support the
ledger's decisions. It carries those paths and declared digests as metadata;
it does not open those referenced files. A matching sidecar digest or selected
pytest node/call index is not proof of physical truth or of a test having run.
Likewise, a recorded pytest exit code is an unsigned claim, not execution
attestation. The JSON and TOML digests bind bytes; they do not establish who
produced those bytes or whether their claims are scientifically correct.

The compiler also does not reconstruct the complete Euclidean horizon
candidate set, check that every supplied pair is within the horizon, or prove
that no eligible pair is missing. Its candidate list and horizon are declared
ledger content. Internal checksums and count agreement cannot establish
horizon completeness.

## Required native checks before initialization

The future pinned native fork must independently load and validate the actual
JSON and TOML bytes, bind the policy to the actual source CAD hash, point
payload hash, point-row order, coordinate frame, unit conversion, and horizon,
then rebuild the **complete** raw Euclidean candidate-pair set using the
pinned native convention. It must compare every canonical pair and its
visibility state against the policy. Missing, extra, reordered-inconsistently,
or differently classified pairs must block the whole asset. Comparing only
pair counts or a supplied checksum is insufficient. Any `UNKNOWN`, empty raw
candidate set, empty included mask, or zero-degree point must block whole-asset
initialization; filtering unknown rows and continuing is not qualification.

After that comparison, the native constructor must verify its actual graph
readback against the included material pairs and rebuild every graph-derived
quantity before initialization, including neighbor counts, weighted volumes,
dilatation support, surface correction, and any halo or kernel data. Geometry
exclusions are initial-domain exclusions, not physical fracture events. Later
physical bond failures must remain separately represented and counted. A
horizon cutoff is a candidate-generation rule, not a phase or fracture event.

Only a future native readback and these checks can change the native gate.
This Python transport cannot produce a native pass or a qualified asset.
