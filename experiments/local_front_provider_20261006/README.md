# Pinned local-front MLS provider

This source-only experiment maps the pinned 100-source role manifest into the existing scalar `u_z` registry and continuous 3-D MLS shapes. It accepts only the exact support case, classifier, manifest, registry, MLS modules, and manufactured-field fixture whose hashes are embedded in `local_front_provider.py`.

The provider uses all 100 sources without a visibility mask: 100 BASE shapes, 30 shifted `F1/sqrt(R0)` shapes, and 10 shifted Heaviside shapes, for 140 canonical semantic IDs. The local support radius is `h=0.007 m`, `R0=0.01 m`; shifted sources must have nonzero `y`. Crack-face rows require an explicit PLUS or MINUS side at exact `y=0` points in the open finite patch. Front, ambiguous, finite-end, outside-domain, and MLS support-switch rows return empty `UNKNOWN` samples. It does not normalize enriched weights or assemble an operator.

Run the single bounded provider diagnostic with the configured Python runtime and `run_provider_diagnostic.py`. It caps each frozen batch at 256 points and writes `provider_diagnostic.json` plus `local_front_section.png`. The JSON includes actual section samples, trace rows, status records, results, and source pins. The section image is rendered directly from those saved samples.

The measured whole-box `z=0` section error is descriptive. Exact shifted-F1 reproduction is checked only at four sample points with `r<0.001 m`; this does not qualify the field in the collar or couple the provider to a mechanics solve.
