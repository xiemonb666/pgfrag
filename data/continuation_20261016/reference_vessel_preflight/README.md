# Reference vessel zero-step preflight

The native preflight preserved its inputs and stopped before Akantu initialization. The frozen v5 mesh declares `cut_interface` and `source_exterior`; the unchanged wrapper requires physical groups named `interface` and, for load-patch construction, `exterior`. The run returned `FAILED_PRESERVED_NO_RETRY` with process exit code 2. No geometry, mesh labels, material values, or scientific parameters were changed, and no retry or solver evolution occurred.

The source reports 56,139 T4 elements, 18,236 nodes, and 541 interface Tri3 faces, giving 1,623 native quadrature points under the wrapper's three-point-per-face contract. The reported exterior source-face ancestry remains unbound; this preflight grants no fragment or surface-paint qualification.

A total 50 μs history forecast is unavailable because native insertion did not produce a stable time step or a feasible balanced patch. The wrapper's history arrays require at least 181,776 bytes per time frame for the 1,623 QPs plus 320 bytes for the energy, work, and ledger arrays, before patch displacement or saved node states. Each saved node-state frame adds 1,750,656 bytes. With the wrapper's 2,000,000-byte safety margin, the 64 MiB cap admits at most 357 time frames from QP and time arrays alone; actual capacity is lower after patch and node-state arrays.

A rerun would require a main-agent decision on a byte-audited physical-name alias copy (`cut_interface` → `interface`, `source_exterior` → `exterior`). That operation would change only names in `$PhysicalNames`, but this run did not make that change.

See `execution_record.json` for the exact argv, hashes, exit code, logs, runner receipt, and forecast formula. Its SHA-256 is stored in `execution_record.sha256`.