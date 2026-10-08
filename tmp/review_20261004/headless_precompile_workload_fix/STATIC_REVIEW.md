# Static review: headless precompile workload candidate

Status: static candidate only. No Julia process, package import, solver, pytest, or runtime precompile was run.

## Source of truth and baseline identity

The reviewed source is the pinned PGFrag resource at pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/auxiliary/precompile_workload.jl. Its SHA-256 is recorded in candidate_freeze_manifest.json. The pinned Peridynamics package is UUID 4dc47793-80f3-4232-b30e-ca78ca9d621b, version 0.5.4+pgfrag.1; the resource Manifest entry points to vendor/Peridynamics.jl. The baseline Project and Manifest hashes are also recorded in the freeze manifest.

The raw run01 stderr artifact at pgfrag/tmp/review_20261004/pd_point_event_probe/receipts/zero_load_control_run01.stderr.log (SHA-256: 212e034693d4035986be206b89e1c89de402bc2b8ee0f0a0ce326451d0853bee) says precompilation errored and names ArgumentError: Peridynamics headless backend rejects file-output Job options, with the include location at this workload line 1 and the catch location in src/Peridynamics.jl. That message remains an error; this candidate does not change the include try/catch or reinterpret its status.

## Cause and compatibility path

In src/core/job.jl:56-60, Job(spatial_setup, time_solver; kwargs...) collects the keywords and calls get_job_options. In the pinned headless implementation at src/core/job.jl:147-150, any nonempty keyword dictionary throws the observed ArgumentError, while an empty dictionary returns JobOptions(spatial_setup).

The no-export constructors in src/core/job.jl:138-145 set export_allowed=false for both a body and a multibody setup. src/auxiliary/logs.jl:26-35 returns before creating output folders or a logfile in that state; add_to_logfile at lines 46-52 also returns. log_it at lines 153-157 still uses the normal logger, and submit accepts quiet at src/core/submit.jl:1-18. With quiet=true, the workload retains the package's normal quiet submit route.

The VelocityVerlet reference and per-step export calls route through src/core/threads_body_data_handler.jl:195-209 and the multibody path in src/core/threads_multibody_data_handler.jl:105-112; those entry points return when export_allowed is false. Thus this workload does not reach the headless _export_results error at src/auxiliary/io.jl:124-126.

## Candidate delta and preservation

The candidate removes only root = mktempdir() and the four path variables that were used solely for VTK file output. It changes the four submits to Job(b1, vv), Job(b2, vv), Job(b3, vv), and Job(ms, vv), each retaining quiet=true.

The two point/volume arrays, all four body/material setups, numerical material parameters, point sets, precrack/no-failure setup, initial conditions, boundary conditions, contact setup, and VelocityVerlet(steps=2) remain byte-for-byte unchanged. The mpi_run() safety guard and the package-level catch remain unchanged. No constitutive model or submit implementation is modified.

Expected behavior from source inspection: all four constructors reach the no-export options path and avoid the rejecting nonempty-options branch. This is a source-based prediction only, not a runtime result. It does not establish nonzero force, material response, crack growth, or physical validity.