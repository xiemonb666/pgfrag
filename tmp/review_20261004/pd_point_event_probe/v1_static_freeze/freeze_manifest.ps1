$ErrorActionPreference = 'Stop'
$Workspace = 'C:\Users\xiemo\Desktop\ShardForge'
$Experiment = 'pgfrag/experiments/pd_point_event_probe_20261004'
$FileList = @(
    @{ path = "$Experiment/README.md"; role = 'experiment_contract' },
    @{ path = "$Experiment/run_pd_point_event_probe.jl"; role = 'native_event_observer' },
    @{ path = "$Experiment/make_point_inputs.ps1"; role = 'deterministic_fixture_generator' },
    @{ path = "$Experiment/freeze_manifest.ps1"; role = 'manifest_builder' },
    @{ path = "$Experiment/launch_zero_load_probe.ps1"; role = 'frozen_zero_load_launcher' },
    @{ path = "$Experiment/common_case.toml"; role = 'common_native_solver_and_material_config' },
    @{ path = "$Experiment/case_zero_load_control.toml"; role = 'first_runtime_gate_only' },
    @{ path = "$Experiment/case_tinyload_event_diagnostic.toml"; role = 'future_separately_authorized_case' },
    @{ path = "$Experiment/input/coupon_uniform_voxel_5cube_125.points.csv"; role = 'primary_frozen_points' },
    @{ path = "$Experiment/input/coupon_uniform_voxel_5cube_125.metadata.toml"; role = 'primary_point_contract' },
    @{ path = "$Experiment/input/coupon_cell_jitter_5cube_seed17.points.csv"; role = 'jitter_candidate_frozen_points' },
    @{ path = "$Experiment/input/coupon_cell_jitter_5cube_seed17.metadata.toml"; role = 'jitter_candidate_point_contract' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/Project.toml'; role = 'native_julia_project' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/Manifest.toml'; role = 'native_julia_lockfile' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/backend_provenance.json'; role = 'native_backend_provenance' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/run_fixture.jl'; role = 'frozen_baseline_runner_reference' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/Peridynamics.jl'; role = 'native_package_module' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/submit.jl'; role = 'native_submit_sequence' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/threads_body_data_handler.jl'; role = 'native_step_output_hooks' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/job.jl'; role = 'native_job_and_options_contract' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/data_handler.jl'; role = 'native_initial_force_call' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/systems.jl'; role = 'native_system_logger_path' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/time_solvers.jl'; role = 'native_generic_time_solver_logger' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/time_solvers/velocity_verlet.jl'; role = 'native_velocity_verlet_clock_and_step' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/discretization/body.jl'; role = 'native_arbitrary_point_body' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/discretization/point_sets.jl'; role = 'native_point_identity_sets' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/boundary_conditions.jl'; role = 'native_velocity_boundary_conditions' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/discretization/bond_system.jl'; role = 'native_bond_failure_dispatch' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/ordinary_state_based.jl'; role = 'native_osb_critical_stretch_failure' },
    @{ path = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/auxiliary/logs.jl'; role = 'native_options_logger_access' }
)

$FileHashes = @()
$HashByPath = @{}
foreach ($Entry in $FileList) {
    $AbsolutePath = Join-Path $Workspace ($Entry.path -replace '/', '\')
    $Hash = (Get-FileHash -LiteralPath $AbsolutePath -Algorithm SHA256).Hash.ToLowerInvariant()
    $FileHashes += [ordered]@{ path = $Entry.path; role = $Entry.role; sha256 = $Hash }
    $HashByPath[$Entry.path] = $Hash
}

$RuntimePath = 'C:\Users\xiemo\Desktop\ShardForge\third_party\pgfrag-runtime\extract\julia-1.10.12\bin\julia.exe'
$RuntimeHash = (Get-FileHash -LiteralPath $RuntimePath -Algorithm SHA256).Hash.ToLowerInvariant()
$NativeRoot = 'pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src'
$Manifest = [ordered]@{
    schema = 'pgfrag.pd.point_event_probe.manifest.v1'
    experiment_id = 'pd_point_event_probe_20261004'
    prepared_utc_date = '2026-10-04'
    status = 'SOURCE_PREPARED_STATIC_REVIEW_ONLY_AWAITING_ROOT_FINAL_FREEZE'
    scope = [ordered]@{
        baseline = 'Existing PD baseline stays frozen; this is one isolated native OSB point-event probe.'
        max_material_points = 512
        no_tetrahedral_crack_facets = $true
        no_solver_rewrite_or_restart_per_step = $true
        native_submit_used_once = $true
    }
    fixture = [ordered]@{
        primary_input_id = 'coupon_uniform_voxel_5cube_125_v1'
        primary_points = "$Experiment/input/coupon_uniform_voxel_5cube_125.points.csv"
        primary_metadata = "$Experiment/input/coupon_uniform_voxel_5cube_125.metadata.toml"
        point_count = 125
        declared_cube_extent_m = @(0.005, 0.005, 0.005)
        point_volume_sum_m3 = 1.25e-7
        volume_weighted_com_m = @(0.0025, 0.0025, 0.0025)
        primary_quadrature = 'Uniform voxel cell-center; each point has a known positive 1e-9 m^3 cell weight.'
        jitter_candidate_input_id = 'coupon_cell_jitter_5cube_seed17_v1'
        jitter_candidate_points = "$Experiment/input/coupon_cell_jitter_5cube_seed17.points.csv"
        jitter_semantics = 'Deterministic within-cell positional jitter with retained cell quadrature weights; not exact arbitrary Voronoi volume and not a real shell.'
        external_input_contract = 'Frozen row-order CSV point IDs plus SI xyz/positive volume and required TOML source, units, hash, point count, volume sum, COM, quadrature semantics, and optional driven point-ID sets.'
    }
    material = [ordered]@{
        model = 'native OSBMaterial / CriticalStretch'
        rho_kg_m3 = 2100.0
        E_pa = 5000000000.0
        nu = 0.24
        Gc_j_m2 = 45.0
        calibration_status = 'UNVALIDATED_PROTOTYPE'
    }
    native_observer = [ordered]@{
        runtime_route = 'Job(body, VelocityVerlet, EventProbeOptions); Peridynamics.submit(job; quiet=true) exactly once.'
        event_hooks = @(
            [ordered]@{
                signature = 'export_reference_results(dh::ThreadsBodyDataHandler, options::EventProbeOptions)'
                native_source = "$NativeRoot/time_solvers/velocity_verlet.jl"
                source_sha256 = $HashByPath["$NativeRoot/time_solvers/velocity_verlet.jl"]
                semantics = 'Native solve! invokes this before step 1; observer records initialized state at step 0 and t=0.'
            },
            [ordered]@{
                signature = 'export_results(dh::ThreadsBodyDataHandler, options::EventProbeOptions, chunk_id::Int, timestep::Int, time::Float64)'
                native_source = "$NativeRoot/core/threads_body_data_handler.jl"
                source_sha256 = $HashByPath["$NativeRoot/core/threads_body_data_handler.jl"]
                semantics = 'Native verlet_timestep! calls this after position, force, and velocity updates at every native step.'
            }
        )
        native_flow_signatures = @(
            [ordered]@{ signature = 'submit(job::Job; kwargs...)'; source = "$NativeRoot/core/submit.jl"; sha256 = $HashByPath["$NativeRoot/core/submit.jl"] },
            [ordered]@{ signature = '_submit_threads(job::Job, n_chunks::Int)'; source = "$NativeRoot/core/submit.jl"; sha256 = $HashByPath["$NativeRoot/core/submit.jl"] },
            [ordered]@{ signature = 'solve!(dh::AbstractDataHandler, vv::VelocityVerlet, options::AbstractJobOptions)'; source = "$NativeRoot/time_solvers/velocity_verlet.jl"; sha256 = $HashByPath["$NativeRoot/time_solvers/velocity_verlet.jl"] },
            [ordered]@{ signature = 'verlet_timestep!(dh::AbstractThreadsBodyDataHandler, options::AbstractJobOptions, Δt::Float64, Δt½::Float64, n::Int)'; source = "$NativeRoot/time_solvers/velocity_verlet.jl"; sha256 = $HashByPath["$NativeRoot/time_solvers/velocity_verlet.jl"] },
            [ordered]@{ signature = 'Body(mat::M, position::AbstractMatrix, volume::AbstractVector)'; source = "$NativeRoot/discretization/body.jl"; sha256 = $HashByPath["$NativeRoot/discretization/body.jl"] },
            [ordered]@{ signature = 'calc_failure!(storage::OSBStorage, system::BondSystem, ::OSBMaterial, ::CriticalStretch, paramsetup::AbstractParameterSetup, i)'; source = "$NativeRoot/physics/ordinary_state_based.jl"; sha256 = $HashByPath["$NativeRoot/physics/ordinary_state_based.jl"] }
        )
        output_fields = 'Point reference/current position, displacement, velocity, acceleration, b_int/b_ext force densities in N/m^3, damage, volume/mass; every directed bond state and canonical undirected pair state per sample; transitions are observed native state differences.'
        logger_compatibility = 'EventProbeOptions <: AbstractJobOptions and export_allowed=false. Native init_logs/add_to_logfile short-circuit before vtk/logfile fields; experiment defines Peridynamics.log_job_options(::EventProbeOptions) because pinned source specializes that logger to JobOptions; log_data_handler/system and native timesolver messages use log_it/add_to_logfile, safe with export_allowed=false. Hook specializations still capture despite export_allowed=false.'
        hook_contract_status = 'Internal native output generics, not a documented public user callback; exact dispatch must be proven by the first separately released no-damage preflight.'
        internal_init_or_step_override = $false
        event_clock = 'Step 0 after native initialize!; later samples follow native verlet_timestep! at n × accepted native Δt.'
        energy = 'UNKNOWN when no complete term and definition exist; absent quantities are never replaced with zero.'
        graph_interpretation = 'Native bond graph only; no physical fragment, crack-face, contact, R2, or training claim.'
    }
    stability = [ordered]@{
        threads = 1
        max_chunks = 1
        dt_mode = 'AUTO_NATIVE_STABLE_BOUND'
        retained_check = 'Compute native calc_stable_timestep on a single-chunk handler before submit, then require accepted solver.Δt to equal that bound within the frozen baseline tolerance.'
    }
    cases = @(
        [ordered]@{ case_id = 'zero_load_control'; config = "$Experiment/case_zero_load_control.toml"; steps = 1; runtime_status = 'FIRST_RUNTIME_GATE_PENDING_ROOT_FINAL_FREEZE'; zero_load = $true },
        [ordered]@{ case_id = 'tinyload_event_diagnostic'; config = "$Experiment/case_tinyload_event_diagnostic.toml"; steps = 5; speed_m_s = 0.0001; runtime_status = 'SEPARATE_RUNTIME_APPROVAL_REQUIRED'; zero_load = $false }
    )
    launcher = [ordered]@{
        script = "$Experiment/launch_zero_load_probe.ps1"
        fixed_case = 'zero_load_control'
        tinyload_reachable = $false
        records_numeric_julia_process_exit_code = $true
        timeout_seconds = 120
    }
    runtime = [ordered]@{
        executable = $RuntimePath
        sha256 = $RuntimeHash
        version_path_component = 'julia-1.10.12; version not runtime-probed under the static-only instruction'
        project_path = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\src\pgfrag\resources\julia'
        depot_path = 'C:\Users\xiemo\Desktop\ShardForge\third_party\pgfrag-runtime\depot'
        threads = 1
        precompile_tasks = 1
    }
    output_paths = [ordered]@{
        out_directory = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\tmp\review_20261004\pd_point_event_probe\OUT'
        out_must_be_empty_before_launch = $true
        snapshots_sibling = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\tmp\review_20261004\pd_point_event_probe\snapshots'
        receipts_sibling = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\tmp\review_20261004\pd_point_event_probe\receipts'
    }
    runtime_gate = [ordered]@{
        first_action = 'One invocation, <=512 points, one native zero-load step, no damage requirement.'
        timeout_seconds = 120
        status = 'NOT_AUTHORIZED_TO_RUN_UNTIL_ROOT_FINAL_FREEZE_REVIEW'
        acceptance = @(
            'two samples at steps 0 and 1',
            'one native chunk and one thread',
            'frozen point IDs complete and unique',
            'positive point mass closes against rho times point volume sum',
            'accepted native dt matches native stable bound',
            'complete reciprocal directed bond states agree with canonical states at both samples',
            'zero-load produces no bond-activity transition',
            'OUT remains empty; snapshots and receipt stay in sibling directories'
        )
    }
    later_gates = @(
        [ordered]@{ gate = 'tinyload_event_diagnostic'; status = 'PENDING_SEPARATE_RUNTIME_APPROVAL' },
        [ordered]@{ gate = 'stock_native_zero_load_terminal_comparison'; status = 'OPTIONAL_FUTURE_GATE_PENDING_REVIEW' }
    )
    file_hashes = $FileHashes
}

$ManifestPath = Join-Path $Workspace 'pgfrag/experiments/pd_point_event_probe_20261004/manifest.json'
$Json = ($Manifest | ConvertTo-Json -Depth 16) + [Environment]::NewLine
$Utf8NoBom = New-Object System.Text.UTF8Encoding($false)
[System.IO.File]::WriteAllText($ManifestPath, $Json, $Utf8NoBom)
Write-Output "Wrote manifest with $($FileHashes.Count) pinned files and runtime SHA-256 $RuntimeHash"
