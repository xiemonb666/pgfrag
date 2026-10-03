using Peridynamics
using TOML
using SHA

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const JULIA_PROJECT = joinpath(REPO_ROOT, "src", "pgfrag", "resources", "julia")
const VENDORED_PD = joinpath(JULIA_PROJECT, "vendor", "Peridynamics.jl")
const DX = 2.0e-3
const HORIZON = 3.0e-3
const RHO = 1000.0
const YOUNG = 1.0e6
const GC = 1.0
const NU = 0.25
const CONTACT_RADIUS = 2.0e-3
const PENALTY = 1.0e6
const SPEED = 0.2
const N_PER_BODY = 8

function sha256_file(path)
    return bytes2hex(open(sha256, path))
end

function make_body(center_z, vz)
    position = Matrix{Float64}(undef, 3, N_PER_BODY)
    k = 0
    for z in (-DX / 2, DX / 2), y in (-DX / 2, DX / 2), x in (-DX / 2, DX / 2)
        k += 1
        position[:, k] .= (x, y, center_z + z)
    end
    body = Body(BBMaterial(), position, fill(DX^3, N_PER_BODY))
    material!(body; horizon=HORIZON, rho=RHO, E=YOUNG, Gc=GC, nu=NU)
    no_failure!(body)
    velocity_ic!(body, :all_points, :z, vz)
    return body, copy(position)
end

function write_points(path, dh, rests, initial_velocities)
    headers = ["point_label", "body_name", "local_point_id", "native_point_id", "export_point_id",
               "rest_x_m", "rest_y_m", "rest_z_m", "final_x_m", "final_y_m", "final_z_m",
               "displacement_x_m", "displacement_y_m", "displacement_z_m", "initial_vx_m_per_s",
               "initial_vy_m_per_s", "initial_vz_m_per_s", "final_vx_m_per_s", "final_vy_m_per_s",
               "final_vz_m_per_s", "mass_kg", "b_int_x_N_per_m3", "b_int_y_N_per_m3",
               "b_int_z_N_per_m3", "active_peridynamic_bonds", "broken_peridynamic_bonds"]
    open(path, "w") do io
        println(io, join(headers, ','))
        for body_idx in Peridynamics.each_body_idx(dh)
            body_name = dh.body_names[body_idx]
            body_dh = Peridynamics.get_body_dh(dh, body_idx)
            length(body_dh.chunks) == 1 || error("expected one native chunk per body")
            chunk = only(body_dh.chunks)
            system, storage = chunk.system, chunk.storage
            native_ids = system.chunk_handler.point_ids
            native_ids == collect(1:N_PER_BODY) || error("unexpected native point ordering for body $body_name")
            for local_id in 1:N_PER_BODY
                native_id = native_ids[local_id]
                export_id = (body_idx - 1) * N_PER_BODY + native_id
                label = string(body_name, native_id)
                active = storage.n_active_bonds[local_id]
                broken = length(Peridynamics.each_bond_idx(system, local_id)) - active
                row = Any[label, body_name, local_id, native_id, export_id,
                          rests[body_idx][1, native_id], rests[body_idx][2, native_id], rests[body_idx][3, native_id],
                          storage.position[1, local_id], storage.position[2, local_id], storage.position[3, local_id],
                          storage.displacement[1, local_id], storage.displacement[2, local_id], storage.displacement[3, local_id],
                          0.0, 0.0, initial_velocities[body_idx],
                          storage.velocity[1, local_id], storage.velocity[2, local_id], storage.velocity[3, local_id],
                          RHO * system.volume[local_id], storage.b_int[1, local_id], storage.b_int[2, local_id],
                          storage.b_int[3, local_id], active, broken]
                println(io, join(string.(row), ','))
            end
        end
    end
    return nothing
end

const ARG_OFFSET = length(ARGS) == 3 ? 1 : 0
length(ARGS) == ARG_OFFSET + 2 || error("usage: run_native_contact_minimal.jl CASE OUTPUT_DIR")
const CASE = ARGS[ARG_OFFSET + 1]
const OUTDIR = abspath(ARGS[ARG_OFFSET + 2])
const VALID_CASES = ("preflight", "contact_baseline", "contact_halfdt", "nocontact_control")
CASE in VALID_CASES || error("unknown case: $CASE")
mkpath(OUTDIR)

body_a, rest_a = make_body(2.25e-3, -SPEED)
body_b, rest_b = make_body(-2.25e-3, SPEED)
ms = MultibodySetup(:a => body_a, :b => body_b)
uses_contact = CASE != "nocontact_control"
uses_contact && contact!(ms, :a, :b; radius=CONTACT_RADIUS, penalty_factor=PENALTY)
length(ms.srf_contacts) == (uses_contact ? 2 : 0) || error("unexpected native contact registration")
all(!, body_a.fail_permit) && all(!, body_b.fail_permit) || error("no_failure! was not applied")

bound_solver = VelocityVerlet(steps=240, safety_factor=0.5)
bound_dh = Peridynamics.threads_data_handler(ms, bound_solver, 1)
native_bound = Peridynamics.calc_stable_timestep(bound_dh, 1.0)
isfinite(native_bound) && native_bound > 0 || error("native material stability bound is not positive finite")

steps = CASE == "preflight" ? 1 : CASE == "contact_halfdt" ? 480 : 240
dt = CASE == "contact_halfdt" ? 0.25 * native_bound : 0.5 * native_bound
initial_vz = [-SPEED, SPEED]
rests = [rest_a, rest_b]
solver = VelocityVerlet(steps=steps, stepsize=dt, safety_factor=0.5)
dh = submit(Job(ms, solver); quiet=true)
length(dh.body_dhs) == 2 || error("native submit did not return both bodies")
write_points(joinpath(OUTDIR, "points.csv"), dh, rests, initial_vz)

bond_counts = Dict{String,Any}()
max_abs_external = 0.0
for body_idx in Peridynamics.each_body_idx(dh)
    body_dh = Peridynamics.get_body_dh(dh, body_idx)
    chunk = only(body_dh.chunks)
    storage = chunk.storage
    bonds = storage.bond_active
    bond_counts[string(dh.body_names[body_idx])] = Dict(
        "active" => count(identity, bonds), "broken" => count(!, bonds), "directed_total" => length(bonds))
    global max_abs_external = max(max_abs_external, maximum(abs, storage.b_ext))
end

physics_path = joinpath(VENDORED_PD, "src", "physics", "short_range_force_contact.jl")
unit_path = joinpath(REPO_ROOT, "vendor", "Peridynamics.jl", "test", "unit", "physics", "test_short_range_force_contact.jl")
summary = Dict{String,Any}(
    "schema" => "native_contact_minimal.case.v1", "case" => CASE,
    "native_source_path" => abspath(@__FILE__), "native_case_argv" => copy(ARGS),
    "process_exit_code" => 0,
    "submit_completed" => true, "preflight_is_contact_success" => false,
    "native_submit_api" => "submit(Job(ms, VelocityVerlet(...)); quiet=true)",
    "native_data_handler" => "ThreadsMultibodyDataHandler; submit return; one chunk/body",
    "step_count" => solver.n_steps, "solver_dt_s" => solver.Δt,
    "end_time_s" => solver.end_time, "native_material_stability_bound_s" => native_bound,
    "stability_bound_multiplier" => 1.0, "point_spacing_m" => DX,
    "cube_edge_m" => 4.0e-3, "points_per_body" => N_PER_BODY, "total_points" => 2 * N_PER_BODY,
    "initial_center_a_z_m" => 2.25e-3, "initial_center_b_z_m" => -2.25e-3,
    "initial_min_point_distance_m" => 2.5e-3, "horizon_m" => HORIZON,
    "density_kg_per_m3" => RHO, "young_modulus_Pa" => YOUNG,
    "fracture_energy_J_per_m2" => GC, "poisson_ratio" => NU,
    "no_failure" => true, "contact_enabled" => uses_contact,
    "contact_radius_m" => uses_contact ? CONTACT_RADIUS : 0.0,
    "contact_penalty_factor" => uses_contact ? PENALTY : 0.0,
    "initial_speed_m_per_s" => SPEED, "boundary_condition_count" => 0,
    "external_force_max_abs_N_per_m3" => max_abs_external,
    "active_bond_counts_by_body" => bond_counts,
    "point_state_units" => Dict("rest_position" => "m", "final_position" => "m",
        "displacement" => "m", "velocity" => "m/s", "mass" => "kg", "native_final_b_int" => "N/m^3"),
    "source_hashes_sha256" => Dict(
        "native_contact_script" => sha256_file(@__FILE__),
        "vendored_contact_physics" => sha256_file(physics_path),
        "upstream_contact_unit_test" => sha256_file(unit_path),
        "julia_project" => sha256_file(joinpath(JULIA_PROJECT, "Project.toml")),
        "julia_manifest" => sha256_file(joinpath(JULIA_PROJECT, "Manifest.toml"))))
open(joinpath(OUTDIR, "summary.toml"), "w") do io
    TOML.print(io, summary)
end
