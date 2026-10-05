using Peridynamics
using LinearAlgebra
using SHA
using TOML

function require(condition, message)
    condition || error(message)
end

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const JULIA_PROJECT = joinpath(REPO_ROOT, "pgfrag", "src", "pgfrag", "resources", "julia")
const SOURCE_FREEZE_PATH = joinpath(@__DIR__, "source_freeze.toml")

sha256_file(path::String) = bytes2hex(sha256(read(path)))
arithmetic_mean(values) = sum(values) / length(values)

function load_and_verify_source_freeze()
    isfile(SOURCE_FREEZE_PATH) || error("source_freeze.toml is missing")
    freeze = TOML.parsefile(SOURCE_FREEZE_PATH)
    require(get(freeze, "schema", "") == "pgfrag.free_evolution_coupon.source_freeze.v1",
            "unsupported source-freeze schema")
    require(get(freeze, "freeze_status", "") == "FROZEN_SOURCE_ONLY",
            "source freeze is not marked FROZEN_SOURCE_ONLY")
    pins = get(freeze, "pins", Any[])
    require(!isempty(pins), "source freeze contains no pinned files")
    actual_hashes = Dict{String,String}()
    for pin in pins
        relative_path = String(pin["path"])
        path = normpath(joinpath(REPO_ROOT, split(relative_path, '/')...))
        require(isfile(path), "pinned source is missing: $relative_path")
        actual_size = filesize(path)
        actual_hash = sha256_file(path)
        require(actual_size == Int(pin["bytes"]),
                "pinned source size changed: $relative_path")
        require(actual_hash == lowercase(String(pin["sha256"])),
                "pinned source SHA-256 changed: $relative_path")
        actual_hashes[relative_path] = actual_hash
    end
    return (; freeze, actual_hashes, freeze_sha256=sha256_file(SOURCE_FREEZE_PATH))
end

mutable struct CouponRecorder
    displacement_m::Array{Float64,3} # point × xyz × captured state
    velocity_m_s::Array{Float64,3}
    velocity_half_m_s::Array{Float64,3}
    acceleration_m_s2::Array{Float64,3}
    b_int_N_m3::Array{Float64,3}
    b_ext_N_m3::Array{Float64,3}
    strain_energy_density_J_m3::Matrix{Float64} # point × captured state
    damage::Matrix{Float64}
    active_bonds_per_point::Matrix{Int}
    seen::BitVector
    last_step::Int
    directed_bond_count::Int
    active_directed_bonds_at_reference::Int
end

function CouponRecorder(n_points::Int, n_steps::Int)
    n_states = n_steps + 1
    return CouponRecorder(zeros(n_points, 3, n_states), zeros(n_points, 3, n_states),
                          zeros(n_points, 3, n_states), zeros(n_points, 3, n_states),
                          zeros(n_points, 3, n_states), zeros(n_points, 3, n_states),
                          zeros(n_points, n_states), zeros(n_points, n_states),
                          zeros(Int, n_points, n_states), falses(n_states), -1, 0, 0)
end

mutable struct CouponOptions <: Peridynamics.AbstractJobOptions
    export_allowed::Bool
    recorder::CouponRecorder
    rest_position_m::Matrix{Float64}
    point_volume_m3::Vector{Float64}
    density_kg_m3::Float64
    solver::Peridynamics.VelocityVerlet
end

# This job writes only the explicit PG-Frag CSV/TOML files below; it does not use VTK/log export.
Peridynamics.log_job_options(::CouponOptions) = nothing

function collect_native_state!(options::CouponOptions,
                               dh::Peridynamics.ThreadsBodyDataHandler,
                               step::Int, time_s::Float64)
    recorder = options.recorder
    n_points = size(options.rest_position_m, 2)
    n_steps = options.solver.n_steps
    require(length(dh.chunks) == 1, "coupon recorder requires exactly one native chunk")
    require(step == recorder.last_step + 1,
            "native hook states are missing, duplicated, or out of sequence")
    require(0 <= step <= n_steps, "native hook returned a step outside the requested range")
    dt_s = options.solver.Δt
    require(isfinite(dt_s) && dt_s > 0, "native solver did not retain a positive time step")
    require(isapprox(time_s, step * dt_s; rtol=1e-12, atol=128 * eps(Float64) * dt_s),
            "native hook time disagrees with step × accepted dt")

    chunk = only(dh.chunks)
    system = chunk.system
    storage = chunk.storage
    native_point_ids = system.chunk_handler.point_ids
    require(length(native_point_ids) == n_points &&
            length(unique(native_point_ids)) == n_points &&
            Set(native_point_ids) == Set(1:n_points),
            "native point indices do not map one-to-one to the analytic-box point IDs")
    require(system.chunk_handler.n_loc_points == n_points,
            "single chunk contains a halo or omits a specimen point")

    state_index = step + 1
    for local_i in 1:n_points
        point_id = native_point_ids[local_i]
        recorder.displacement_m[point_id, :, state_index] .= storage.displacement[:, local_i]
        recorder.velocity_m_s[point_id, :, state_index] .= storage.velocity[:, local_i]
        recorder.velocity_half_m_s[point_id, :, state_index] .= storage.velocity_half[:, local_i]
        recorder.acceleration_m_s2[point_id, :, state_index] .= storage.acceleration[:, local_i]
        recorder.b_int_N_m3[point_id, :, state_index] .= storage.b_int[:, local_i]
        recorder.b_ext_N_m3[point_id, :, state_index] .= storage.b_ext[:, local_i]
        recorder.damage[point_id, state_index] = storage.damage[local_i]
        recorder.active_bonds_per_point[point_id, state_index] = storage.n_active_bonds[local_i]
    end

    # This is the pinned native OSB field exporter; integration by point volumes happens below.
    strain_energy_density = Peridynamics.export_field(
        Val(:strain_energy_density), chunk.mat, system, storage, chunk.paramsetup, time_s)
    require(length(strain_energy_density) == n_points && all(isfinite, strain_energy_density),
            "native OSB strain-energy-density field is missing or non-finite")
    for local_i in 1:n_points
        point_id = native_point_ids[local_i]
        recorder.strain_energy_density_J_m3[point_id, state_index] =
            strain_energy_density[local_i]
    end

    n_directed = length(storage.bond_active)
    n_active = count(identity, storage.bond_active)
    n_broken = n_directed - n_active
    require(n_broken == 0,
            "failure-disabled elastic coupon unexpectedly contains an inactive bond")
    if step == 0
        recorder.directed_bond_count = n_directed
        recorder.active_directed_bonds_at_reference = n_active
    else
        require(n_directed == recorder.directed_bond_count &&
                n_active == recorder.active_directed_bonds_at_reference,
                "native bond count/activity changed during the no-failure run")
    end
    recorder.seen[state_index] = true
    recorder.last_step = step
    return nothing
end

# These exact options overloads observe native state zero and each completed VelocityVerlet state.
function Peridynamics.export_reference_results(dh::Peridynamics.ThreadsBodyDataHandler,
                                                options::CouponOptions)
    collect_native_state!(options, dh, 0, 0.0)
    return nothing
end

function Peridynamics.export_results(dh::Peridynamics.ThreadsBodyDataHandler,
                                     options::CouponOptions, chunk_id::Int,
                                     timestep::Int, time_s::Float64)
    require(chunk_id == 1, "single-chunk recorder received an unexpected chunk ID")
    collect_native_state!(options, dh, timestep, time_s)
    return nothing
end

function target_displacement_m(time_s::Float64, ramp_duration_s::Float64,
                               hold_duration_s::Float64, unload_duration_s::Float64,
                               target_displacement::Float64)
    require(isfinite(time_s), "load schedule time must be finite")
    t = max(time_s, 0.0)
    if t <= ramp_duration_s
        ξ = t / ramp_duration_s
        return target_displacement * (3 * ξ^2 - 2 * ξ^3)
    elseif t <= ramp_duration_s + hold_duration_s
        return target_displacement
    elseif t <= ramp_duration_s + hold_duration_s + unload_duration_s
        ξ = (t - ramp_duration_s - hold_duration_s) / unload_duration_s
        return target_displacement * (1 - 3 * ξ^2 + 2 * ξ^3)
    end
    return 0.0
end

function prescribed_step_velocity_m_s(time_s::Float64, dt_s::Float64,
                                      ramp_duration_s::Float64,
                                      hold_duration_s::Float64,
                                      unload_duration_s::Float64,
                                      target_displacement::Float64)
    require(isfinite(time_s) && time_s >= 0 && isfinite(dt_s) && dt_s > 0,
            "native boundary-condition time and time step must be finite and nonnegative/positive")
    current_u = target_displacement_m(time_s, ramp_duration_s, hold_duration_s,
                                      unload_duration_s, target_displacement)
    previous_u = target_displacement_m(time_s - dt_s, ramp_duration_s,
                                       hold_duration_s, unload_duration_s,
                                       target_displacement)
    # Velocity BC controls velocity_half before the native position update. This per-step
    # secant makes actual grip displacement match one continuous physical-time path at t_n.
    return (current_u - previous_u) / dt_s
end

const POINT_STATE_FIELDS = [
    "ux_m", "uy_m", "uz_m",
    "vx_m_s", "vy_m_s", "vz_m_s",
    "vhalf_x_m_s", "vhalf_y_m_s", "vhalf_z_m_s",
    "accel_x_m_s2", "accel_y_m_s2", "accel_z_m_s2",
    "b_int_x_N_m3", "b_int_y_N_m3", "b_int_z_N_m3",
    "b_ext_x_N_m3", "b_ext_y_N_m3", "b_ext_z_N_m3",
    "osb_strain_energy_density_J_m3", "damage", "n_active_bonds_exact_float64"]

const POINT_STATE_UNITS = [
    "m", "m", "m", "m/s", "m/s", "m/s", "m/s", "m/s", "m/s",
    "m/s^2", "m/s^2", "m/s^2", "N/m^3", "N/m^3", "N/m^3",
    "N/m^3", "N/m^3", "N/m^3", "J/m^3", "1", "count"]

function write_point_state_binary(path::String, layout_path::String,
                                 recorder::CouponRecorder, dt_s::Float64,
                                 max_output_bytes::Int64)
    n_points, _, n_states = size(recorder.displacement_m)
    n_fields = length(POINT_STATE_FIELDS)
    expected_bytes = Int64(sizeof(Float64)) * n_fields * n_points * n_states
    require(expected_bytes <= max_output_bytes,
            "full point-state binary exceeds the fixed output-byte budget")
    open(path, "w") do io
        frame = Matrix{Float64}(undef, n_fields, n_points)
        for state_index in 1:n_states
            for dim in 1:3
                frame[dim, :] .= @view recorder.displacement_m[:, dim, state_index]
                frame[3 + dim, :] .= @view recorder.velocity_m_s[:, dim, state_index]
                frame[6 + dim, :] .= @view recorder.velocity_half_m_s[:, dim, state_index]
                frame[9 + dim, :] .= @view recorder.acceleration_m_s2[:, dim, state_index]
                frame[12 + dim, :] .= @view recorder.b_int_N_m3[:, dim, state_index]
                frame[15 + dim, :] .= @view recorder.b_ext_N_m3[:, dim, state_index]
            end
            frame[19, :] .= @view recorder.strain_energy_density_J_m3[:, state_index]
            frame[20, :] .= @view recorder.damage[:, state_index]
            for point_id in 1:n_points
                active_count = recorder.active_bonds_per_point[point_id, state_index]
                active_count_float = Float64(active_count)
                require(Int(active_count_float) == active_count,
                        "active-bond count cannot be represented exactly as Float64")
                frame[21, point_id] = active_count_float
            end
            little_endian_words = htol.(reinterpret(UInt64, vec(frame)))
            write(io, reinterpret(UInt8, little_endian_words))
        end
    end
    actual_bytes = filesize(path)
    require(actual_bytes == expected_bytes,
            "point-state binary byte count disagrees with the pinned tensor shape")
    binary_sha256 = sha256_file(path)
    time_levels = Dict{String,String}(
        "displacement" => "at t_n",
        "full_velocity" => "at t_n",
        "velocity_half" => "for n>=1 centered at t_(n-1/2), after BC evaluated at t_n; n=0 is the initialized zero placeholder",
        "acceleration" => "at t_n after native force update",
        "internal_force_density" => "at t_n",
        "external_force_density" => "at t_n",
        "osb_strain_energy_density" => "at t_n from native export_field",
        "damage_and_active_bond_count" => "at t_n; failure is disabled in this case")
    layout = Dict{String,Any}(
        "schema" => "pgfrag.free_evolution_coupon.point_state_binary.v1",
        "data_file" => basename(path),
        "format" => "IEEE-754 Float64, little-endian",
        "storage_order" => "Julia column-major tensor [field_index, point_id, state_index]; field index varies fastest, then ascending point_id, then ascending state_index",
        "array_shape" => [n_fields, n_points, n_states],
        "field_names" => POINT_STATE_FIELDS,
        "field_units" => POINT_STATE_UNITS,
        "point_id_mapping" => "point_id is 1..point_count; analytic reference coordinates and volumes are in point_reference.csv",
        "state_index_mapping" => "state_index = step + 1; step 0 is the reference state, step 1..n_steps are native VelocityVerlet output states",
        "time_step_s" => dt_s,
        "state_time_levels" => time_levels,
        "active_bond_count_encoding" => "Float64 with exact integer representability checked before writing",
        "point_count" => n_points,
        "state_count" => n_states,
        "field_count" => n_fields,
        "byte_count" => actual_bytes,
        "sha256" => binary_sha256)
    write_summary(layout_path, layout)
    return (; byte_count=actual_bytes, sha256=binary_sha256,
            layout_path=basename(layout_path))
end

function write_point_reference_csv(path::String, recorder::CouponRecorder,
                                   rest_position_m::Matrix{Float64},
                                   volumes_m3::Vector{Float64},
                                   density_kg_m3::Float64,
                                   complete_horizon_mask::BitVector,
                                   horizon_margin_m::Vector{Float64},
                                   left_grip_ids::Vector{Int},
                                   right_grip_ids::Vector{Int})
    n_points = size(rest_position_m, 2)
    left = falses(n_points)
    right = falses(n_points)
    left[left_grip_ids] .= true
    right[right_grip_ids] .= true
    open(path, "w") do io
        println(io, "point_id,x_ref_m,y_ref_m,z_ref_m,volume_m3,mass_kg,complete_spherical_horizon,minimum_face_clearance_minus_horizon_m,left_grip,right_grip,reference_active_bonds")
        for point_id in 1:n_points
            row = (point_id, rest_position_m[1, point_id], rest_position_m[2, point_id],
                   rest_position_m[3, point_id], volumes_m3[point_id],
                   density_kg_m3 * volumes_m3[point_id], complete_horizon_mask[point_id],
                   horizon_margin_m[point_id], left[point_id], right[point_id],
                   recorder.active_bonds_per_point[point_id, 1])
            println(io, join(string.(row), ','))
        end
    end
    return nothing
end

function output_directory_bytes(directory::String)
    total = Int64(0)
    for (root, _, files) in walkdir(directory), file in files
        total += filesize(joinpath(root, file))
    end
    return total
end

function require_output_budget(directory::String, max_output_bytes::Int64)
    actual_bytes = output_directory_bytes(directory)
    require(actual_bytes <= max_output_bytes,
            "coupon output directory exceeded the fixed byte budget")
    return actual_bytes
end

function write_energy_history(path::String, recorder::CouponRecorder,
                              rest_position_m::Matrix{Float64}, volumes_m3::Vector{Float64},
                              left_grip_ids::Vector{Int}, right_grip_ids::Vector{Int},
                              density_kg_m3::Float64, dt_s::Float64)
    n_states = size(recorder.displacement_m, 3)
    mass_kg = density_kg_m3 .* volumes_m3
    total_mass_kg = sum(mass_kg)
    impulses_left_N_s = zeros(3, n_states)
    impulses_right_N_s = zeros(3, n_states)
    reaction_force_average_left_N = zeros(3, n_states)
    reaction_force_average_right_N = zeros(3, n_states)
    discrete_bc_kinetic_work_left_J = zeros(n_states)
    discrete_bc_kinetic_work_right_J = zeros(n_states)
    max_free_dof_impulse_residual_N_s = zeros(n_states)
    momentum_update_residual_N_s = zeros(3, n_states)
    actual_left_grip_u_x = zeros(n_states)
    actual_right_grip_u_x = zeros(n_states)
    kinetic_energy_J = zeros(n_states)
    osb_strain_energy_J = zeros(n_states)
    momentum_kg_m_s = zeros(3, n_states)
    native_net_force_N = zeros(3, n_states)
    is_left_grip = falses(length(volumes_m3))
    is_right_grip = falses(length(volumes_m3))
    is_constrained = falses(length(volumes_m3), 3)
    is_left_grip[left_grip_ids] .= true
    is_right_grip[right_grip_ids] .= true
    is_constrained[left_grip_ids, :] .= true
    is_constrained[right_grip_ids, :] .= true

    for state_index in 1:n_states
        actual_left_grip_u_x[state_index] = arithmetic_mean(
            recorder.displacement_m[left_grip_ids, 1, state_index])
        actual_right_grip_u_x[state_index] = arithmetic_mean(
            recorder.displacement_m[right_grip_ids, 1, state_index])
        osb_strain_energy_J[state_index] = sum(
            recorder.strain_energy_density_J_m3[:, state_index] .* volumes_m3)
        for point_id in eachindex(mass_kg)
            velocity = recorder.velocity_m_s[point_id, :, state_index]
            kinetic_energy_J[state_index] += 0.5 * mass_kg[point_id] * sum(abs2, velocity)
            momentum_kg_m_s[:, state_index] .+= mass_kg[point_id] .* velocity
            native_net_force_N[:, state_index] .+=
                (recorder.b_int_N_m3[point_id, :, state_index] .+
                 recorder.b_ext_N_m3[point_id, :, state_index]) .* volumes_m3[point_id]
        end
    end

    discrete_bc_kinetic_work_step_J = zeros(n_states)
    for state_index in 2:n_states
        previous_state = state_index - 1
        v_before_half = recorder.velocity_m_s[:, :, previous_state] .+
                         recorder.acceleration_m_s2[:, :, previous_state] .* (dt_s / 2)
        v_after_half = recorder.velocity_half_m_s[:, :, state_index]
        max_free_impulse = 0.0
        for point_id in eachindex(mass_kg), dim in 1:3
            Δv = v_after_half[point_id, dim] - v_before_half[point_id, dim]
            impulse = mass_kg[point_id] * Δv
            if is_constrained[point_id, dim]
                if is_left_grip[point_id]
                    impulses_left_N_s[dim, state_index] += impulse
                    discrete_bc_kinetic_work_left_J[state_index] +=
                        0.5 * mass_kg[point_id] *
                        (v_after_half[point_id, dim]^2 - v_before_half[point_id, dim]^2)
                elseif is_right_grip[point_id]
                    impulses_right_N_s[dim, state_index] += impulse
                    discrete_bc_kinetic_work_right_J[state_index] +=
                        0.5 * mass_kg[point_id] *
                        (v_after_half[point_id, dim]^2 - v_before_half[point_id, dim]^2)
                end
            else
                max_free_impulse = max(max_free_impulse, abs(impulse))
            end
        end
        reaction_force_average_left_N[:, state_index] .= impulses_left_N_s[:, state_index] ./ dt_s
        reaction_force_average_right_N[:, state_index] .= impulses_right_N_s[:, state_index] ./ dt_s
        max_free_dof_impulse_residual_N_s[state_index] = max_free_impulse
        discrete_bc_kinetic_work_step_J[state_index] =
            discrete_bc_kinetic_work_left_J[state_index] +
            discrete_bc_kinetic_work_right_J[state_index]

        impulse_scale = max(maximum(abs, impulses_left_N_s[:, state_index]),
                            maximum(abs, impulses_right_N_s[:, state_index]), 1e-30)
        free_tolerance = 256 * eps(Float64) * impulse_scale + 1e-30
        require(max_free_impulse <= free_tolerance,
                "unconstrained degree of freedom received a non-zero native half-step BC impulse")

        force_impulse_N_s = (dt_s / 2) .* (
            native_net_force_N[:, previous_state] .+ native_net_force_N[:, state_index])
        total_bc_impulse_N_s = impulses_left_N_s[:, state_index] .+
                               impulses_right_N_s[:, state_index]
        momentum_update_residual_N_s[:, state_index] =
            momentum_kg_m_s[:, state_index] .- momentum_kg_m_s[:, previous_state] .-
            force_impulse_N_s .- total_bc_impulse_N_s
        momentum_scale = max(maximum(abs, momentum_kg_m_s[:, state_index]),
                             maximum(abs, momentum_kg_m_s[:, previous_state]),
                             maximum(abs, force_impulse_N_s),
                             maximum(abs, total_bc_impulse_N_s), 1e-30)
        momentum_tolerance = 1024 * eps(Float64) * momentum_scale + 1e-30
        require(maximum(abs, momentum_update_residual_N_s[:, state_index]) <= momentum_tolerance,
                "native full-step momentum change disagrees with force and BC-impulse time levels")
    end

    actuator_work_trapezoid_surrogate_J = zeros(n_states)
    actuator_work_trapezoid_surrogate_step_J = zeros(n_states)
    for state_index in 2:n_states
        force_before_N = reaction_force_average_right_N[1, state_index - 1]
        force_after_N = reaction_force_average_right_N[1, state_index]
        delta_u_m = actual_right_grip_u_x[state_index] - actual_right_grip_u_x[state_index - 1]
        actuator_work_trapezoid_surrogate_step_J[state_index] =
            0.5 * (force_before_N + force_after_N) * delta_u_m
        actuator_work_trapezoid_surrogate_J[state_index] =
            actuator_work_trapezoid_surrogate_J[state_index - 1] +
            actuator_work_trapezoid_surrogate_step_J[state_index]
    end
    cumulative_discrete_bc_kinetic_work_J = cumsum(discrete_bc_kinetic_work_step_J)
    initial_total_energy_J = osb_strain_energy_J[1] + kinetic_energy_J[1]
    energy_balance_residual_diagnostic_J = osb_strain_energy_J .+ kinetic_energy_J .-
                                            initial_total_energy_J .-
                                            cumulative_discrete_bc_kinetic_work_J
    open(path, "w") do io
        println(io, "step,time_s,bc_impulse_interval_start_s,bc_impulse_interval_end_s,mass_kg,momentum_x_kg_m_s,momentum_y_kg_m_s,momentum_z_kg_m_s,native_net_force_x_N,native_net_force_y_N,native_net_force_z_N,osb_strain_energy_J,kinetic_energy_J,left_grip_u_x_mean_m,right_grip_u_x_mean_m,left_bc_impulse_x_N_s,left_bc_impulse_y_N_s,left_bc_impulse_z_N_s,left_bc_force_average_x_N,left_bc_force_average_y_N,left_bc_force_average_z_N,right_bc_impulse_x_N_s,right_bc_impulse_y_N_s,right_bc_impulse_z_N_s,right_bc_force_average_x_N,right_bc_force_average_y_N,right_bc_force_average_z_N,left_bc_discrete_kinetic_work_J,right_bc_discrete_kinetic_work_J,total_bc_discrete_kinetic_work_step_J,actuator_work_trapezoid_surrogate_step_J,actuator_work_trapezoid_surrogate_cumulative_J,max_free_dof_impulse_residual_N_s,momentum_update_residual_x_N_s,momentum_update_residual_y_N_s,momentum_update_residual_z_N_s,energy_balance_residual_diagnostic_J")
        for state_index in 1:n_states
            step = state_index - 1
            row = (step, step * dt_s, max(step - 1, 0) * dt_s, step * dt_s,
                   total_mass_kg,
                   momentum_kg_m_s[1, state_index], momentum_kg_m_s[2, state_index],
                   momentum_kg_m_s[3, state_index], native_net_force_N[1, state_index],
                   native_net_force_N[2, state_index], native_net_force_N[3, state_index],
                   osb_strain_energy_J[state_index], kinetic_energy_J[state_index],
                   actual_left_grip_u_x[state_index], actual_right_grip_u_x[state_index],
                   impulses_left_N_s[1, state_index], impulses_left_N_s[2, state_index],
                   impulses_left_N_s[3, state_index], reaction_force_average_left_N[1, state_index],
                   reaction_force_average_left_N[2, state_index], reaction_force_average_left_N[3, state_index],
                   impulses_right_N_s[1, state_index], impulses_right_N_s[2, state_index],
                   impulses_right_N_s[3, state_index], reaction_force_average_right_N[1, state_index],
                   reaction_force_average_right_N[2, state_index], reaction_force_average_right_N[3, state_index],
                   discrete_bc_kinetic_work_left_J[state_index],
                   discrete_bc_kinetic_work_right_J[state_index],
                   discrete_bc_kinetic_work_step_J[state_index],
                   actuator_work_trapezoid_surrogate_step_J[state_index],
                   actuator_work_trapezoid_surrogate_J[state_index],
                   max_free_dof_impulse_residual_N_s[state_index],
                   momentum_update_residual_N_s[1, state_index],
                   momentum_update_residual_N_s[2, state_index],
                   momentum_update_residual_N_s[3, state_index],
                   energy_balance_residual_diagnostic_J[state_index])
            println(io, join(string.(row), ','))
        end
    end
    return (; total_mass_kg, initial_total_energy_J, final_osb_strain_energy_J=osb_strain_energy_J[end],
            final_kinetic_energy_J=kinetic_energy_J[end],
            final_cumulative_discrete_bc_kinetic_work_J=cumulative_discrete_bc_kinetic_work_J[end],
            final_actuator_work_trapezoid_surrogate_J=actuator_work_trapezoid_surrogate_J[end],
            final_energy_balance_residual_diagnostic_J=energy_balance_residual_diagnostic_J[end],
            initial_momentum_kg_m_s=collect(momentum_kg_m_s[:, 1]),
            final_momentum_kg_m_s=collect(momentum_kg_m_s[:, end]),
            final_left_grip_u_x_mean_m=actual_left_grip_u_x[end],
            final_right_grip_u_x_mean_m=actual_right_grip_u_x[end],
            max_abs_right_grip_u_x_m=maximum(abs, actual_right_grip_u_x),
            max_free_dof_impulse_residual_N_s=maximum(max_free_dof_impulse_residual_N_s),
            max_momentum_update_residual_N_s=maximum(abs, momentum_update_residual_N_s))
end

function write_summary(path::String, values::Dict{String,Any})
    open(path, "w") do io
        TOML.print(io, values; sorted=true)
    end
    return nothing
end

function run_coupon(args::Vector{String})
    require(length(args) == 2, "expected arguments: pinned_case.toml new_output_directory")
    case_path = abspath(args[1])
    output_dir = abspath(args[2])
    require(isfile(case_path), "fixed coupon case file is missing")
    require(!ispath(output_dir), "output directory already exists; refusing to overwrite results")
    require(VERSION == v"1.10.12", "coupon source is pinned to Julia 1.10.12")
    active_project = Base.active_project()
    require(active_project !== nothing &&
            normpath(active_project) == normpath(joinpath(JULIA_PROJECT, "Project.toml")),
            "run under the pinned PG-Frag Julia project")
    freeze_record = load_and_verify_source_freeze()
    case = TOML.parsefile(case_path)
    require(get(case, "schema", "") == "pgfrag.free_evolution_coupon.case.v1",
            "unsupported coupon case schema")
    require(String(case["run_status"]) == "SOURCE_ONLY_NOT_RUN",
            "fixed case status must be SOURCE_ONLY_NOT_RUN")
    case_relative_path = replace(relpath(case_path, REPO_ROOT), '\\' => '/')
    require(haskey(freeze_record.actual_hashes, case_relative_path),
            "case file must be in the repository and explicitly pinned by source_freeze.toml")

    geometry = case["geometry"]
    material_config = case["material"]
    loading = case["loading"]
    solver_config = case["solver"]
    lengths_m = Float64[geometry["length_x_m"], geometry["length_y_m"], geometry["length_z_m"]]
    h_m = Float64(geometry["point_spacing_h_m"])
    horizon_m = Float64(geometry["horizon_delta_m"])
    n_expected = Int(geometry["expected_point_count"])
    expected_grid = Int.(geometry["expected_grid_shape"])
    expected_interior = Int(geometry["expected_complete_horizon_interior_points"])
    max_points = Int(case["limits"]["max_points"])
    max_steps = Int(case["limits"]["max_steps"])
    min_steps_per_phase = Int(case["limits"]["min_steps_per_phase"])
    max_output_bytes = Int64(case["limits"]["max_output_bytes"])
    require(0 < max_output_bytes <= 200 * 1024^2,
            "fixed coupon output budget may not exceed 200 MiB")
    require(String(geometry["shape"]) == "analytic_box", "geometry must remain analytic_box")
    require(0 < n_expected <= max_points <= 5000, "fixed analytic coupon exceeds the 5000-point budget")
    require(h_m > 0 && horizon_m > h_m && all(>(0), lengths_m),
            "invalid independent spacing, horizon, or analytic-box dimensions")
    require(length(expected_grid) == 3 &&
            all(isapprox.(lengths_m ./ h_m, expected_grid; rtol=0.0, atol=1e-10)),
            "analytic-box dimensions must divide exactly into the frozen point spacing")
    require(Float64(geometry["horizon_to_spacing_ratio"]) ≈ horizon_m / h_m,
            "declared δ/h ratio disagrees with independently stored h and δ")
    require(String(material_config["model"]) == "OSB" &&
            String(material_config["surface_correction"]) == "NoCorrection",
            "fixed case must remain native OSB with explicit NoCorrection")
    require(String(material_config["calibration_status"]) == "UNCALIBRATED",
            "prototype material values must remain explicitly uncalibrated")
    require(material_config["failure_enabled"] == false,
            "the fixed engineering diagnostic must keep material failure disabled")
    require(String(loading["control_mode"]) == "native_velocity_grip_bc",
            "coupon must use native time-dependent velocity grip conditions")

    Peridynamics.force_threads_run!()
    require(Threads.nthreads() == 1,
            "this fixed source uses one Julia thread for deterministic point-state capture")
    rest_position_m, point_volume_m3 = uniform_box(lengths_m..., h_m; center=(0.0, 0.0, 0.0))
    n_points = size(rest_position_m, 2)
    require(n_points == n_expected && length(expected_grid) == 3 &&
            n_points == prod(expected_grid),
            "uniform_box point count/grid does not match the frozen case")
    grid_shape_readback = Int[length(unique(view(rest_position_m, dim, :))) for dim in 1:3]
    require(grid_shape_readback == expected_grid,
            "uniform_box readback axis grid shape does not match the frozen case")
    require(all(isfinite, rest_position_m) && all(isfinite, point_volume_m3) &&
            all(>(0), point_volume_m3), "analytic-box points or volumes are invalid")
    require(all(isapprox.(point_volume_m3, h_m^3; rtol=1e-12, atol=0.0)),
            "uniform_box point volumes do not equal h³")
    half_lengths_m = lengths_m ./ 2
    complete_horizon_ids = Int[]
    complete_horizon_mask = falses(n_points)
    horizon_margin_m = zeros(n_points)
    for point_id in 1:n_points
        clearance_m = minimum(half_lengths_m[dim] - abs(rest_position_m[dim, point_id]) for dim in 1:3)
        horizon_margin_m[point_id] = clearance_m - horizon_m
        if clearance_m + 64 * eps(horizon_m) >= horizon_m
            push!(complete_horizon_ids, point_id)
            complete_horizon_mask[point_id] = true
        end
    end
    n_complete_horizon = length(complete_horizon_ids)
    require(n_complete_horizon == expected_interior,
            "complete-horizon interior count disagrees with the fixed analytic geometry")
    require(isapprox(Float64(geometry["expected_complete_horizon_interior_ratio"]),
                     n_complete_horizon / n_points; rtol=0.0, atol=1e-12),
            "complete-horizon interior ratio disagrees with its count")
    margin_in_h = horizon_margin_m ./ h_m
    horizon_margin_histogram_by_h = Dict{String,Int}(
        "negative" => count(x -> x < 0, margin_in_h),
        "zero_to_0_5" => count(x -> 0 <= x < 0.5, margin_in_h),
        "0_5_to_1" => count(x -> 0.5 <= x < 1.0, margin_in_h),
        "1_to_2" => count(x -> 1.0 <= x < 2.0, margin_in_h),
        "2_or_more" => count(x -> x >= 2.0, margin_in_h))

    x_min = minimum(@view rest_position_m[1, :])
    x_max = maximum(@view rest_position_m[1, :])
    grip_tol_m = h_m * 1e-8
    left_grip_ids = findall(abs.(@view(rest_position_m[1, :]) .- x_min) .<= grip_tol_m)
    right_grip_ids = findall(abs.(@view(rest_position_m[1, :]) .- x_max) .<= grip_tol_m)
    require(!isempty(left_grip_ids) && !isempty(right_grip_ids) &&
            isempty(intersect(left_grip_ids, right_grip_ids)),
            "analytic coupon does not have disjoint non-empty end grip layers")

    rho = Float64(material_config["rho_kg_m3"])
    young = Float64(material_config["E_pa"])
    nu = Float64(material_config["nu"])
    gc = Float64(material_config["Gc_j_m2"])
    require(isfinite(rho) && rho > 0 && isfinite(young) && young > 0 &&
            isfinite(nu) && 0 < nu < 0.5 && isfinite(gc) && gc > 0,
            "invalid uncalibrated prototype material parameters")
    body = Body(OSBMaterial{NoCorrection}(), rest_position_m, point_volume_m3)
    point_set!(body, :coupon_left_grip, left_grip_ids)
    point_set!(body, :coupon_right_grip, right_grip_ids)
    material!(body; horizon=horizon_m, rho=rho, E=young, nu=nu, Gc=gc)
    no_failure!(body) # material! precedes this call so failure permission cannot be overwritten.
    all(!, body.fail_permit) || error("native failure permission was not disabled for every point")
    for dim in 1:3
        velocity_ic!(body, :all_points, dim, 0.0)
        velocity_bc!(t -> 0.0, body, :coupon_left_grip, dim)
        dim == 1 || velocity_bc!(t -> 0.0, body, :coupon_right_grip, dim)
    end

    raw_native_bound_solver = VelocityVerlet(steps=1, safety_factor=0.5)
    stability_dh = Peridynamics.threads_data_handler(body, raw_native_bound_solver, 1)
    require(length(stability_dh.chunks) == 1,
            "native stable-step preflight did not produce exactly one chunk")
    native_raw_bound_s = Peridynamics.calc_stable_timestep(stability_dh, 1.0)
    require(isfinite(native_raw_bound_s) && native_raw_bound_s > 0,
            "native OSB stability bound is not finite and positive")
    requested_dt_fraction = Float64(solver_config["dt_fraction_of_native_raw_bound"])
    require(0 < requested_dt_fraction <= 0.5,
            "requested manual time step fraction must not exceed half the raw native bound")
    requested_dt_s = native_raw_bound_s * requested_dt_fraction
    ramp_duration_s = Float64(loading["ramp_duration_s"])
    hold_duration_s = Float64(loading["hold_duration_s"])
    unload_duration_s = Float64(loading["unload_duration_s"])
    post_unload_dwell_duration_s = Float64(loading["post_unload_dwell_duration_s"])
    requested_durations_s = (ramp_duration_s, hold_duration_s,
                             unload_duration_s, post_unload_dwell_duration_s)
    require(all(isfinite, requested_durations_s) && all(>=(0), requested_durations_s) &&
            ramp_duration_s > 0 && unload_duration_s > 0 &&
            post_unload_dwell_duration_s > 0,
            "load and post-unload dwell durations must be finite with positive ramp/unload/dwell times")
    positive_phase_durations_s = [duration for duration in requested_durations_s if duration > 0]
    minimum_phase_duration_s = minimum(positive_phase_durations_s)
    phase_resolution_cap_s = minimum_phase_duration_s / min_steps_per_phase
    dt_s = min(requested_dt_s, phase_resolution_cap_s)
    require(isfinite(dt_s) && dt_s > 0,
            "resolved time step is not finite and positive")
    ramp_resolution_intervals = ceil(Int, ramp_duration_s / dt_s)
    hold_resolution_intervals = ceil(Int, hold_duration_s / dt_s)
    unload_resolution_intervals = ceil(Int, unload_duration_s / dt_s)
    dwell_resolution_intervals = ceil(Int, post_unload_dwell_duration_s / dt_s)
    phase_resolution_intervals = (ramp_resolution_intervals, hold_resolution_intervals,
                                  unload_resolution_intervals, dwell_resolution_intervals)
    positive_phase_resolution_intervals = [phase_resolution_intervals[i]
                                           for i in eachindex(phase_resolution_intervals)
                                           if requested_durations_s[i] > 0]
    require(minimum(positive_phase_resolution_intervals) >= min_steps_per_phase,
            "time-step resolution cap failed to provide the declared minimum steps per phase")
    total_physical_duration_s = sum(requested_durations_s)
    n_steps = ceil(Int, total_physical_duration_s / dt_s)
    require(0 < n_steps <= max_steps,
            "fixed physical load history exceeds the declared total-step budget")
    predicted_point_state_bytes = Int64(8) * length(POINT_STATE_FIELDS) * n_points * (n_steps + 1)
    require(predicted_point_state_bytes + 5 * 1024^2 <= max_output_bytes,
            "full point-state binary plus output reserve exceeds the declared output budget")
    solver = VelocityVerlet(steps=n_steps, stepsize=dt_s, safety_factor=0.5)
    target_strain = Float64(loading["target_engineering_strain"])
    require(isfinite(target_strain) && target_strain > 0 && target_strain < 1e-2,
            "fixed reversible engineering strain must be small, positive, and below 1%")
    grip_span_m = x_max - x_min
    target_grip_displacement_m = target_strain * grip_span_m
    velocity_bc!(t -> prescribed_step_velocity_m_s(t, dt_s, ramp_duration_s,
                                                    hold_duration_s, unload_duration_s,
                                                    target_grip_displacement_m),
                 body, :coupon_right_grip, :x)

    recorder = CouponRecorder(n_points, n_steps)
    options = CouponOptions(false, recorder, rest_position_m, point_volume_m3, rho, solver)
    job = Job(body, solver, options)
    dh = submit(job; quiet=true)
    require(length(dh.chunks) == 1, "native submit returned more than one coupon chunk")
    require(all(recorder.seen) && recorder.last_step == n_steps,
            "native hook did not capture every state from 0 through n_steps")
    require(recorder.directed_bond_count > 0 &&
            recorder.active_directed_bonds_at_reference == recorder.directed_bond_count,
            "initial native graph is empty or contains inactive directed bonds")
    require(all(isfinite, recorder.displacement_m) && all(isfinite, recorder.velocity_m_s) &&
            all(isfinite, recorder.velocity_half_m_s) && all(isfinite, recorder.acceleration_m_s2) &&
            all(isfinite, recorder.b_int_N_m3) && all(isfinite, recorder.b_ext_N_m3) &&
            all(isfinite, recorder.strain_energy_density_J_m3) && all(isfinite, recorder.damage),
            "captured native state contains non-finite values")
    require(all(iszero, recorder.damage) &&
            all(recorder.active_bonds_per_point .> 0),
            "failure-disabled coupon developed nonzero damage or a zero-active-bond point")
    max_right_grip_path_error_m = 0.0
    max_right_grip_point_spread_m = 0.0
    max_fixed_grip_displacement_m = 0.0
    for state_index in axes(recorder.displacement_m, 3)
        step = state_index - 1
        time_s = step * dt_s
        target_u_m = target_displacement_m(time_s, ramp_duration_s,
                                           hold_duration_s, unload_duration_s,
                                           target_grip_displacement_m)
        right_u_x = recorder.displacement_m[right_grip_ids, 1, state_index]
        measured_right_u_x = arithmetic_mean(right_u_x)
        max_right_grip_path_error_m = max(max_right_grip_path_error_m,
                                          abs(measured_right_u_x - target_u_m))
        max_right_grip_point_spread_m = max(max_right_grip_point_spread_m,
                                             maximum(abs, right_u_x .- measured_right_u_x))
        max_fixed_grip_displacement_m = max(
            max_fixed_grip_displacement_m,
            maximum(abs, recorder.displacement_m[left_grip_ids, :, state_index]),
            maximum(abs, recorder.displacement_m[right_grip_ids, 2:3, state_index]))
    end
    grip_path_tolerance_m = max(256 * eps(target_grip_displacement_m), 1e-14)
    require(max_right_grip_path_error_m <= grip_path_tolerance_m &&
            max_right_grip_point_spread_m <= grip_path_tolerance_m,
            "native right grip positions do not reproduce the continuous physical-time target path")
    require(max_fixed_grip_displacement_m <= grip_path_tolerance_m,
            "a zero-velocity grip component moved during native drift")
    actual_final_right_u_m = arithmetic_mean(recorder.displacement_m[right_grip_ids, 1, end])
    require(abs(actual_final_right_u_m) <= grip_path_tolerance_m,
            "native velocity grip did not return to the prescribed zero-displacement state")

    mkpath(output_dir)
    point_state_binary = write_point_state_binary(
        joinpath(output_dir, "point_state_history.f64le"),
        joinpath(output_dir, "point_state_layout.toml"), recorder, dt_s,
        max_output_bytes)
    require_output_budget(output_dir, max_output_bytes)
    write_point_reference_csv(joinpath(output_dir, "point_reference.csv"), recorder,
                              rest_position_m, point_volume_m3, rho,
                              complete_horizon_mask, horizon_margin_m,
                              left_grip_ids, right_grip_ids)
    require_output_budget(output_dir, max_output_bytes)
    ledger = write_energy_history(joinpath(output_dir, "energy_reaction_history.csv"),
                                  recorder, rest_position_m, point_volume_m3,
                                  left_grip_ids, right_grip_ids, rho, dt_s)
    require_output_budget(output_dir, max_output_bytes)
    summary = Dict{String,Any}(
        "schema" => "pgfrag.free_evolution_coupon.result.v1",
        "status" => "ENGINEERING_DIAGNOSTIC_COMPLETED",
        "experiment_kind" => "single_body_native_OSB_reversible_coupon_with_free_interior_dwell",
        "run_status_before_launch" => String(case["run_status"]),
        "geometry" => "uniform_box analytic cell-centre grid; no CAD or imported artifact",
        "shape_m" => lengths_m,
        "grip_span_m" => grip_span_m,
        "grid_shape_readback" => grid_shape_readback,
        "point_spacing_h_m" => h_m,
        "horizon_delta_m" => horizon_m,
        "horizon_to_spacing_ratio_delta_over_h" => horizon_m / h_m,
        "point_count" => n_points,
        "point_state_binary_bytes" => point_state_binary.byte_count,
        "point_state_binary_sha256" => point_state_binary.sha256,
        "point_state_layout_toml" => point_state_binary.layout_path,
        "point_reference_csv" => "point_reference.csv",
        "point_state_encoding" => "self-describing little-endian Float64 tensor; see point_state_layout.toml",
        "output_byte_budget" => max_output_bytes,
        "point_state_bytes_predicted_before_solver" => predicted_point_state_bytes,
        "output_bytes_measured_before_summary" => output_directory_bytes(output_dir),
        "point_volume_m3" => h_m^3,
        "complete_spherical_horizon_interior_point_count" => n_complete_horizon,
        "complete_spherical_horizon_interior_point_ratio" => n_complete_horizon / n_points,
        "complete_horizon_margin_histogram_in_units_of_h" => horizon_margin_histogram_by_h,
        "left_grip_point_count" => length(left_grip_ids),
        "right_grip_point_count" => length(right_grip_ids),
        "material_model" => "OSBMaterial{NoCorrection} with native CriticalStretch damage model",
        "surface_correction" => "NoCorrection",
        "surface_correction_gap" => "surface support truncation is not corrected; no surface-response qualification is implied",
        "E_pa" => young,
        "nu" => nu,
        "rho_kg_m3" => rho,
        "Gc_j_m2" => gc,
        "calibration_status" => "UNCALIBRATED prototype values; no material calibration was performed",
        "failure_enabled" => false,
        "precrack_or_manual_bond_edit" => false,
        "all_point_affine_prescription" => false,
        "target_engineering_strain" => target_strain,
        "target_grip_displacement_m" => target_grip_displacement_m,
        "sampled_peak_right_grip_displacement_m" => ledger.max_abs_right_grip_u_x_m,
        "max_right_grip_path_error_m" => max_right_grip_path_error_m,
        "max_right_grip_point_spread_m" => max_right_grip_point_spread_m,
        "max_zero_velocity_grip_displacement_m" => max_fixed_grip_displacement_m,
        "interior_kinematics" => "free under the native VelocityVerlet force and acceleration updates",
        "grip_condition" => "native velocity_bc! on x-end grip sets; left xyz fixed, right yz fixed and x velocity is the secant of a continuous fixed-physical-time displacement ramp/hold/unload path; both grip sets remain at zero velocity during the post-unload dwell",
        "post_unload_dwell_boundary_status" => "fixed-grip dwell; interior degrees of freedom remain free, specimen is not free-boundary",
        "ramp_duration_s" => ramp_duration_s,
        "hold_duration_s" => hold_duration_s,
        "unload_duration_s" => unload_duration_s,
        "post_unload_dwell_duration_s" => post_unload_dwell_duration_s,
        "total_physical_duration_s" => total_physical_duration_s,
        "solver_end_time_s" => solver.n_steps * solver.Δt,
        "solver_end_time_overrun_s" => solver.n_steps * solver.Δt - total_physical_duration_s,
        "loading_path_sampling_note" => "ramp, hold, unload, and dwell durations do not change with dt; the measured peak on discrete t_n samples may be below the continuous target peak if no solver time lands on the hold interval",
        "ramp_resolution_intervals_ceiling" => ramp_resolution_intervals,
        "hold_resolution_intervals_ceiling" => hold_resolution_intervals,
        "unload_resolution_intervals_ceiling" => unload_resolution_intervals,
        "post_unload_dwell_resolution_intervals_ceiling" => dwell_resolution_intervals,
        "step_count" => solver.n_steps,
        "dt_s" => solver.Δt,
        "native_raw_stability_bound_s" => native_raw_bound_s,
        "requested_dt_fraction_of_native_raw_bound" => requested_dt_fraction,
        "accepted_dt_fraction_of_native_raw_bound" => solver.Δt / native_raw_bound_s,
        "requested_dt_s_before_phase_resolution_cap" => requested_dt_s,
        "phase_resolution_cap_s" => phase_resolution_cap_s,
        "directed_native_bond_count" => recorder.directed_bond_count,
        "active_directed_native_bond_count_at_reference" => recorder.active_directed_bonds_at_reference,
        "inactive_native_bond_count_at_end" => count(!, only(dh.chunks).storage.bond_active),
        "contact_status" => "NOT_RUN; one body and no contact registrations",
        "grip_reaction_method" => "native discrete velocity-constraint impulse on the specimen: m*(v_half_after_BC - (v_full_previous + a_previous*dt/2)); impulse/dt is the interval-average constrained force, not a continuum endpoint reaction; support reaction has the opposite sign",
        "state_time_convention" => "U, full V, acceleration, force density, and OSB strain energy are at t_n; for n>=1 native velocity_half is centered at (n-1/2)*dt after BC evaluation at t_n; at n=0 velocity_half is the initialized zero placeholder",
        "discrete_boundary_work_method" => "sum over constrained grip DOFs of 0.5*m*(v_half_after_BC^2 - v_half_before_BC^2) at the native half-step velocity overwrite",
        "actuator_work_trapezoid_status" => "SURROGATE_DIAGNOSTIC_ONLY; trapezoidal average-impulse-force times actual grip displacement; excluded from the energy residual",
        "free_dof_impulse_check" => "each unconstrained native DOF checks m*(v_half_current - (v_full_previous + a_previous*dt/2)) against floating-point tolerance",
        "momentum_time_level" => "full-step momentum at t_n; each row's BC impulse spans [t_(n-1),t_n] and is checked with the trapezoidal endpoint force impulse from native b_int/b_ext",
        "strain_energy_method" => "native Peridynamics.export_field(Val(:strain_energy_density), OSBMaterial, BondSystem, storage, paramsetup, t), integrated as sum(Psi_i * point_volume_i)",
        "energy_balance_status" => "DIAGNOSTIC_ONLY; native OSB strain energy and full-step kinetic energy are compared against accumulated half-step BC kinetic-energy jumps; VelocityVerlet state times are staggered, so this is not a qualified energy-closure test",
        "fracture_dissipation_inference" => "NOT_PERMITTED; failure was disabled and this reversible elastic diagnostic has no fracture-energy ledger",
        "energy_ledger" => Dict("total_mass_kg" => ledger.total_mass_kg,
                                "initial_native_strain_plus_kinetic_energy_J" => ledger.initial_total_energy_J,
                                "final_native_strain_energy_J" => ledger.final_osb_strain_energy_J,
                                "final_kinetic_energy_J" => ledger.final_kinetic_energy_J,
                                "final_cumulative_discrete_bc_kinetic_work_J" => ledger.final_cumulative_discrete_bc_kinetic_work_J,
                                "final_actuator_work_trapezoid_surrogate_J" => ledger.final_actuator_work_trapezoid_surrogate_J,
                                "final_energy_balance_residual_diagnostic_J" => ledger.final_energy_balance_residual_diagnostic_J,
                                "max_free_dof_impulse_residual_N_s" => ledger.max_free_dof_impulse_residual_N_s,
                                "max_momentum_update_residual_N_s" => ledger.max_momentum_update_residual_N_s,
                                "initial_momentum_kg_m_s" => ledger.initial_momentum_kg_m_s,
                                "final_momentum_kg_m_s" => ledger.final_momentum_kg_m_s,
                                "final_left_grip_u_x_mean_m" => ledger.final_left_grip_u_x_mean_m,
                                "final_right_grip_u_x_mean_m" => ledger.final_right_grip_u_x_mean_m,
                                "maximum_absolute_right_grip_u_x_m" => ledger.max_abs_right_grip_u_x_m),
        "real_object_physics_qualification" => "UNKNOWN; no measured coupon or artifact data were supplied or fitted",
        "source_freeze_path" => relpath(SOURCE_FREEZE_PATH, REPO_ROOT),
        "source_freeze_sha256" => freeze_record.freeze_sha256,
        "verified_source_sha256" => freeze_record.actual_hashes,
        "point_state_binary" => "point_state_history.f64le",
        "energy_reaction_history_csv" => "energy_reaction_history.csv")
    write_summary(joinpath(output_dir, "summary.toml"), summary)
    require_output_budget(output_dir, max_output_bytes)
    println("FREE_EVOLUTION_COUPON_COMPLETED points=$n_points interior=$n_complete_horizon/$n_points steps=$(solver.n_steps) dt_s=$(solver.Δt) bonds=$(recorder.directed_bond_count)")
    return nothing
end

try
    run_coupon(ARGS)
catch err
    showerror(stderr, err, catch_backtrace())
    println(stderr)
    exit(1)
end
