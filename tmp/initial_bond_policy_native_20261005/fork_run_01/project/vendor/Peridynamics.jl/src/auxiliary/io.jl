@inline default_export_fields() = [:displacement, :damage]

@inline function default_fields_spec(ms::AbstractMultibodySetup)
    fields_spec = Dict{Symbol,Vector{Symbol}}()
    for body_name in each_body_name(ms)
        fields_spec[body_name] = default_export_fields()
    end
    return fields_spec
end

function get_export_fields(body::AbstractBody, solver::AbstractTimeSolver,
                           o::Dict{Symbol,Any})
    if haskey(o, :fields)
        fields = extract_export_fields(o[:fields])
    else
        fields = default_export_fields()
    end

    check_export_fields(body, solver, fields)

    return fields
end

function get_export_fields(ms::AbstractMultibodySetup, solver::AbstractTimeSolver,
                           o::Dict{Symbol,Any})
    if haskey(o, :fields)
        fields_spec = extract_fields_spec(ms, o[:fields])
    else
        fields_spec = default_fields_spec(ms)
    end

    check_export_fields(ms, solver, fields_spec)

    return fields_spec
end

@inline function extract_export_fields(o)
    fields::Vector{Symbol} = _extract_export_fields(o)
    return fields
end

@inline _extract_export_fields(fields::Symbol) = [fields]
@inline _extract_export_fields(fields::NTuple{N,Symbol}) where {N} = [f for f in fields]
@inline _extract_export_fields(fields::Vector{Symbol}) = fields
@inline function _extract_export_fields(arbitrary_fields::T) where {T}
    fields::Vector{Symbol} = try
        convert(Vector{Symbol}, arbitrary_fields)
    catch
        msg = "wrong type of specification for keyword `fields`!\n"
        msg *= "Cannot convert input of type `$(T)` to `Vector{Symbol}`!\n"
        throw(ArgumentError(msg))
    end
    return fields
end

@inline function extract_fields_spec(ms::AbstractMultibodySetup, o)
    fields_spec::Dict{Symbol,Vector{Symbol}} = _extract_fields_spec(ms, o)
    for body_name in each_body_name(ms)
        if !haskey(fields_spec, body_name)
            fields_spec[body_name] = default_export_fields()
        end
    end
    return fields_spec
end

function _extract_fields_spec(ms::AbstractMultibodySetup, o::Dict{Symbol,T}) where {T}
    fields_spec = Dict{Symbol,Vector{Symbol}}()
    for (body_name, fields) in o
        check_if_bodyname_is_defined(ms, body_name)
        fields_spec[body_name] = extract_export_fields(fields)
    end
    return fields_spec
end

function _extract_fields_spec(ms::AbstractMultibodySetup, o)
    fields = extract_export_fields(o)
    fields_spec = Dict{Symbol,Vector{Symbol}}()
    for body_name in each_body_name(ms)
        fields_spec[body_name] = fields
    end
    return fields_spec
end

function check_export_fields(body::AbstractBody, solver::AbstractTimeSolver,
                             fields::Vector{Symbol})
    S = storage_type(body)
    check_export_fields(S, fields)
    return nothing
end

function check_export_fields(ms::AbstractMultibodySetup, solver::AbstractTimeSolver,
                             fields_spec::Dict{Symbol,Vector{Symbol}})
    for body_name in each_body_name(ms)
        body = get_body(ms, body_name)
        S = storage_type(body)
        fields = fields_spec[body_name]
        check_export_fields(S, fields)
    end
    return nothing
end

function check_export_fields(::Type{S}, fields::Vector{Symbol}) where {S}
    allowed_fields = point_data_fields(S)
    for f in fields
        if !in(f, allowed_fields) && !custom_field(S, f)
            msg = "unknown point data field `:$(f)` specified for export!\n"
            msg *= "If you intend to export a custom field, please define a method:"
            msg *= "\n   `custom_field(::Type{MyStorage}, ::Val{:$(f)}) = true`\n"
            msg *= "Otherwise, see here all available point data fields of $S:\n"
            for allowed_name in allowed_fields
                msg *= "  - $allowed_name\n"
            end
            throw(ArgumentError(msg))
        end
    end
    return nothing
end

custom_field(S::Type{<:AbstractStorage}, field::Symbol) = custom_field(S, Val(field))
# this function can be specialized to indicate custom fields
# if this is not done, the export_field function will error out!
custom_field(::Type{<:AbstractStorage}, ::Val{field}) where {field} = false

function _export_results(options::AbstractJobOptions, chunk::AbstractBodyChunk,
                         chunk_id::Int, n_chunks::Int, n::Int, t::Float64)
    throw(ArgumentError("VTK export is disabled in the PGFrag headless backend; use the PGFrag CSV/NPZ result adapter"))
end

function msg_export_fields(fields::Vector{Symbol}; indentation=2)
    return msg_vec("exported fields", fields; continuation_label=" "^15, indentation)
end

function msg_export_fields(fields_spec::Dict{Symbol,Vector{Symbol}}; indentation=2)
    msg = " "^indentation * "EXPORTED FIELDS PER BODY\n"
    for (body_name, fields) in fields_spec
        body_name_str = string(body_name)
        n = length(body_name_str)
        ind = indentation + 2
        msg *= msg_vec(body_name_str, fields; continuation_label=" "^n, indentation=ind)
    end
    return msg
end
