using InvertedIndices
using JuMP
include("stream.jl")

#=
  _separator, separator       – phase split constraints
  same_stream                 – composition equality between streams
  pressure_reduction          – throttle valve (Joule-Thomson)
  create_vap, create_liq      – fix stream phase
  same_pressure               – pressure equality between two streams
  total_streams               – allocate N stream variables on model
  total_same_pressure_idx     – pressure equality across an index list
  all_vapor_streams, all_liquid_streams
  simple_thermal_exergy       – exergy destruction for one HX
  to_print                    – post-solve stream summary
=#


function _separator(
    m::Model,
    stream_2P_in::create_stream,
    stream_vap_out::create_stream,
    stream_liq_out::create_stream
)
    n = length(stream_2P_in.x)
    @constraint(m, stream_2P_in.T == stream_vap_out.T)
    @constraint(m, stream_2P_in.T == stream_liq_out.T)

    JuMP.fix(stream_vap_out.q, 0.0; force = true)
    JuMP.fix(stream_liq_out.q, 1.0; force = true)

    for i in 1:(n - 1)
        @constraint(m, stream_liq_out.z[i] == stream_2P_in.x[i])
        @constraint(m, stream_vap_out.z[i] == stream_2P_in.y[i])
    end

    return (stream_2P_in, stream_vap_out, stream_liq_out)
end

struct separator
    stream_2P_in::create_stream
    stream_vap_out::create_stream
    stream_liq_out::create_stream

    function separator(
        m::Model,
        stream_2P_in::create_stream,
        stream_vap_out::create_stream,
        stream_liq_out::create_stream
    )
        data = _separator(m, stream_2P_in, stream_vap_out, stream_liq_out)
        return new(data...)
    end
end

function same_stream(
    m::Model,
    stream_1::create_stream,
    stream_2::create_stream
)
    n = length(stream_1.z)
    for i in 1:(n - 1)
        @constraint(m, stream_1.z[i] == stream_2.z[i])
    end
end

function all_same_stream(
    m::Model,
    streams::Vector{create_stream},
    idx::Vector{Int}
)
    any(x -> x > length(streams), idx) && error("stream index out of range: max valid is $(length(streams)), got $(filter(x -> x > length(streams), idx))")
    for i in 1:(length(idx) - 1)
        same_stream(m, streams[idx[i]], streams[idx[i + 1]])
    end
end

function _press_reduc(
    m::Model,
    stream::create_stream
)
    pRef = m[:prop].pRef
    P_low = @variable(m, lower_bound = 0.1, upper_bound = pRef)
    P_cons = @constraint(m, stream.P == P_low)
    return (stream, P_low, P_cons)
end

struct pressure_reduction
    stream::create_stream
    P_low::VariableRef
    P_cons::ConstraintRef

    function pressure_reduction(m::Model, stream::create_stream)
        data = _press_reduc(m, stream)
        return new(data...)
    end
end


function create_vap(
    m::Model,
    prop::Properties,
    stream::create_stream
)
    n = prop.n
    α = prop.α

    @constraint(m, stream.q == 0.0)
    @constraint(m, stream.q_LH == 0)
    @constraint(m, stream.β * sum(stream.z[i] / α[i] for i in 1:n) == 1)

    return stream
end

function create_liq(
    m::Model,
    prop::Properties,
    stream::create_stream
)
    n = prop.n
    α = prop.α

    @constraint(m, stream.q == 1)
    @constraint(m, stream.q_LH == 1)
    @constraint(m, stream.β == sum(α[i] * stream.z[i] for i in 1:n))

    return stream 
end


function same_pressure(
    m::Model,
    stream1::create_stream,
    stream2::create_stream
)
    @constraint(m, stream1.P == stream2.P)
end

function to_print(stream_name)
    n = length(stream_name.z)
    println("Temperature: ", round(value(stream_name.T), digits = 4))
    println("z: ", [round(value(stream_name.z[i]); digits = 4) for i in 1:n])
    println("Pressure: ", round(value(stream_name.P), digits = 4))
    println("Liquid Fraction: ", round(value(stream_name.q), digits = 4))
end

function total_streams(
    m::Model,
    nop::Int,
)
    streams = Vector{create_stream}(undef, nop)
    for i in 1:nop
        streams[i] = create_stream(m)
    end

    return streams
end

function total_same_pressure_idx(
    m::Model,
    streams::Vector{create_stream},
    idx::Vector{Int}
)
    any(x -> x > length(streams), idx) && error("same-pressure index out of range: max valid is $(length(streams)), got $(filter(x -> x > length(streams), idx))")
    for i in 1:(length(idx) - 1)
        @constraint(m, streams[idx[i]].P == streams[idx[i + 1]].P)
    end
end

function all_vapor_streams(
    m::Model,
    streams::Vector{create_stream},
    idx::Vector{Int}
)
    prop = m[:prop]
    any(x -> x > length(streams), idx) && error("vapor stream index out of range: max valid is $(length(streams)), got $(filter(x -> x > length(streams), idx))")
    for i in 1:length(idx)
        streams[idx[i]] = create_vap(m, prop, streams[idx[i]])
    end
end

function all_liquid_streams(
    m::Model,
    streams::Vector{create_stream},
    idx::Vector{Int}
)
    prop = m[:prop]
    any(x -> x > length(streams), idx) && error("liquid stream index out of range: max valid is $(length(streams)), got $(filter(x -> x > length(streams), idx))")
    for i in 1:length(idx)
        streams[idx[i]] = create_liq(m, prop, streams[idx[i]])
    end
end

# function detailed_thermal_exergy(
#     model::Model,
#     T_hot_high::Union{VariableRef,AffExpr,NonlinearExpr,Float64},  # T1_in  > T1_out
#     T_hot_low::Union{VariableRef,AffExpr,NonlinearExpr,Float64},   # T1_out
#     T_cold_low::Union{VariableRef,AffExpr,NonlinearExpr,Float64},  # T2_in  < T2_out
#     T_cold_high::Union{VariableRef,AffExpr,NonlinearExpr,Float64}, # T2_out
#     Q::Union{VariableRef,AffExpr,NonlinearExpr,Float64}
# )
#     T0 = 25.0 + 273.15
#     ϵ      = @variable(model, lower_bound = 0.0,    upper_bound = 5.0)
#     # s_hot  = ln(T1_in/T1_out)  / (T1_in - T1_out)  — log-mean inverse temperature, hot side
#     # s_cold = ln(T2_out/T2_in) / (T2_out - T2_in)  — log-mean inverse temperature, cold side
#     s_hot  = @variable(model, lower_bound = 1/600,  upper_bound = 1/273.15)
#     s_cold = @variable(model, lower_bound = 1/600,  upper_bound = 1/273.15)
#     @constraint(model, s_hot  * (T_hot_high  - T_hot_low)  == log(T_hot_high)  - log(T_hot_low))
#     @constraint(model, s_cold * (T_cold_high - T_cold_low) == log(T_cold_high) - log(T_cold_low))
#     @constraint(model, ϵ == T0 * Q * (s_cold - s_hot))
#     return ϵ
# end

function simple_thermal_exergy(
    model::Model,
    T_hot_high::Union{VariableRef,AffExpr, NonlinearExpr, Float64}, # T_hot_high > T_hot_low
    T_hot_low::Union{VariableRef,AffExpr, NonlinearExpr, Float64}, 
    T_cold_low::Union{VariableRef,AffExpr, NonlinearExpr, Float64}, # T_cold_low > T_cold_high
    T_cold_high::Union{VariableRef,AffExpr, NonlinearExpr, Float64},
    Q::Union{VariableRef,AffExpr, NonlinearExpr, Float64}
)
    ϵ = @variable(model, lower_bound = 0.0, upper_bound = 5.0)
    b1 = @variable(model, lower_bound = 1/600, upper_bound = 1/273.15)
    b2 = @variable(model, lower_bound = 1/600, upper_bound = 1/273.15)
    b3 = @variable(model, lower_bound = 1/600, upper_bound = 1/273.15)
    b4 = @variable(model, lower_bound = 1/600, upper_bound = 1/273.15)
    @constraints(
        model, begin
        b1 * T_cold_high - 1 == 0
        b2 * T_hot_high - 1 == 0
        b3 * T_cold_low - 1 == 0
        b4 * T_hot_low - 1 == 0
        ϵ == 0.5 * Q * (25 + 273.15) * ((b1 - b2) + (b3 - b4))
    end)
    return ϵ
end
