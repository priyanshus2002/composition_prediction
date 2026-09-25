using JuMP
# Requires utils.jl (which includes stream.jl) to be loaded first.

abstract type HeatPumpConfig end

# ─── 1 Source / 1 Sink ────────────────────────────────────────────────────────

struct InternalCascade_1so1si <: HeatPumpConfig
    duty::Float64
    source_temp::Vector{Float64}     # [high, low]
    sink_temp::Vector{Float64}       # [low, high]
    min_TP::Float64

    function InternalCascade_1so1si(;
        duty::Real = 1.0,
        source_temp::Vector{Float64} = [300.0, 300.0],
        sink_temp::Vector{Float64} = [400.0, 400.0],
        min_TP::Real = 1.0
    )
        length(source_temp) == 2 || error("source_temp must have 2 elements: [high, low]")
        length(sink_temp)   == 2 || error("sink_temp must have 2 elements: [low, high]")
        new(Float64(duty), source_temp, sink_temp, Float64(min_TP))
    end
end

n_streams(::InternalCascade_1so1si) = 9

function build!(m::Model, config::InternalCascade_1so1si)
    prop = m[:prop]
    streams = total_streams(m, n_streams(config))
    (; duty, source_temp, sink_temp, min_TP) = config

    total_same_pressure_idx(m, streams, [1, 2, 3, 4, 5])
    total_same_pressure_idx(m, streams, [6, 7, 8, 9])

    ϵ1 = @variable(m, lower_bound = -duty, upper_bound = duty)
    ϵ2 = @variable(m, lower_bound = -duty, upper_bound = duty)
    ϵ3 = @variable(m, lower_bound = -duty, upper_bound = duty)

    streams[1] = create_vap(m, prop, streams[1])
    same_stream(m, streams[1], streams[2])
    @constraint(m, streams[2].q_LH == 0.5)
    separator(m, streams[2], streams[3], streams[4])

    streams[5] = create_liq(m, prop, streams[5])
    same_stream(m, streams[3], streams[5])

    same_stream(m, streams[5], streams[6])
    throttle = pressure_reduction(m, streams[6])
    streams[6] = create_liq(m, prop, streams[6])
    streams[7] = create_vap(m, prop, streams[7])
    same_stream(m, streams[6], streams[7])

    same_stream(m, streams[2], streams[8])
    same_stream(m, streams[8], streams[9])
    streams[9] = create_vap(m, prop, streams[9])
    @constraint(m, streams[8].q_LH == 0.5)

    @constraints(m, begin
        streams[1].T >= sink_temp[2] + 5
        streams[2].T >= sink_temp[1] + 5
        streams[6].T <= source_temp[2] - 5
        streams[7].T <= source_temp[1] - 5
        streams[5].T >= streams[8].T + 5
    end)

    @constraints(m, begin
        ϵ1 == simple_thermal_exergy(m, streams[1].T, streams[2].T, sink_temp[1], sink_temp[2], duty)
        ϵ2 == simple_thermal_exergy(m, source_temp[1], source_temp[2], streams[6].T, streams[7].T, duty)
        ϵ3 == simple_thermal_exergy(m, streams[3].T, streams[5].T, streams[8].T, streams[9].T, duty)
    end)

    @constraint(m, throttle.P_low >= min_TP)

    m[:streams] = streams
    return ϵ1 + ϵ2 + ϵ3
end

# ─── 1 Source / 2 Sinks ───────────────────────────────────────────────────────

struct InternalCascade_1so2si <: HeatPumpConfig
    duty::Vector{Float64}                  # [duty_sink1, duty_sink2]
    source_temp::Vector{Float64}           # [high, low]
    sink_temp::Vector{Vector{Float64}}     # [[low1, high1], [low2, high2]]
    min_TP::Float64

    function InternalCascade_1so2si(;
        duty::Vector{Float64} = [1.0, 1.0],
        source_temp::Vector{Float64} = [300.0, 300.0],
        sink_temp::Vector{Vector{Float64}} = [[440.0, 450.0], [420.0, 440.0]],
        min_TP::Real = 1.0
    )
        length(duty)       == 2 || error("duty must have 2 elements: [duty_sink1, duty_sink2]")
        length(source_temp) == 2 || error("source_temp must have 2 elements: [high, low]")
        length(sink_temp)  == 2 || error("sink_temp must have 2 entries")
        all(length(st) == 2 for st in sink_temp) || error("each sink_temp entry must be [low, high]")
        new(duty, source_temp, sink_temp, Float64(min_TP))
    end
end

n_streams(::InternalCascade_1so2si) = 10

function build!(m::Model, config::InternalCascade_1so2si)
    prop = m[:prop]
    streams = total_streams(m, n_streams(config))
    (; duty, source_temp, sink_temp, min_TP) = config
    d1, d2 = duty[1], duty[2]
    total_d = d1 + d2

    total_same_pressure_idx(m, streams, [1, 2, 3, 4, 5, 6])
    total_same_pressure_idx(m, streams, [7, 8, 9, 10])

    ϵ1 = @variable(m, lower_bound = -total_d, upper_bound = total_d)
    ϵ2 = @variable(m, lower_bound = -total_d, upper_bound = total_d)
    ϵ3 = @variable(m, lower_bound = -total_d, upper_bound = total_d)
    ϵ4 = @variable(m, lower_bound = -total_d, upper_bound = total_d)

    streams[1] = create_vap(m, prop, streams[1])
    same_stream(m, streams[1], streams[2])
    @constraint(m, streams[2].q_LH == d1 / (2 * d1 + d2))
    separator(m, streams[2], streams[3], streams[4])

    @constraint(m, streams[5].q_LH == d2 / total_d)
    same_stream(m, streams[3], streams[5])

    streams[6] = create_liq(m, prop, streams[6])
    same_stream(m, streams[5], streams[6])
    throttle = pressure_reduction(m, streams[7])
    streams[7] = create_liq(m, prop, streams[7])
    streams[8] = create_vap(m, prop, streams[8])
    same_stream(m, streams[6], streams[7])
    same_stream(m, streams[7], streams[8])

    same_stream(m, streams[2], streams[9])
    same_stream(m, streams[9], streams[10])
    streams[10] = create_vap(m, prop, streams[10])
    @constraint(m, streams[9].q_LH == d1 / total_d)

    @constraints(m, begin
        streams[1].T >= sink_temp[1][2] + 5
        streams[2].T >= sink_temp[1][1] + 5
        streams[3].T >= sink_temp[2][2] + 5
        streams[5].T >= sink_temp[2][1] + 5
        streams[5].T >= streams[10].T + 5
        streams[6].T >= streams[9].T + 5
        streams[7].T <= source_temp[2] - 5
        streams[8].T <= source_temp[1] - 5
    end)

    @constraints(m, begin
        ϵ1 == simple_thermal_exergy(m, streams[1].T, streams[2].T, sink_temp[1][1], sink_temp[1][2], d1)
        ϵ2 == simple_thermal_exergy(m, streams[3].T, streams[5].T, sink_temp[2][1], sink_temp[2][2], d2)
        ϵ3 == simple_thermal_exergy(m, source_temp[1], source_temp[2], streams[7].T, streams[8].T, total_d)
        ϵ4 == simple_thermal_exergy(m, streams[5].T, streams[6].T, streams[9].T, streams[10].T, d1)
    end)

    @constraint(m, throttle.P_low >= min_TP)

    m[:streams] = streams
    return ϵ1 + ϵ2 + ϵ3 + ϵ4
end


# ─── 1 Source / 2 Sinks Alternate ───────────────────────────────────────────────────────

struct InternalCascade_1so2si_new <: HeatPumpConfig
    duty::Vector{Float64}                  # [duty_sink1, duty_sink2]
    source_temp::Vector{Float64}           # [high, low]
    sink_temp::Vector{Vector{Float64}}     # [[low1, high1], [low2, high2]]
    min_TP::Float64

    function InternalCascade_1so2si_new(;
        duty::Vector{Float64} = [1.0, 1.0],
        source_temp::Vector{Float64} = [300.0, 300.0],
        sink_temp::Vector{Vector{Float64}} = [[440.0, 450.0], [420.0, 430.0]],
        min_TP::Real = 1.0
    )
        length(duty)       == 2 || error("duty must have 2 elements: [duty_sink1, duty_sink2]")
        length(source_temp) == 2 || error("source_temp must have 2 elements: [high, low]")
        length(sink_temp)  == 2 || error("sink_temp must have 2 entries")
        all(length(st) == 2 for st in sink_temp) || error("each sink_temp entry must be [low, high]")
        new(duty, source_temp, sink_temp, Float64(min_TP))
    end
end

n_streams(::InternalCascade_1so2si_new) = 10

function build!(m::Model, config::InternalCascade_1so2si_new)
    prop = m[:prop]
    streams = total_streams(m, n_streams(config))
    (; duty, source_temp, sink_temp, min_TP) = config
    d1, d2 = duty[1], duty[2]
    total_d = d1 + d2

    total_same_pressure_idx(m, streams, [1, 2, 3, 4, 5, 6])
    total_same_pressure_idx(m, streams, [7, 8, 9, 10])

    ϵ1 = @variable(m, lower_bound = -total_d, upper_bound = total_d)
    ϵ2 = @variable(m, lower_bound = -total_d, upper_bound = total_d)
    ϵ3 = @variable(m, lower_bound = -total_d, upper_bound = total_d)
    ϵ4 = @variable(m, lower_bound = -total_d, upper_bound = total_d)

    streams[1] = create_vap(m, prop, streams[1])
    same_stream(m, streams[1], streams[2])
    @constraint(m, streams[2].q_LH == d1 / (2 * d1 + d2))
    separator(m, streams[2], streams[3], streams[4])

    @constraint(m, streams[5].q_LH == d1 / total_d)
    same_stream(m, streams[3], streams[5])

    streams[6] = create_liq(m, prop, streams[6])
    same_stream(m, streams[5], streams[6])
    throttle = pressure_reduction(m, streams[7])
    streams[7] = create_liq(m, prop, streams[7])
    streams[8] = create_vap(m, prop, streams[8])
    same_stream(m, streams[6], streams[7])
    same_stream(m, streams[7], streams[8])

    same_stream(m, streams[2], streams[9])
    same_stream(m, streams[9], streams[10])
    streams[10] = create_vap(m, prop, streams[10])
    @constraint(m, streams[9].q_LH == d1 / total_d)

    @constraints(m, begin
        streams[1].T >= sink_temp[1][2] + 5
        streams[2].T >= sink_temp[1][1] + 5
        streams[5].T >= sink_temp[2][2] + 5
        streams[6].T >= sink_temp[2][1] + 5
        streams[3].T >= streams[10].T + 5
        streams[5].T >= streams[9].T + 5
        streams[7].T <= source_temp[2] - 5
        streams[8].T <= source_temp[1] - 5
    end)

    @constraints(m, begin
        ϵ1 == simple_thermal_exergy(m, streams[1].T, streams[2].T, sink_temp[1][1], sink_temp[1][2], d1)
        ϵ2 == simple_thermal_exergy(m, streams[5].T, streams[6].T, sink_temp[2][1], sink_temp[2][2], d2)
        ϵ3 == simple_thermal_exergy(m, source_temp[1], source_temp[2], streams[7].T, streams[8].T, total_d)
        ϵ4 == simple_thermal_exergy(m, streams[3].T, streams[5].T, streams[9].T, streams[10].T, d1)
    end)

    @constraint(m, throttle.P_low >= min_TP)

    m[:streams] = streams
    return ϵ1 + ϵ2 + ϵ3 + ϵ4
end


# ─── 2 Sources / 2 Sinks ──────────────────────────────────────────────────────

struct InternalCascade_2so2si <: HeatPumpConfig
    duty::Vector{Float64}                  # [duty_sink, duty_source]
    source_temp::Vector{Vector{Float64}}   # [[high, low] source2, [high, low] source1]
    sink_temp::Vector{Vector{Float64}}     # [[low, high] sink1,   [low, high] sink2]
    min_TP::Float64

    function InternalCascade_2so2si(;
        duty::Vector{Float64} = [1.0, 1.0],
        source_temp::Vector{Vector{Float64}} = [[330.0, 320.0], [320.0, 300.0]],
        sink_temp::Vector{Vector{Float64}} = [[440.0, 450.0], [420.0, 440.0]],
        min_TP::Real = 1.0
    )
        length(duty)       == 2 || error("duty must have 2 elements: [duty_sink, duty_source]")
        length(source_temp) == 2 || error("source_temp must have 2 entries")
        length(sink_temp)  == 2 || error("sink_temp must have 2 entries")
        all(length(st) == 2 for st in source_temp) || error("each source_temp entry must be [high, low]")
        all(length(st) == 2 for st in sink_temp)   || error("each sink_temp entry must be [low, high]")
        new(duty, source_temp, sink_temp, Float64(min_TP))
    end
end

n_streams(::InternalCascade_2so2si) = 9

function build!(m::Model, config::InternalCascade_2so2si)
    prop = m[:prop]
    streams = total_streams(m, n_streams(config))
    (; duty, source_temp, sink_temp, min_TP) = config
    d1, d2 = duty[1], duty[2]
    total_d = d1 + d2

    total_same_pressure_idx(m, streams, [1, 2, 3, 4, 5])
    total_same_pressure_idx(m, streams, [6, 7, 8, 9])

    ϵ1 = @variable(m, lower_bound = -d1, upper_bound = d1)
    ϵ2 = @variable(m, lower_bound = -d1, upper_bound = d1)
    ϵ3 = @variable(m, lower_bound = -d2, upper_bound = d2)
    ϵ4 = @variable(m, lower_bound = -d2, upper_bound = d2)

    streams[1] = create_vap(m, prop, streams[1])
    same_stream(m, streams[1], streams[2])
    @constraint(m, streams[2].q_LH == d1 / total_d)
    separator(m, streams[2], streams[3], streams[4])

    streams[5] = create_liq(m, prop, streams[5])
    same_stream(m, streams[3], streams[5])

    same_stream(m, streams[5], streams[6])
    throttle = pressure_reduction(m, streams[6])
    streams[6] = create_liq(m, prop, streams[6])
    streams[7] = create_vap(m, prop, streams[7])
    same_stream(m, streams[6], streams[7])

    same_stream(m, streams[2], streams[8])
    same_stream(m, streams[8], streams[9])
    streams[9] = create_vap(m, prop, streams[9])
    @constraint(m, streams[8].q_LH == d1 / total_d)

    @constraints(m, begin
        streams[1].T >= sink_temp[1][2] + 5
        streams[2].T >= sink_temp[1][1] + 5
        streams[3].T >= sink_temp[2][2] + 5
        streams[5].T >= sink_temp[2][1] + 5
        streams[6].T <= source_temp[2][2] - 5
        streams[7].T <= source_temp[2][1] - 5
        streams[8].T <= source_temp[1][2] - 5
        streams[9].T <= source_temp[1][1] - 5
    end)

    @constraints(m, begin
        ϵ1 == simple_thermal_exergy(m, streams[1].T, streams[2].T, sink_temp[1][1], sink_temp[1][2], d1)
        ϵ2 == simple_thermal_exergy(m, streams[3].T, streams[5].T, sink_temp[2][1], sink_temp[2][2], d1)
        ϵ3 == simple_thermal_exergy(m, source_temp[1][1], source_temp[1][2], streams[8].T, streams[9].T, d2)
        ϵ4 == simple_thermal_exergy(m, source_temp[2][1], source_temp[2][2], streams[6].T, streams[7].T, d2)
    end)

    @constraint(m, throttle.P_low >= min_TP)

    m[:streams] = streams
    return ϵ1 + ϵ2 + ϵ3 + ϵ4
end

# ─── 3 Sources / 3 Sinks ──────────────────────────────────────────────────────

struct InternalCascade_3so3si <: HeatPumpConfig
    duty::Vector{Float64}                  # [duty_sink1, duty_sink2, duty_sink3]
    source_temp::Vector{Vector{Float64}}   # [[high, low] source3, source2, source1] (highest temp first)
    sink_temp::Vector{Vector{Float64}}     # [[low, high] sink1, sink2, sink3]
    min_TP::Float64

    function InternalCascade_3so3si(;
        duty::Vector{Float64} = [1.0, 1.0, 1.0],
        source_temp::Vector{Vector{Float64}} = [[330.0, 320.0], [320.0, 300.0], [310.0, 290.0]],
        sink_temp::Vector{Vector{Float64}} = [[440.0, 450.0], [420.0, 440.0], [400.0, 420.0]],
        min_TP::Real = 1.0
    )
        length(duty)       == 3 || error("duty must have 3 elements: [duty_sink1, duty_sink2, duty_sink3]")
        length(source_temp) == 3 || error("source_temp must have 3 entries")
        length(sink_temp)  == 3 || error("sink_temp must have 3 entries")
        all(length(st) == 2 for st in source_temp) || error("each source_temp entry must be [high, low]")
        all(length(st) == 2 for st in sink_temp)   || error("each sink_temp entry must be [low, high]")
        new(duty, source_temp, sink_temp, Float64(min_TP))
    end
end

n_streams(::InternalCascade_3so3si) = 14

function build!(m::Model, config::InternalCascade_3so3si)
    streams = total_streams(m, n_streams(config))
    (; duty, source_temp, sink_temp, min_TP) = config
    d1, d2, d3 = duty[1], duty[2], duty[3]
    total_d = d1 + d2 + d3

    total_same_pressure_idx(m, streams, [1, 2, 3, 4, 5, 6, 7, 8])
    total_same_pressure_idx(m, streams, [9, 10, 11, 12, 13, 14])
    all_vapor_streams(m, streams, [1, 10, 12, 14])
    all_liquid_streams(m, streams, [8, 9])

    ϵ = [@variable(m, lower_bound = -maximum(duty), upper_bound = maximum(duty)) for _ in 1:6]

    same_stream(m, streams[1],  streams[2])
    same_stream(m, streams[3],  streams[5])
    same_stream(m, streams[6],  streams[8])
    same_stream(m, streams[8],  streams[9])
    same_stream(m, streams[9],  streams[10])
    same_stream(m, streams[11], streams[5])
    same_stream(m, streams[11], streams[12])
    same_stream(m, streams[13], streams[2])
    same_stream(m, streams[13], streams[14])

    separator(m, streams[2], streams[3], streams[4])
    separator(m, streams[5], streams[6], streams[7])

    @constraint(m, streams[2].q_LH  == d1 / total_d)
    @constraint(m, streams[13].q_LH == streams[2].q_LH)
    @constraint(m, streams[5].q_LH  == d2 / (d2 + d3))
    @constraint(m, streams[11].q_LH == streams[5].q_LH)

    throttle = pressure_reduction(m, streams[9])
    @constraint(m, throttle.P_low >= min_TP)

    @constraints(m, begin
        streams[1].T  >= sink_temp[1][2] + 5
        streams[2].T  >= sink_temp[1][1] + 5
        streams[3].T  >= sink_temp[2][2] + 5
        streams[5].T  >= sink_temp[2][1] + 5
        streams[6].T  >= sink_temp[3][2] + 5
        streams[8].T  >= sink_temp[3][1] + 5
        streams[9].T  <= source_temp[3][2] - 5
        streams[10].T <= source_temp[3][1] - 5
        streams[11].T <= source_temp[2][2] - 5
        streams[12].T <= source_temp[2][1] - 5
        streams[13].T <= source_temp[1][2] - 5
        streams[14].T <= source_temp[1][1] - 5
    end)

    @constraints(m, begin
        ϵ[1] == simple_thermal_exergy(m, streams[1].T,  streams[2].T,  sink_temp[1][1], sink_temp[1][2], d1)
        ϵ[2] == simple_thermal_exergy(m, streams[3].T,  streams[5].T,  sink_temp[2][1], sink_temp[2][2], d2)
        ϵ[3] == simple_thermal_exergy(m, streams[6].T,  streams[8].T,  sink_temp[3][1], sink_temp[3][2], d3)
        ϵ[4] == simple_thermal_exergy(m, source_temp[3][1], source_temp[3][2], streams[9].T,  streams[10].T, d3)
        ϵ[5] == simple_thermal_exergy(m, source_temp[2][1], source_temp[2][2], streams[11].T, streams[12].T, d2)
        ϵ[6] == simple_thermal_exergy(m, source_temp[1][1], source_temp[1][2], streams[13].T, streams[14].T, d1)
    end)

    m[:streams] = streams
    return sum(ϵ)
end


# ─── 1 Sources / 2 Sinks (External Cascade) ──────────────────────────────────────────────────────

struct ExternalCascade_1so2si <: HeatPumpConfig
    duty::Vector{Float64}                  # [duty_sink1, duty_sink2]
    source_temp::Vector{Float64}           # [high, low]
    sink_temp::Vector{Vector{Float64}}     # [[low1, high1], [low2, high2]]
    min_TP::Float64

    function ExternalCascade_1so2si(;
        duty::Vector{Float64} = [0.5, 0.5],
        source_temp::Vector{Float64} = [300.0, 300.0],
        sink_temp::Vector{Vector{Float64}} = [[440.0, 450.0], [420.0, 440.0]],
        min_TP::Real = 1.0
    )
        length(duty)       == 2 || error("duty must have 2 elements: [duty_sink1, duty_sink2]")
        length(source_temp) == 2 || error("source_temp must have 2 elements: [high, low]")
        length(sink_temp)  == 2 || error("sink_temp must have 2 entries")
        all(length(st) == 2 for st in sink_temp) || error("each sink_temp entry must be [low, high]")
        new(duty, source_temp, sink_temp, Float64(min_TP))
    end
end

n_streams(::ExternalCascade_1so2si) = 6

function build!(m::Model, config::ExternalCascade_1so2si)
    prop = m[:prop]
    streams = total_streams(m, n_streams(config))
    (; duty, source_temp, sink_temp, min_TP) = config
    d1, d2 = duty[1], duty[2]
    total_d = d1 + d2

    total_same_pressure_idx(m, streams, [2, 3])
    total_same_pressure_idx(m, streams, [1, 4])
    total_same_pressure_idx(m, streams, [5, 6])
    all_vapor_streams(m, streams, [1, 2, 6])
    all_liquid_streams(m, streams, [3, 4, 5])

    ϵ = [@variable(m, lower_bound = -maximum(duty), upper_bound = maximum(duty)) for _ in 1:3]
    all_same_stream(m, streams, [1, 2, 3, 4, 5, 6])

    throttle = pressure_reduction(m, streams[5])
    @constraint(m, throttle.P_low >= min_TP)

    @constraints(m, begin
        streams[2].T >= sink_temp[1][2] + 5
        streams[3].T >= sink_temp[1][1] + 5
        streams[1].T >= sink_temp[2][2] + 5
        streams[4].T >= sink_temp[2][1] + 5
        streams[5].T <= source_temp[2] - 5
        streams[6].T <= source_temp[1] - 5
    end)

    @constraints(m, begin
        ϵ[1] == simple_thermal_exergy(m, streams[2].T, streams[3].T, sink_temp[1][1], sink_temp[1][2], d1)
        ϵ[2] == simple_thermal_exergy(m, streams[1].T, streams[4].T, sink_temp[2][1], sink_temp[2][2], d2)
        ϵ[3] == simple_thermal_exergy(m, source_temp[1], source_temp[2], streams[5].T, streams[6].T, total_d)
    end)

    m[:streams] = streams
    return sum(ϵ)
end


# ─── 1 Sources / 1 Sinks (n-Stage External Cascade) ──────────────────────────────────────────────────────

struct ExternalCascade_1so1si <: HeatPumpConfig
    duty::Float64
    source_temp::Vector{Float64}     # [high, low]
    sink_temp::Vector{Float64}       # [low, high]
    min_TP::Float64
    N::Int

    function ExternalCascade_1so1si(;
        duty::Real = 1.0,
        source_temp::Vector{Float64} = [300.0, 300.0],
        sink_temp::Vector{Float64} = [400.0, 400.0],
        min_TP::Real = 1.0,
        N::Int = 1
    )
        length(source_temp) == 2 || error("source_temp must have 2 elements: [high, low]")
        length(sink_temp)   == 2 || error("sink_temp must have 2 elements: [low, high]")
        new(Float64(duty), source_temp, sink_temp, Float64(min_TP), N)
    end
end

n_streams(config::ExternalCascade_1so1si) = 4 * config.N

function build!(m::Model, config::ExternalCascade_1so1si)
    prop = m[:prop]
    streams = total_streams(m, n_streams(config))
    (; duty, source_temp, sink_temp, min_TP, N) = config
    throttle = Vector{pressure_reduction}(undef, N)
    for i in 1:N
        same_pressure(m, streams[1 + 4 * (i - 1)], streams[2 + 4 * (i - 1)])
        same_pressure(m, streams[3 + 4 * (i - 1)], streams[4 * (i)])
        streams[1 + 4 * (i - 1)] = create_vap(m, prop, streams[1 + 4 * (i - 1)])
        streams[4 * i] = create_vap(m, prop, streams[4 * i])
        streams[2 + 4 * (i - 1)] = create_liq(m, prop, streams[2 + 4 * (i - 1)])
        streams[3 + 4 * (i - 1)] = create_liq(m, prop, streams[3 + 4 * (i - 1)])
        same_stream(m, streams[1 + 4 * (i - 1)], streams[2 + 4 * (i - 1)])
        same_stream(m, streams[2 + 4 * (i - 1)], streams[3 + 4 * (i - 1)])
        same_stream(m, streams[3 + 4 * (i - 1)], streams[4 * i])
        throttle[i] = pressure_reduction(m, streams[3 + 4 * (i - 1)])
        @constraint(m, throttle[i].P_low >= min_TP)
    end

    ϵ = [@variable(m, lower_bound = -duty, upper_bound = duty) for _ in 1:(N + 1)]

    @constraints(m, begin
        streams[1].T >= sink_temp[2] + 5
        streams[2].T >= sink_temp[1] + 5
        streams[4 * N - 1].T <= source_temp[2] - 5
        streams[4 * N].T <= source_temp[1] - 5
        [i = 2:N], streams[4 * (i - 1) + 1].T >= streams[4 * (i - 1)].T + 5
        [i = 2:N], streams[4 * (i - 1) + 2].T >= streams[4 * (i - 2) + 3].T + 5
    end)

    @constraints(m, begin
        ϵ[1] == simple_thermal_exergy(m, streams[1].T, streams[2].T, sink_temp[1], sink_temp[2], duty)
        [i = 2:N], ϵ[i] == simple_thermal_exergy(m, streams[4 * (i - 1) + 1].T, streams[4 * (i - 1) + 2].T, streams[4 * (i - 2) + 3].T, streams[4 * (i - 1)].T, duty)
        ϵ[N + 1] == simple_thermal_exergy(m, source_temp[1], source_temp[2], streams[4 * (N - 1) + 3].T, streams[4 * N].T, duty)
    end)

    m[:streams] = streams
    return sum(ϵ)
end




