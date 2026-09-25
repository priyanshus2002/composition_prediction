using CairoMakie

struct HXData
    label::String
    T_hot::Vector{Float64}   # [T_in, T_out] — high to low (hot side enters hot)
    T_cold::Vector{Float64}  # [T_in, T_out] — low to high (cold side enters cold)
    Q::Float64
end

# ── Data extraction (one method per config) ────────────────────────────────────

function tq_data(m::Model, config::InternalCascade_1so1si)
    s = m[:streams]
    (; duty, source_temp, sink_temp) = config
    return [
        HXData("Sink Condenser",
            [value(s[1].T), value(s[2].T)],
            [sink_temp[1], sink_temp[2]],
            duty),
        HXData("Source Evaporator",
            [source_temp[1], source_temp[2]],
            [value(s[6].T), value(s[7].T)],
            duty),
        HXData("Cascade HX",
            [value(s[3].T), value(s[5].T)],
            [value(s[8].T), value(s[9].T)],
            duty),
    ]
end

function tq_data(m::Model, config::InternalCascade_1so2si)
    s = m[:streams]
    (; duty, source_temp, sink_temp) = config
    d1, d2 = duty[1], duty[2]
    return [
        HXData("Sink 1 Condenser",
            [value(s[1].T), value(s[2].T)],
            [sink_temp[1][1], sink_temp[1][2]],
            d1),
        HXData("Sink 2 Condenser",
            [value(s[3].T), value(s[5].T)],
            [sink_temp[2][1], sink_temp[2][2]],
            d2),
        HXData("Source Evaporator",
            [source_temp[1], source_temp[2]],
            [value(s[7].T), value(s[8].T)],
            d1 + d2),
        HXData("Cascade HX",
            [value(s[5].T), value(s[6].T)],
            [value(s[9].T), value(s[10].T)],
            d1),
    ]
end

function tq_data(m::Model, config::InternalCascade_1so2si_new)
    s = m[:streams]
    (; duty, source_temp, sink_temp) = config
    d1, d2 = duty[1], duty[2]
    return [
        HXData("Sink 1 Condenser",
            [value(s[1].T), value(s[2].T)],
            [sink_temp[1][1], sink_temp[1][2]],
            d1),
        HXData("Sink 2 Condenser",
            [value(s[5].T), value(s[6].T)],
            [sink_temp[2][1], sink_temp[2][2]],
            d2),
        HXData("Source Evaporator",
            [source_temp[1], source_temp[2]],
            [value(s[7].T), value(s[8].T)],
            d1 + d2),
        HXData("Cascade HX",
            [value(s[3].T), value(s[5].T)],
            [value(s[9].T), value(s[10].T)],
            d1),
    ]
end

function tq_data(m::Model, config::InternalCascade_2so2si)
    s = m[:streams]
    (; duty, source_temp, sink_temp) = config
    d1, d2 = duty[1], duty[2]
    return [
        HXData("Sink 1 Condenser",
            [value(s[1].T), value(s[2].T)],
            [sink_temp[1][1], sink_temp[1][2]],
            d1),
        HXData("Sink 2 Condenser",
            [value(s[3].T), value(s[5].T)],
            [sink_temp[2][1], sink_temp[2][2]],
            d1),
        HXData("Source 1 Evaporator",
            [source_temp[1][1], source_temp[1][2]],
            [value(s[8].T), value(s[9].T)],
            d2),
        HXData("Source 2 Evaporator",
            [source_temp[2][1], source_temp[2][2]],
            [value(s[6].T), value(s[7].T)],
            d2),
    ]
end

function tq_data(m::Model, config::InternalCascade_3so3si)
    s = m[:streams]
    (; duty, source_temp, sink_temp) = config
    d1, d2, d3 = duty[1], duty[2], duty[3]
    return [
        HXData("Sink 1 Condenser",
            [value(s[1].T),  value(s[2].T)],
            [sink_temp[1][1], sink_temp[1][2]],
            d1),
        HXData("Sink 2 Condenser",
            [value(s[3].T),  value(s[5].T)],
            [sink_temp[2][1], sink_temp[2][2]],
            d2),
        HXData("Sink 3 Condenser",
            [value(s[6].T),  value(s[8].T)],
            [sink_temp[3][1], sink_temp[3][2]],
            d3),
        HXData("Source 3 Evaporator",
            [source_temp[3][1], source_temp[3][2]],
            [value(s[9].T),  value(s[10].T)],
            d3),
        HXData("Source 2 Evaporator",
            [source_temp[2][1], source_temp[2][2]],
            [value(s[11].T), value(s[12].T)],
            d2),
        HXData("Source 1 Evaporator",
            [source_temp[1][1], source_temp[1][2]],
            [value(s[13].T), value(s[14].T)],
            d1),
    ]
end

function tq_data(m::Model, config::ExternalCascade_1so2si)
    s = m[:streams]
    (; duty, source_temp, sink_temp) = config
    d1, d2 = duty[1], duty[2]
    return [
        HXData("Sink 1 Condenser",
            [value(s[2].T), value(s[3].T)],
            [sink_temp[1][1], sink_temp[1][2]],
            d1),
        HXData("Sink 2 Condenser",
            [value(s[1].T), value(s[4].T)],
            [sink_temp[2][1], sink_temp[2][2]],
            d2),
        HXData("Source Evaporator",
            [source_temp[1], source_temp[2]],
            [value(s[5].T), value(s[6].T)],
            d1 + d2),
    ]
end

function tq_data(m::Model, config::ExternalCascade_1so1si)
    s = m[:streams]
    (; duty, source_temp, sink_temp, N) = config
    hxs = HXData[]
    push!(hxs, HXData("Stage 1 Condenser (Sink)",
        [value(s[1].T), value(s[2].T)],
        [sink_temp[1], sink_temp[2]],
        duty))
    for i in 2:N
        push!(hxs, HXData("Stage $i Cascade HX",
            [value(s[4*(i-1)+1].T), value(s[4*(i-1)+2].T)],
            [value(s[4*(i-2)+3].T), value(s[4*(i-1)].T)],
            duty))
    end
    push!(hxs, HXData("Stage $N Evaporator (Source)",
        [source_temp[1], source_temp[2]],
        [value(s[4*(N-1)+3].T), value(s[4*N].T)],
        duty))
    return hxs
end

# ── Plot ───────────────────────────────────────────────────────────────────────

function tq_plot(m::Model, config::HeatPumpConfig)
    hxs = tq_data(m, config)
    n = length(hxs)
    ncols = min(n, 3)
    nrows = cld(n, ncols)

    fig = Figure(size = (420 * ncols, 360 * nrows))

    for (k, hx) in enumerate(hxs)
        row = cld(k, ncols)
        col = mod1(k, ncols)
        ax = Axis(fig[row, col];
            title  = hx.label,
            xlabel = "Heat Duty",
            ylabel = "Temperature (K)")
        x = [0.0, hx.Q]
        T_hot = [hx.T_hot[2], hx.T_hot[1]]
        T_cold = [hx.T_cold[1], hx.T_cold[2]]

        # Shade the region between the hot and cold curves
        band!(ax, x, T_cold, T_hot; color = (:orange, 0.2))

        # Hot side: enters at Q=duty (high T), exits at Q=0 (low T)
        lines!(ax, x, T_hot; color = :red, linewidth = 2, label = "Hot side")
        # Cold side: enters at Q=0 (low T), exits at Q=duty (high T)
        lines!(ax, x, T_cold; color = :blue, linewidth = 2, label = "Cold side")
        axislegend(ax; position = :rc)
    end

    return fig
end
