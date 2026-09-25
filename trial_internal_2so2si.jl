using Revise
include("load_model.jl")

m = Model()
m[:prop] = Properties(mixture_property("Mixture_33"))

config = InternalCascade_2so2si(
    duty        = [1.3, 2.8],                       # [duty_sink, duty_source]
    source_temp = [[80.0 + 273.0, 75.0 + 273.0], [75.0 + 273.0, 70.0 + 273.0]], # [[high, low] source2, [high, low] source1]
    sink_temp   = [[186.0 + 273.0, 224.0 + 273.0], [120.0 + 273.0, 164.0 + 273.0]], # [[low, high] sink1, [low, high] sink2]
    min_TP      = 1.0
)

obj = build!(m, config)

optimize_model!(m, obj)

# 9 streams total
for i in 1:n_streams(config)
    println("=== Stream $i ===")
    to_print(m[:streams][i])
    println()
end

fig = tq_plot(m, config)
# save("tq_internal_2so2si.png", fig)
