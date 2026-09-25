using Revise
include("load_model.jl")

m = Model()
m[:prop] = Properties(mixture_property("Mixture_30"))

config = InternalCascade_1so2si_new(
    duty        = [1.3, 2.8],                                     # [duty_sink1, duty_sink2]
    source_temp = [80.0 + 273.0, 70.0 + 273.0],                   # [high, low]
    sink_temp   = [[186.0 + 273.0, 224.0 + 273.0], [120.0 + 273.0, 164.0 + 273.0]], # [[low1, high1], [low2, high2]]
    min_TP      = 1.0
)

obj = build!(m, config)

optimize_model!(m, obj)

# 10 streams total
for i in 1:n_streams(config)
    println("=== Stream $i ===")
    to_print(m[:streams][i])
    println()
end

fig = tq_plot(m, config)
# save("tq_internal_1so2si.png", fig)
