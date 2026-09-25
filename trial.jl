using Revise
include("load_model.jl")

m = Model()
m[:prop] = Properties(mixture_property("Mixture_49"))

# config = ExternalCascade_1so1si(
#     duty        = 1.0,
#     source_temp = [68.0 + 273.0, 45.0 + 273.0],   # [high, low]
#     sink_temp   = [58.0 + 273.0, 85.6 + 273.0],   # [low, high]
#     min_TP      = 1.0,
#     N           = 1
# )

config = ExternalCascade_1so2si(
    duty        = [1.0, 1.0],
    source_temp = [68.0 + 273.0, 45.0 + 273.0],                     # [high, low]
    sink_temp   = [[58.0 + 50.0 + 273.0, 85.6 + 50.0 + 273.0],                     # sink 1 (lower temp): [low, high]
                   [58.0 + 273.0, 85.6 + 273.0]],       # sink 2 (higher temp): [low, high]
    min_TP      = 1.0
)


obj = build!(m, config)

optimize_model!(m, obj)

# 4 streams total: loop 1 → streams 1-4
for i in 1:n_streams(config)
    println("=== Stream $i ===")
    to_print(m[:streams][i])
    println()
end

fig = tq_plot(m, config)
display(fig)
# save("tq_external_1so1si_3.png", fig)
