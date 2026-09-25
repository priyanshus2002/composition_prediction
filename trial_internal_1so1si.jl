using Revise
include("load_model.jl")

m = Model()
m[:prop] = Properties(mixture_property("Mixture_48"))

config = InternalCascade_1so1si(
    duty        = 1.0,
    source_temp = [363.0, 326.0],   # [high, low]
    sink_temp   = [451.0, 463.0],   # [low, high]
    min_TP      = 1.0
)



obj = build!(m, config)

# Fix mole fraction of stream 1 
# JuMP.fix(m[:streams][1].z[1], 0.2; force = true)
# JuMP.fix(m[:streams][1].z[2], 0.2; force = true)
# JuMP.fix(m[:streams][1].z[3], 0.1; force = true)
# # JuMP.fix(m[:streams][1].z[4], 0.4; force = true)
# # JuMP.fix(m[:streams][6].P, 1.1; force = true)
# JuMP.fix(m[:streams][1].P, 13.0; force = true)

optimize_model!(m, obj)

# 9 streams total
for i in 1:n_streams(config)
    println("=== Stream $i ===")
    to_print(m[:streams][i])
    println()
end

fig = tq_plot(m, config)
# save("tq_internal_1so1si.png", fig)
