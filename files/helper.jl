using JuMP
using Gurobi

Properties(mix::MixtureProperties) =
    Properties(α=mix.α, λ=mix.λ, λ_P=mix.λ_P, antoine=mix.antoine, pRef=mix.pRef, n=mix.n)

function optimize_model!(
    m::Model,
    objective;
    sense = MOI.MIN_SENSE,
    time_limit::Real = 600,
    mip_gap::Real = 0.005
)
    set_optimizer(m, Gurobi.Optimizer)
    set_optimizer_attribute(m, "TimeLimit", time_limit)
    set_optimizer_attribute(m, "MIPGap", mip_gap)
    set_objective(m, sense, objective)
    optimize!(m)
    return m
end

