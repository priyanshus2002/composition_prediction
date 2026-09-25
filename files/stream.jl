using JuMP

struct Properties
    α::Vector{<:Number}
    λ::Vector{<:Number}
    λ_P::Vector{Vector{<:Number}}
    antoine::Vector{<:Number}
    pRef::Number
    n::Int


    function Properties(; 
        α::Vector{<:Number},
        antoine::Vector{<:Number},
        λ::Vector{Float64},
        λ_P::Vector{Vector{Float64}},
        pRef::Number = 1.0,
        n::Int = 2
        )

        if λ_P ==[]
            λ_P = [[λ[i], 0] for i in 1:n]
        end

        if length(α) != n
            error("Length of α must be equal to n")
        end

        if length(antoine) != 3
            error("Length of antoine must be 3")
        end

        if length(λ) != n
            error("Length of λ must be equal to n")
        end

        if length(λ_P) != n
            error("Length of λ_P must be equal to n")
        end

        if any(x -> length(x) != 2, λ_P)
            error("Each element of λ_P must be a vector of length 2")
        end

        return new(α, λ, λ_P, antoine, pRef, n)
    end
end

function _stream_create(
    m::Model,
    prop::Properties
)
    A, B, C = prop.antoine
    α = prop.α
    λ = prop.λ
    λ_P = prop.λ_P
    pRef = prop.pRef
    n = prop.n

    T_func(ρ) = B / (A + log(ρ)) - C

    x        = @variable(m, [1:n], lower_bound = 0, upper_bound = 1)
    y        = @variable(m, [1:n], lower_bound = 0, upper_bound = 1)
    z        = @variable(m, [1:n], lower_bound = 0, upper_bound = 1)
    q        = @variable(m, lower_bound = 0, upper_bound = 1)
    q_LH     = @variable(m, lower_bound = 0, upper_bound = 1)
    β        = @variable(m, lower_bound = minimum(α), upper_bound = maximum(α))
    ρ        = @variable(m, lower_bound = minimum(α), upper_bound = pRef * maximum(α) / 0.5)
    T        = @variable(m, lower_bound = T_func(pRef * maximum(α) / 0.5), upper_bound = T_func(minimum(α)))
    mean_λ   = @variable(m, lower_bound = 0, upper_bound = maximum(λ_P)[1])
    mean_λ_x = @variable(m, lower_bound = 0, upper_bound = maximum(λ_P)[1])
    r_1      = @variable(m, [1:n], lower_bound = 0, upper_bound = maximum(α) / minimum(α))
    r_2      = @variable(m, [1:n], lower_bound = 0, upper_bound = maximum(α) / minimum(α))
    bq       = @variable(m, lower_bound = 0, upper_bound = maximum(α))
    P        = @variable(m, lower_bound = 0.1, upper_bound = pRef)

    @constraints(m, begin
        sum(x[i] for i in 1:n) == 1
        sum(y[i] for i in 1:n) == 1
        sum(z[i] for i in 1:n) == 1

        [i = 1:(n - 1)], z[i] == q * x[i] + (1 - q) * y[i]

        mean_λ == sum((λ_P[i][1] + λ_P[i][2] * P) * z[i] for i in 1:n)
        mean_λ_x == sum((λ_P[i][1] + λ_P[i][2] * P) * x[i] for i in 1:n)

        q_LH * mean_λ == q * mean_λ_x

        sum(r_1[i] for i in 1:n) == 1
        sum(r_2[i] for i in 1:n) == 1
        bq == β * q
        [i = 1:n], α[i] * z[i] == r_1[i] * (bq + (1 - q) * α[i])
        [i = 1:n], β * z[i] == r_2[i] * (bq + (1 - q) * α[i])

        [j = 1:(n)], y[j] * β == α[j] * x[j]

        ρ * P == pRef * β

        (T + C) * (A + log(ρ)) == B
    end)

    return (x, y, z, q, q_LH, mean_λ, mean_λ_x, r_1, r_2, bq, β, ρ, T, P)
end

struct create_stream
    x::Vector{VariableRef}
    y::Vector{VariableRef}
    z::Vector{VariableRef}
    q::VariableRef
    q_LH::VariableRef
    mean_λ::VariableRef
    mean_λ_x::VariableRef
    r_1::Vector{VariableRef}
    r_2::Vector{VariableRef}
    bq::VariableRef
    β::VariableRef
    ρ::VariableRef
    T::VariableRef
    P::VariableRef

    function create_stream(m::Model)
        prop = m[:prop]
        data = _stream_create(m, prop)
        return new(data...)
    end
end





