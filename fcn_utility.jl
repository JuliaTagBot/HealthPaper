"""
    function calc_utility(c::Float64, l::Float64, a′, P::Params)

Description
"""
function calc_utility(c::Float64, l::Float64, h::Float64, P::Params)
    γ = P.P.γ
    σ = P.P.σ
    ν = P.P.ν
    ψ = P.P.ψ
    les = 1 - l
    u = ν*(c^γ*les^(1-γ))^ψ
    u = u + (1-ν)*h^ψ
    u = u^((1-σ)/ψ)/(1-σ)
    return u
end

"""
    function calc_prob_survival(h′, age)

Description
"""
function calc_prob_survival(age::Int64,h::Float64,P::Params)
    n_age = P.D.n_age
    pd_vec = P.D.pd_vec

    if age == 0
        p_surv = 1.0
    elseif age >= n_age
        p_surv = 0.0
    else
        dh = (1-pd_vec[age])*0.5*exp(-0.24*h)
        p_surv = pd_vec[age] - dh
    end
    return p_surv
end
"""
    function calc_warm_glow(a′, P)

Description
"""
function calc_warm_glow(a′::Float64,P::Params)
    θᵦ = P.P.θᵦ
    κ = P.P.κ
    σ = P.P.σ

    b = θᵦ*(a′+κ)^(1-σ)/(1-σ)
end
"""
    function calc_warm_glow(a′, P)

Description
"""
function calc_current_health(m::Float64,h::Float64,zh::Float64,age::Int64,P::Params)
    δₕ = P.P.δₕ
    Q = P.P.Q
    z = P.P.z

    h′ = (1-δₕ[age])*h*exp(zh) + Q[age]*m^z[age]
end
