
"""
    struct CParameters{T<:AbstractFloat}

Description
"""
m                 = readdlm("input_vectors.txt")
pd_vec      = m[: ,2]
e_vec0       = m[: ,3]
#@load "src\\results.jld2" mc vfi
pd_vec24 = zeros(12)
e_vec24 = zeros(12)
for i=1:12
    pd_vec24[i] = pd_vec[1+6*i-6]
    e_vec24[i] = e_vec0[1+6*i-6]
end
pd_vec= pd_vec24
e_vec0 = e_vec24
#
Kv_initial = [ 0.01958164631147754
             0.7520826034598947
             1.9456940599353418
             3.2198568969161316
             4.313816725882135
             5.118029008613562
             6.208340888874423
             6.326671675562466
             5.335424053447406
             4.1696216792756005
             2.6816085899936515
             0.7991155739521087]
Lv_initial = [ 0.4826031725043264
             0.5791238412539973
             0.6088865008592285
             0.5982649764550019
             0.5658970879502183
             0.5250512503385936
             0.45599227547990673
             0.283791593638212
             0.0023326853140481748
             0.00019902306429751116
             0.0002703079757560597
             6.285105305725038e-5]
Mv_initial = [ 0.20990128122439222
             0.2508361669400409
             0.2709010722724846
             0.3152901721321519
             0.3738116583487679
             0.4388172154143019
             0.5086232525157377
             0.5872670542382546
             0.6913191037172397
             0.6877433602613067
             0.6870268438507426
             0.8755549283549331]
#


@with_kw struct InitialGuess{F<:Float64}
    K::F =4.003851661979106
    L::F = 0.3925936870533928
    T::F = 0.06802416941178728
    τw::F = 0.3076344271799764
    Kv::Array{F} = Kv_initial
    Lv::Array{F} = Lv_initial
    Mv::Array{F} = Mv_initial
end

@with_kw struct Demographics{T<:Float64,I<:Int64}
    n_age::I = 12
    R_age::I = 7
    n_pop::T = (1+0.03)^1 -1
    NN::I = 10000
    pd_vec::Array{T} = pd_vec
end

demographics = Demographics()

@with_kw struct Preferences{T<:Float64}
    σ::T = 1.5
    ψ::T = -0.5
    γ::T = 0.45
    ν::T = 0.75
    θᵦ::T = 0.01
    κ::T = 0.0
    # Constant in the utility function
    β::T = 0.975^1
    δₕ::Array{T} = collect(range(0.1,stop=0.9,length = demographics.n_age))
    Q::Array{T} = 2*ones(demographics.n_age)
    z::Array{T} = 0.6*ones(demographics.n_age)

end

preferences = Preferences()

@with_kw struct Grids{I<:Int64,T<:Float64}

    #Assets
    n_a::I = 30
    a_min::T = 0.0
    a_max::T = 50.0
    a_scale::T =3.0
    a_grid::Array{T} = expspace(a_min,a_max,a_scale,n_a)

    # Health
    n_h::I = 20
    h_min::T = 0.1
    h_max::T = 5.0
    h_grid::Array{T} = collect(range(h_min,stop=h_max,length=n_h))


    # Health Shock
    n_zh::I = 1
    ρₕ::T = 0.86
    σₕ::T = 0.045
    μₕ::T = 2.5
    Nₕ::I = 2
    zh::Array{T} = collect(tauchen(n_zh, ρₕ, σₕ,0.0, Nₕ)[1])
    Pzh::Array{T} = tauchen(n_zh, ρₕ, σₕ,0.0, Nₕ)[2]
    Pzh_CDF::Array{T} = tauchen(n_zh, ρₕ, σₕ,0.0, Nₕ)[3]
    # Initial health endowment
    zh1::Array{T} = collect(tauchen(n_zh, 0.0, σₕ,μₕ, Nₕ)[1])
    Pzh1::Array{T} = tauchen(n_zh, 0.0, σₕ,μₕ, Nₕ)[2]
    Pzh1_CDF::Array{T} = tauchen(n_zh, 0.0, σₕ,0.0, Nₕ)[3]


    # Productivity Shock
    n_zl::I = 1
    ρₗ::T = 0.96
    σₗ::T = 0.045
    μₗ::T = 0.0
    Nₗ::I = 2
    zl::Array{T} = collect(tauchen(n_zl, ρₗ, σₗ,0.0, Nₗ)[1])
    Pzl::Array{T} = tauchen(n_zl, ρₗ, σₗ,0.0, Nₗ)[2]
    Pzl_CDF::Array{T} = tauchen(n_zl, ρₗ, σₗ,0.0, Nₗ)[3]

    zl1::Array{T} = collect(tauchen(n_zl, 0.0, σₗ,μₗ, Nₗ)[1])
    Pzl1::Array{T} = tauchen(n_zl, 0.0, σₗ,μₗ, Nₗ)[2]
    Pzl1_CDF::Array{T} = tauchen(n_zl, 0.0, σₗ,μₗ, Nₗ)[3]
    e_vec::Array{T} = e_vec0
end

grids = Grids()



@with_kw struct Government{T<:Float64}
    g_ratio::T = 0.194
    ss_rate::T = 0.335
    τc::T = 0.05
    τa::T = 0.36
    c̲::T = 0.01
end

government = Government()

@with_kw struct Firms{T<:Float64}
    A::T = 5.0
    α::T = 0.36
    δ::T = (1+0.06)^3 - 1
end

firms = Firms()

@with_kw struct Options{I<:Int64}
    SEPARABLE::I = 1
    CALIBRATE::I = 0
end

options = Options()

@with_kw struct Params
    D::Demographics = demographics
    P::Preferences = preferences
    grid::Grids = grids
    G::Government = government
    F::Firms = firms
    O::Options = options
end


function Base.show(io::IO, n::Options)
    printstyled("----    Options     -----\n",color=:red)
    printstyled("Separable Utility:      \t|",color=:green)
    println(io, " SEPARABLE = $(n.SEPARABLE)")
    printstyled("Calibration:            \t|",color=:green)
    println(io, " CALIBRATE = $(n.CALIBRATE)")
    return
end

function Base.show(io::IO, n::Demographics)
    printstyled("---- Demographics   -----\n",color=:red)
    printstyled("Max Age:                \t|",color=:green)
    println(io, " n_age = $(n.n_age)")
    printstyled("Retirement Age:         \t|",color=:green)
    println(io, " R_age = $(n.R_age)")
    printstyled("Population Growth:      \t|",color=:green)
    println(io, " n_pop = $(n.n_pop)")
    printstyled("N⁰ Simulations:         \t|",color=:green)
    println(io, " NN = $(n.NN)")
    return
end

function Base.show(io::IO, n::Firms)
    printstyled("----     Firms      -----\n",color=:red)
    printstyled("TFP:                    \t|",color=:green)
    println(io, " A = $(n.A)")
    printstyled("Share of Capital:       \t|",color=:green)
    println(io, " α = $(n.α)")
    printstyled("Capital Depreciation:   \t|",color=:green)
    println(io, " δ = $(n.δ)")
    return
end

function Base.show(io::IO, n::Government)
    printstyled("----   Government   -----\n",color=:red)
    printstyled("Spending Ratio:         \t|",color=:green)
    println(io, " g_ratio= $(n.g_ratio)")
    printstyled("Replacement Rate:       \t|",color=:green)
    println(io, " ss_rate = $(n.ss_rate)")
    printstyled("Consumption Tax:        \t|",color=:green)
    println(io, " τc = $(n.τc)")
    printstyled("Capital Income Tax:     \t|",color=:green)
    println(io, " τa = $(n.τa)")
    return
end


function Base.show(io::IO, n::Preferences)
    printstyled("----  Preferences   -----\n",color=:red)
    printstyled("CRRA Consumption:       \t|",color=:green)
    println(io, " ψ = $(n.ψ)")
    printstyled("CRRA Leisure:           \t|",color=:green)
    println(io, " γ = $(n.γ)")
    printstyled("Discount Factor:        \t|",color=:green)
    println(io, " β = $(n.β)")
    return
end

function Base.show(io::IO, n::Grids)
    printstyled("----    Grids       -----\n",color=:red)
    printstyled("----    Assets      -----\n",color=:yellow)
    printstyled("N⁰ Grid Points:         \t|",color=:green)
    println(io, " n_a = $(n.n_a)")
    printstyled("Min:                    \t|",color=:green)
    println(io, " a_min = $(n.a_min)")
    printstyled("Max:                    \t|",color=:green)
    println(io, " a_max = $(n.a_max)")
    printstyled("Scale:                  \t|",color=:green)
    println(io, " a_scale = $(n.a_scale)")
    printstyled("--  Productivity Shock --\n",color=:yellow)
    printstyled("N⁰ Grid Points:         \t|",color=:green)
    println(io, " n_zl = $(n.n_zl)")
    printstyled("Autocorrelation:        \t|",color=:green)
    println(io, " ρₗ = $(n.ρₗ)")
    printstyled("Variance:               \t|",color=:green)
    println(io, " σₗ = $(n.σₗ)")
    return
end

function printparameters(P,G,Grid,F,D,OPTION)
    printstyled("-------------------------\n",color=:blue)
    printstyled("------  Description -----\n",color=:blue)
    printstyled("-------------------------\n",color=:blue)
    println(D)
    println(P)
    println(F)
    println(G)
    println(Grid)
    println(OPTION)
    printstyled("-------------------------\n",color=:blue)
    printstyled("-------------------------\n",color=:blue)
    printstyled("-------------------------\n",color=:blue)
end


# Print the parameters
#printparameters(demographics,preferences,grids,government,firms,options)
