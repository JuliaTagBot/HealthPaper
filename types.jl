struct VFIResults{T<:AbstractFloat}
    v::Array{T}
    g_c::Array{T}
    g_l::Array{T}
    g_a′::Array{T}
    g_m::Array{T}
    g_h′::Array{T}
end

struct MonteCarloResults{T<:AbstractFloat}
    K::T
    L::T
    C::T
    M::T
    Le::T
    T::T
    Tw::T
    τw::T
    μ::Array{T}
    ψ::Array{T}
    Kv::Array{T}
    Lv::Array{T}
    Lev::Array{T}
    Cv::Array{T}
    Mv::Array{T}
    Twv::Array{T}
    Bqv::Array{T}
    sim_a′::Array{T}
    sim_c::Array{T}
    sim_m::Array{T}
    sim_h′::Array{T}
    sim_l::Array{T}
    sim_le::Array{T}
    sim_wt::Array{T}
end

struct Solution{T<:AbstractFloat}
    v::T
    l::T
    m::T
    a′::T
    c::T
    h′::T
end


function Base.show(io::IO, n::Solution)
    print(io, "val =   $(n.v)\n")
    print(io, "m   =   $(n.m)\n")
    print(io, "l   =   $(n.l)\n")
    print(io, "a′  =   $(n.a′)\n")
    print(io, "c   =   $(n.c)\n")
    print(io, "h   =   $(n.h′)\n")
    return
end
