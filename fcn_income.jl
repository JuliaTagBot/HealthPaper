"""
    function calc_non_labour_income(τa::Float64,r::Float64,a_grid::Array{Float64},i_a::Int64,T::Float64,b::Float64)

Description
"""
function calc_non_labour_income(τa::Float64,r::Float64,a::Float64,T::Float64,b::Float64)
    y0 = (1.0 + (1.0 - τa)*r)*a + T + b
end

"""
    function calc_initial_guess_l(age, w, y0, τw, m, a′)

Description
"""
function calc_initial_guess_l(zl::Float64,h::Float64,age::Int64, w::Float64, y0::Float64, τw::Float64, a′::Float64, P::Params)
    # compute minimum labor productivity
    e = calc_labour_productivity(zl,h,age, P) # e if h = min
    # set n to ensure that it's feasible for (c,f) to be positive
    n0 = 1.50*(a′ - y0) / ( (1.0-τw)*w*e )
    n0 = max( n0, 0.01)  # ensure that n0 is strictly positive
    n0 = min(n0,0.99)
    return n0
end
"""
    function calc_labor_productivity(zl, age, e_vec)

Description
"""
function calc_labour_productivity(zl::Float64,h::Float64,age::Int64, P::Params)
    e_vec = P.grid.e_vec

    el = e_vec[age]*exp(zl+0.001*(h-5))

    return el
end

"""
    function calc_labor_income(zl, age, e_vec)

Description
"""
function calc_labour_income(zl::Float64,h::Float64, age::Int64, w::Float64, l::Float64,τw::Float64,P::Params)
    el = calc_labour_productivity(zl, h, age, P)
    lab_inc = (1-τw)*l*w*el
end

"""
    function calc_transfers(y, P)

Description
"""
function calc_transfers(y::Float64,m::Float64,P::Params)
        c̲ = P.G.c̲

        wt = max(0,c̲ - y)
end
"""
    function calc_net_inc(i_a::Int64, zl::Float64, age::Int64, w::Float64, y0::Float64, l::Float64, a′::Float64,τw::Float64,P::Params)

Description
"""
function calc_net_inc(zl::Float64,h::Float64, age::Int64, w::Float64, y0::Float64, l::Float64,m::Float64, a′::Float64,τw::Float64,P::Params)
    lab_inc =  calc_labour_income(zl, h, age, w, l, τw, P)
    y = y0 + lab_inc
    wt = calc_transfers(y,m,P)
    y = y + wt
end
"""
    function calc_other_inputs(i_a, i_h, i_zl, i_zh, age,w, y0, m_sol, l_sol, a′_sol,τw)

Description
"""
function calc_other_inputs(i_a::Int64, i_h, i_zl::Int64,i_zh::Int64, age::Int64, w::Float64, y0::Float64, l::Float64,m::Float64, a′::Float64,τw::Float64,P::Params)
    τc = P.G.τc
    zl = P.grid.zl[i_zl]
    zh = P.grid.zh[i_zh]
    c̲ = P.G.c̲
    h = P.grid.h_grid[i_h]
    h′ = calc_current_health(m,h,zh,age,P)
    y = calc_net_inc( zl,h′, age, w, y0, l,m,a′,τw,P)

    c = (y - a′-m)/(1+τc)

    return c, h′
end
