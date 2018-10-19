"""
    function VFI(w, r, T, ss, τw, P)

Description
"""
function VFI(IG::InitialGuess,w::Float64, r::Float64, T::Float64, ss::Float64, τw::Float64,P::Params)

    n_a    = P.grid.n_a
    n_h    = P.grid.n_h
    n_zl   = P.grid.n_zl
    n_zh = P.grid.n_zh
    a_grid = P.grid.a_grid
    n_age  = P.D.n_age
    R_age  = P.D.R_age
    τa     = P.G.τa

    v0   = SharedArray(zeros(n_a, n_h, n_zl, n_zh))

    v    =  SharedArray(zeros(n_a, n_h, n_zl, n_zh, n_age))
    g_c  =  SharedArray(zeros(n_a, n_h, n_zl, n_zh, n_age))
    g_l  =  SharedArray(zeros(n_a, n_h, n_zl, n_zh, n_age))
    g_a′ =  SharedArray(zeros(n_a, n_h, n_zl, n_zh, n_age))
    g_m  =  SharedArray(zeros(n_a, n_h, n_zl, n_zh, n_age))
    g_h′ =  SharedArray(zeros(n_a, n_h, n_zl, n_zh, n_age))

    for age = n_age:-1:1

        if age<R_age
            b = 0.0
        else
            b = ss
        end

        #Threads.@threads

        @sync begin
              @distributed for ind = 1:(n_a*n_h)

              i_a      = convert(Int, ceil(ind/n_h));
              i_h      = convert(Int, floor(mod(ind-0.05, n_h))+1);
                @fastmath @inbounds for i_zl = 1:n_zl
                    for i_zh = 1:n_zh
                        y0 = calc_non_labour_income(τa,r,a_grid[i_a],T,b)

                        solution = solve(IG,i_a, i_h, i_zl, i_zh, age, Array(v0), w, τw, y0, P)

                        v[i_a, i_h, i_zl, i_zh, age]    = solution.v
                        g_c[i_a, i_h, i_zl, i_zh, age]  = solution.c
                        g_l[i_a, i_h, i_zl, i_zh, age]  = solution.l
                        g_a′[i_a, i_h, i_zl, i_zh, age] = solution.a′
                        g_m[i_a, i_h, i_zl, i_zh, age]  = solution.m
                        g_h′[i_a, i_h, i_zl, i_zh, age] = solution.h′
                    end
                end
            end
        v0 =  SharedArray(copy(v[:,:,:,:,age]))
    end
    println("Age = ",age)
    end
    return VFIResults(Array(v), Array(g_c), Array(g_l), Array(g_a′), Array(g_m), Array(g_h′))
end
"""
function solve(i_a::Int64, i_zl::Int64, age::Int64, v0::Array{Float64}, w::Float64, τw::Float64, y0::Float64,
    P::Params)

Description
"""
function solve(IG::InitialGuess,i_a::Int64, i_h::Int64, i_zl::Int64, i_zh::Int64 ,age::Int64, v0::Array{Float64}, w::Float64, τw::Float64, y0::Float64, P::Params)

    n_age   = P.D.n_age
    a_min   = P.grid.a_min
    a_max   = P.grid.a_max
    h_min   = P.grid.h_min
    h_max   = P.grid.h_max
    zl      = P.grid.zl[i_zl]
    zh      = P.grid.zh[i_zh]
    guess_a = IG.Kv[age]
    guess_l = IG.Lv[age]
    guess_m = IG.Mv[age]

    func(x) = solver_function(x, i_a, i_h, i_zl, i_zh, age, w, τw, y0, v0, a′, P)

    if age == n_age
        a′_min = 0.0
    else
        a′_min = a_min
    end
    a′_max = a_max

    a′ = 0.01*(a′_min + a′_max)
    #l = calc_initial_guess_l(zl,age, w, y0, τw, a′,P)
    x2 = zeros(3)
    x2[1] = guess_l
    x2[2] = guess_m
    x2[3] = guess_a
    res_val = Optim.optimize(func,x2,Optim.Options(iterations=100))
    val2 = -Optim.minimum(res_val)

    while val2<-1.0e37
        x3      = x2
        x3[1]   = (1+randn())*x2[1]
        x3[2]   = (1+randn())*x2[2]
        x3[3]   = (1+randn())*x2[3]
        res_val = Optim.optimize(func,x3,Optim.Options(iterations=50))
        val2    = -Optim.minimum(res_val)
    end

    values2 = Optim.minimizer(res_val)
    l_sol   = values2[1]
    m_sol   = values2[2]
    a′_sol  = values2[3]
    val     = val2


    c_sol, h′_sol = calc_other_inputs(i_a, i_h, i_zl, i_zh, age,w, y0, l_sol, m_sol, a′_sol,τw, P)
    return Solution(val, l_sol, m_sol, a′_sol, c_sol, h′_sol)
end

"""
    function solver_function(x,i_a, i_h, i_zl, i_zh, age, w, τw, y0,v0,a′)

Description
"""
function solver_function(x, i_a::Int64, i_h::Int64, i_zl::Int64, i_zh::Int64, age::Int64, w::Float64, τw::Float64, y0::Float64, v0::Array{Float64}, a′::Float64, P::Params)
    l  = x[1]
    m  = x[2]
    a′ = x[3]
    y  = - calc_val(l, m, a′, i_a, i_h, i_zl, i_zh, age, w, τw, y0, v0, P)
    return y
end
"""
    function calc_val(m, l, a′, i_a, i_h, i_zl, i_zh, age, w, τw, y0, v0)

Description
"""
function calc_val(l::Float64, m::Float64, a′::Float64, i_a::Int64, i_h::Int64, i_zl::Int64, i_zh::Int64, age::Int64, w::Float64, τw::Float64, y0::Float64, v0::Array{Float64},P::Params)

    if m>0
        τc  = P.G.τc
        Pzl = P.grid.Pzl
        Pzh = P.grid.Pzh
        β   = P.P.β
        a   = P.grid.a_grid[i_a]
        h   = P.grid.h_grid[i_h]
        h_max = P.grid.h_max
        h_min = P.grid.h_min
        zl  = P.grid.zl[i_zl]
        zh =  P.grid.zh[i_zh]
        h′ = calc_current_health(m,h,zh,age,P)
        les = 1 - l
        y   = calc_net_inc(zl, h′,age, w, y0, l, m, a′,τw,P)
        c   = (y - a′- m)/(1+τc)
        if c>0.0 && les>0.0 && les<=1.0 && h′<=h_max && h′>=h_min && a′>=P.grid.a_min
            h′     = calc_current_health(m,h,zh,age,P)
            cv     = calc_cv(a′, h′, v0, Pzl[i_zl,:],Pzh[i_zh,:], P)
            u      = calc_utility(c, l, h′, P)
            p_surv = calc_prob_survival(age,h′,P)
            beq    = calc_warm_glow(a′,P)
            val    = u + β*p_surv*cv +β*(1-p_surv)*beq
        else
            val = -1e+38
        end
    else
        val = -1e+38
    end
    return val
end
"""
    function calc_cv(a′,h′, v0, Pzlvec, Pzhvec)

Description
"""
function calc_cv(a′::Float64, h′::Float64, v0::Array{Float64}, Pzlvec::Array{Float64},Pzhvec::Array{Float64},P::Params)
    n_zl = P.grid.n_zl
    n_zh = P.grid.n_zh
    v0_vals = eval_coefficients(v0, a′,h′,P)

    cv = 0.0
    for i_zl = 1:n_zl
        for i_zh = 1:n_zh
            cv = cv + Pzlvec[i_zl]*Pzhvec[i_zh]*v0_vals[i_zl,i_zh]
        end
    end
    return cv
end
"""
    function eval_coefficients(v0, a′, h′)

Description
"""
function eval_coefficients(v0::Array{Float64}, a′::Float64, h′::Float64,P::Params)
    a_grid = P.grid.a_grid
    h_grid = P.grid.h_grid
    n_zl = P.grid.n_zl
    n_zh = P.grid.n_zh
    a_max = P.grid.a_max
    a_min = P.grid.a_min
    h_max = P.grid.h_max
    h_min = P.grid.h_min

    v0_vals = zeros(n_zl,n_zh)

    for i_zl = 1:n_zl
        for i_zh = 1:n_zh
            interp = interpolate((a_grid, h_grid), v0[:,:,i_zl,i_zh], Gridded(Linear()))
            if a′> a_max
                if h′>h_max
                    v0_vals[i_zl] = interp(a_max,h_max)
                elseif h′<h_min
                    v0_vals[i_zl] = interp(a_max,h_min)
                else
                    v0_vals[i_zl] = interp(a_max,h′)
                end
            elseif a′<a_min
                if h′>h_max
                    v0_vals[i_zl] = interp(a_min,h_max)
                elseif h′<h_min
                    v0_vals[i_zl] = interp(a_min,h_min)
                else
                    v0_vals[i_zl] = interp(a_min,h′)
                end
            else
                if h′>h_max
                    v0_vals[i_zl] = interp(a′,h_max)
                elseif h′<h_min
                    v0_vals[i_zl] = interp(a′,h_min)
                else
                    v0_vals[i_zl] = interp(a′,h′)
                end
            end
        end
    end
    return v0_vals
end
