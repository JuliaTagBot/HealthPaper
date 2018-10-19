function MC(vfi::VFIResults, w::Float64, r::Float64, ss::Float64, τw::Float64, T::Float64, P::Params)
    Random.seed!(1)

    g_c = vfi.g_c
    g_l = vfi.g_l
    g_a′ = vfi.g_a′
    g_m = vfi.g_m
    g_h′ = vfi.g_h′

    n_age = P.D.n_age
    R_age = P.D.R_age
    N  = P.D.NN
    n_pop = P.D.n_pop
    α = P.F.α
    A = P.F.A
    τc = P.G.τc
    τa = P.G.τa
    c̲ = P.G.c̲
    g_ratio = P.G.g_ratio
    Pzl1_cdf = P.grid.Pzl1_CDF
    Pzl_cdf = P.grid.Pzl_CDF
    Pzh1_cdf = P.grid.Pzh1_CDF
    Pzh_cdf = P.grid.Pzh_CDF
    a_grid = P.grid.a_grid
    h_grid = P.grid.h_grid
    h_min = P.grid.h_min
    h_max = P.grid.h_max
    a_min = P.grid.a_min
    a_max = P.grid.a_max
    zl_grid = P.grid.zl
    zh_grid = P.grid.zh
    zh1_cdf = P.grid.Pzh1_CDF
    sim_a=zeros(N,n_age)
    sim_m = zeros(N,n_age)
    sim_h=zeros(N,n_age)
    sim_zl=zeros(Int,N)
    sim_zh=zeros(Int,N)
    sim_a′=zeros(N,n_age)
    sim_h′=zeros(N,n_age)
    sim_zln=zeros(Int,N)
    sim_c=zeros(N,n_age)
    sim_y=zeros(N,n_age)
    sim_wt=zeros(N,n_age)
    sim_l=zeros(N,n_age)
    sim_le=zeros(N,n_age)
    sim_e=zeros(N,n_age)
    sim_lab=zeros(N,n_age)
    sim_ones=ones(N,n_age)
    sim_psurv=zeros(N,n_age)
    sim_acc_beq=zeros(N,n_age)
    sim_death_shock=zeros(N,n_age)
    B_vec=zeros(n_age)
    N_alive = zeros(Int,n_age)
    psihat=zeros(n_age)

    K_vec=zeros(n_age)
    C_vec = zeros(n_age)
    L_vec= zeros(n_age)
    M_vec = zeros(n_age)
    Le_vec = zeros(n_age)
    WT_vec = zeros(n_age)
    acc_beq_vec = zeros(n_age)

    #allocate ( K_vec(n_age), L_vec(n_age), acc_beq_vec(n_age), C_vec(n_age), F_vec(n_age), M_vec(n_age) )
    #allocate ( mu(n_age), psihat(n_age), muhat(n_age), muhatn(n_age) )

    # set weights for each generation, to normalize for population growth
    # this doesn't account for the survival probabilities
    mu=ones(n_age)

    for i = 2:n_age
        mu[i] = mu[i-1] / (1.0 + n_pop)
    end
    # normalize mu
    mu = mu ./ sum(mu)

    # initialize values for simulated agents
    #sim_a  = zeros(N)     # age = 1 agents start with no assets
    sim_zl = simulate_initial_markov_shocks(N, Pzl1_cdf)
    sim_h[:,1] = simulate_initial_health_shocks(N::Int64,P::Params)
    sim_zh = simulate_initial_markov_shocks(N, Pzh1_cdf)
    sim_alive = trues(N,n_age)

      for age = 1:n_age

        # compute number agents still alive in simulation
        N_alive[age] = length( sim_alive[:,age][sim_alive[:,age].==true] )

        if (N_alive[age] == 0)
            println("warning: no agents still alive by age ", age)
            # is there a better way of doing the Monte Carlo simulation, so that the
            # number of simulated agents doesn't drop with age?
        end

        # compute social security benefit for this age
        if age<R_age
            b_val=0.0
        else
            b_val=ss
        end


        # average benefit payments to this generation
        B_vec[age] = b_val

        # fraction of agents sill alive at this age
        # psihat(j) = fraction of agents still alive at age j
        # psihat(1) = 1 by construction
        psihat[age] = N_alive[age] / N

        # now, update distribution for next age group:

        for i = 1:N
            if sim_alive[i,age]
                if sim_h[i,age]< h_min
                    sim_h[i,age]=h_min
                elseif sim_h[i,age]> h_max
                    sim_h[i,age]=h_max
                end
                if sim_a[i,age]<a_min
                    sim_a[i,age]=a_min
                elseif sim_a[i,age]>a_max
                    sim_a[i,age]=a_max
                end

                a′ = interpolate((a_grid, h_grid), g_a′[:,:,sim_zl[i],sim_zh[i],age], Gridded(Linear()))(sim_a[i,age],sim_h[i,age])
                l  =interpolate((a_grid, h_grid), g_l[:,:,sim_zl[i],sim_zh[i],age], Gridded(Linear()))(sim_a[i,age],sim_h[i,age])
                m  =interpolate((a_grid, h_grid), g_m[:,:,sim_zl[i],sim_zh[i],age], Gridded(Linear()))(sim_a[i,age],sim_h[i,age])
                sim_h′[i,age] = calc_current_health(m,sim_h[i,age],zh_grid[sim_zh[i]],age,P)

                y0  = calc_non_labour_income(τa,r,sim_a[i,age],T,b_val)
                el = calc_labour_productivity(zl_grid[sim_zl[i]],sim_h′[i,age],age,P)
                lab_inc =calc_labour_income(zl_grid[sim_zl[i]],sim_h′[i,age], age, w, l, τw, P)
                tot_inc = y0+lab_inc
                wt = calc_transfers(tot_inc,m,P)
                y = calc_net_inc(zl_grid[sim_zl[i]],sim_h′[i,age], age, w, y0, l, m, a′,τw,P)
                c = (y - a′-m)/(1+τc)
                # get other policy choices
                sim_a′[i,age] = a′
                sim_l[i,age] = l
                sim_m[i,age] = m
                sim_y[i,age] = y
                sim_wt[i,age] = wt
                sim_c[i,age] = c
                sim_e[i,age] = calc_labour_productivity(zl_grid[sim_zl[i]],sim_h′[i,age],age,P)

                # labor productivity
                zl_val = zl_grid[sim_zl[i]]
                sim_lab[i,age] = l * el

                # probability of survival
                sim_psurv[i,age] = calc_prob_survival(age,sim_h′[i,age],P)
                # expected accidental bequests
                sim_acc_beq[i,age] = (1.0 - sim_psurv[i,age]) * sim_a′[i,age] * (1.0+r*(1.0-τa))

            end
        end


        if (N_alive[age] > 0)
            # calculate average capital for this age group:
            K_vec[age] = (sum( sim_a′[:,age][sim_alive[:,age].==true] ) / (1.0+n_pop)) / N_alive[age]
            # average consumption for this age group
            C_vec[age] = sum( sim_c[:,age][sim_alive[:,age].==true]) / N_alive[age]
            # average labor
            M_vec[age] = sum( sim_m[:,age][sim_alive[:,age].==true]) / N_alive[age]
            L_vec[age] = sum( sim_l[:,age][sim_alive[:,age].==true]) / N_alive[age]
            Le_vec[age] = sum( sim_lab[:,age][sim_alive[:,age].==true]) / N_alive[age]
            # accidental bequests
            acc_beq_vec[age] = sum(sim_acc_beq[:,age][sim_alive[:,age].==true]) / N_alive[age]
            WT_vec[age] = sum( sim_wt[:,age][sim_alive[:,age].==true]) / N_alive[age]
        else
            K_vec[age] = 0.0
            C_vec[age] = 0.0
            M_vec[age] = 0.0
            L_vec[age] = 0.0
            Le_vec[age] = 0.0
            acc_beq_vec[age] = 0.0
        end

        # now simulate next period's shocks
        sim_zln = simulate_markov_shocks(N, sim_zl, Pzl_cdf)
        sim_zhn = simulate_markov_shocks(N,sim_zh,Pzh_cdf)
        #if mod(age,10)==0
        #    println( "age = ", age)
        #end


        # update initial values
        if age<n_age
            sim_a[:,age+1]  = copy(sim_a′[:,age])
            sim_h[:,age+1] = copy(sim_h′[:,age])
        end
        sim_zl = copy(sim_zln)
        sim_zh = copy(sim_zhn)

        # determine whether each individual survives
        sim_death_shock= rand(N)
        for i = 1:N
            if sim_alive[i,age] && sim_death_shock[i] > sim_psurv[i,age] && age<n_age
                sim_alive[i,age+1] = false
            end
        end

    end

    # compute the endogenous distribution of agents across ages
    muhat = mu.*psihat
    # normalize the distribution -- use this in the calculations below
    # to ensure that everything is in per-capita terms.  However, keep
    # in mind, the total size of the economy's population is endogenous
    muhatn = muhat / sum(muhat)  # normalize

    # now calculate totals
    # K,L are defined in per-capita terms
    K_val = dot( muhatn, K_vec )    # capital
    L_val = dot( muhatn, L_vec)     # labor
    C_val = dot( muhatn, C_vec)
    M_val = dot(muhatn, M_vec)
    Le_val = dot(muhatn, Le_vec)
    # aggregate output (per-capita)
    Y_val = A * (K_val^α) * (Le_val^(1.0-α))

    # calculate accidental bequests (per-capita)
    T_val = dot( muhatn, acc_beq_vec ) / (1.0 + n_pop)
    WT_val = dot(muhatn, WT_vec)
    # government spending (per-capita)
    G_val = g_ratio * Y_val

    # total benefit payments (per-capita)
    B_val = dot( muhatn, B_vec )

    # total consumption (per-capita)
    C_val = dot( muhatn, C_vec )
    # total medical expenditures (per-capita)

    # update tauw
    tauw_val = (G_val + B_val +WT_val - τc*C_val - τa*r*K_val) / (w*Le_val)

    # save results
    #save("results.jld","sim_a", sim_a,"sim_h",sim_h,"sim_a′",sim_a′,"sim_c",sim_c)
    return MonteCarloResults(K_val, L_val, C_val, M_val, Le_val, T_val, WT_val, tauw_val, muhatn, psihat, K_vec, L_vec,Le_vec, C_vec,M_vec, WT_vec, acc_beq_vec,sim_a′,sim_c,sim_m,sim_h′,sim_l,sim_lab,sim_wt)
end


function calc_markov_shock_from_uniform_rv(u::Float64, Pz_cdf::Array{Float64})


    done = false
    i_zn = 0

    # u is uniformly distributed over [0,1)
    # convert it into an index for productivity
    # using the cumulative probability transition matrix
    while done==false
        i_zn = i_zn + 1;

        # stop if u < Pz_cdf(i_zn)
        # u has to be less than 1
        # Pz_cdf(nz) == 1, so eventually this will be true
        done = (u < Pz_cdf[i_zn]);
    end

    return i_zn
end

    #**************************************************************************
    # function simulate_markov_shocks
    #**************************************************************************
function simulate_markov_shocks(N::Int64, sim_z::Array{Int64}, Pz_cdf::Array{Float64})

    sim_zl=zeros(Int,N)
    # allocate memory for u:

    u=rand(N)
    # generate random variables, which is
    # uniformly distributed over [0,1)


    # back out markov shock from u for each z:
    for i = 1:N
        sim_zl[i] = calc_markov_shock_from_uniform_rv(u[i], Pz_cdf[sim_z[i],:]);
    end
    return sim_zl
end

    #**************************************************************************
    # function simulate_initial_markov_shocks
    #**************************************************************************
function simulate_initial_markov_shocks(N::Int64, Pz_initial_cdf::Array{Float64})


    sim_z=zeros(Int,N)
    u=rand(N)
    # generate random variables, which is
    # uniformly distributed over [0,1)


    # back out markov shock from u for each z:
    for i = 1:N
        sim_z[i] = calc_markov_shock_from_uniform_rv(u[i], Pz_initial_cdf);
    end
    return sim_z
end

    #**************************************************************************
    # function simulate_initial_health_shocks
    #**************************************************************************
    #**************************************************************************
function simulate_initial_health_shocks(N::Int64,P::Params)

    μₕ = P.grid.μₕ
    σₕ = P.grid.σₕ
    h_min = P.grid.h_min
    h_max = P.grid.h_max

    # allocate memory for u:

    sim_h0=zeros(N)
    # generate random variables, which is
    # uniformly distributed over [0,1)
    u = rand(N)

    # back out markov shock from u for each z:
    for i = 1:N
        # convert u to normal
        #sim_h0(i) = mu_h0 + sigma_h0 * dinvnr(u(i), 1.0d0-u(i))
        sim_h0[i] = μₕ + σₕ*norminvcdf( u[i] )

        # force sim_h0 to be on grid
        sim_h0[i] = max( sim_h0[i], h_min )
        sim_h0[i] = min( sim_h0[i], h_max )
    end
    return sim_h0
end
