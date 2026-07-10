"""
## Run model with fitted parameters

The default setting for the simulation time corresponds to the experimental time 
frame (58 days). 
For long-term simulations, set t1 (l. 80) to a high value (e.g. 2000 days) and increase the
intervals between model saving times (l. 94, e.g. saveat=collect(t0:10.0:t1)) 
"""

function run_model_tradeoff1(t1=58)
    path = "./"
    file_data = "data/TV_control_long.csv.csv"
    file_params = "results/fit_K-rho_hierarchical_bayesian.csv"
    file_log_params = "results/logistic_fit_params_control.csv"
    
    df_ctrl = CSV.read(path*file_data, DataFrame)
    df_params = CSV.read(path*file_params, DataFrame)
    df_log_params = CSV.read(path*file_log_params, DataFrame)
    
    # fit was done with dimensionless model, now the parameters have to be transformed back 
    
    for idx in eachindex(unique(df_ctrl.replicate_id))
        rep_id = unique(df_ctrl.replicate_id)[idx]
        @load path*"results/"*"params_fit_tradeoff1_$(rep_id).jld2" params_fit
    
        @parameters x
        @variables c(..) d(..)  C_v(..) C_d(..) N(..) integrand(..) x̄(..) x̄c(..) 
    
        Dxx = Differential(x)^2
        Dtt = Differential(t)^2
        Dx = Differential(x)
    
        mu = df_params[666,2]
        K = df_log_params[idx, 2] # individual carrying capacity 
    
        kd_s = params_fit[3]
        s = params_fit[4]
        kcl_s = params_fit[5]
    
        kcl = kcl_s * mu
        kd = kd_s * mu^(1-s)
    
        V_init = df_ctrl[df_ctrl.replicate_id .== rep_id, :tumor_volume_smoothed][1]*0.001
        times = df_ctrl[(df_ctrl.replicate_id .== rep_id), :time_d]
    
        rho_min = 0.0
        rho_max = 0.2 # 1/day
    
        sigma_s = params_fit[2]
        sigma = sigma_s * mu
    
        D_s = params_fit[1]
        D = D_s * mu^3
    
        Ix = Integral(x in DomainSets.ClosedInterval(rho_min, rho_max))
    
        eq  = [Dt(c(t,x)) ~ D*Dxx(c(t, x)) + x*c(t,x)*(1 - N(t)/K) - kd * x^s * c(t,x)
           Dt(d(t,x)) ~ kd * x^s * c(t,x) - kcl * d(t,x)
           C_v(t) ~ Ix(c(t,x))
           C_d(t) ~ Ix(d(t,x))
           N(t) ~ Ix(c(t,x)) + Ix(d(t,x))
           integrand(t,x) ~ x * c(t, x)
           x̄(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)) + Ix(d(t,x)))
           x̄c(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)))]
    
        t0 = times[1]
        t1 = t1
        interval = 1.0
        duration = "short"
        if t1 >= 1000
            interval = 10.0
            duration = "long"
        end
        domains = [t ∈ (t0, t1),
        x ∈ (rho_min, rho_max)]
    
        bcs = [c(0,x) ~ (V_init/(sqrt(2*pi)*sigma))*exp(-((x-mu)^2)/(2*sigma^2)), 
        d(0,x) ~ 0.0, Dx(c(t, rho_min)) ~ 0.0, Dx(c(t, rho_max)) ~ 0, 
        Dx(d(t, rho_min)) ~ 0.0, Dx(d(t, rho_max)) ~ 0]
    
        @named pdesys = PDESystem(eq, bcs, domains, [t,x], 
            [c(t,x), d(t,x), C_v(t), C_d(t), N(t), integrand(t,x), x̄(t), x̄c(t)])
        dx=0.1
        discretization = MOLFiniteDifference([x => 60], t)
    
        prob = discretize(pdesys,discretization)
        sol = solve(prob, Tsit5(), saveat=collect(t0:interval:t1))
    
        vols_vec = df_ctrl[(df_ctrl.replicate_id .== rep_id), :tumor_volume_smoothed].*0.001
    
        # sol = solutions[idx]
        discrete_x = sol[x]
        discrete_t = sol[t]
        solc = sol[c(t,x)]
        sold = sol[d(t,x)]
        sol_Nv = sol[C_v(t)]
        sol_Nd = sol[C_d(t)]
        sol_total = sol[N(t)]
        solx̄ = sol[x̄(t)]
        solx̄c = sol[x̄c(t)]
    
        @save path*"results/sol_fit_tradeoff1_$(duration)_$(rep_id).jld2" discrete_x discrete_t solc sol_total sol_Nv sol_Nd solx̄ solx̄c
    end
end
