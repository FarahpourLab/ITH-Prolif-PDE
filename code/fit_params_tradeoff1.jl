"""
This file fits the model to data from 11 tumor replicates.
"""

# read command line input
idx_str = ARGS[1]

path = "./"
file_data = "data/TV_control_long.csv.csv"
file_params = "results/fit_K-rho_hierarchical_bayesian.csv"
file_log_params = "results/logistic_fit_params_control.csv"

df_ctrl = CSV.read(path*file_data, DataFrame)
df_params = CSV.read(path*file_params, DataFrame)
df_log_params = CSV.read(path*file_log_params, DataFrame)

# get replicate ID from index
idx = parse(Int, idx_str)
rep_id = unique(df_ctrl.replicate_id)[idx]

function fit_trade_off_control(time_array, p)
    @parameters x, D_s, sigma_s, kd_s, s, kcl_s
    @variables c(..) y(..) N(..) integrand(..) x̄(..) 

    Dxx = Differential(x)^2
    Dtt = Differential(t)^2
    Dx = Differential(x)

    mu = df_params[666,2]
    K = df_log_params[idx, 2] # individual carrying capacity 

    V_init = df_ctrl[df_ctrl.replicate_id .== rep_id, :tumor_volume_smoothed][1]*0.001
    times_s = df_ctrl[(df_ctrl.replicate_id .== rep_id), :time_d].*mu

    rho_min = 0.0
    rho_max = 0.2 # 1/day

    # Scaled parameters
    mu_s = mu/mu
    V_init_s = V_init/K

    rho_min_s = rho_min/mu
    rho_max_s = rho_max/mu

    Ix = Integral(x in DomainSets.ClosedInterval(rho_min_s, rho_max_s))

    eq  = [Dt(c(t,x)) ~ D_s*Dxx(c(t, x)) + x*c(t,x)*(1 - N(t)) - kd_s * x^s * c(t,x)
       Dt(y(t,x)) ~ kd_s * x^s * c(t,x) - kcl_s * y(t,x)
       N(t) ~ Ix(c(t,x)) + Ix(y(t,x))
       integrand(t,x) ~ x * c(t, x)
       x̄(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)) + Ix(y(t,x)))]


    t0 = times_s[1]
    t1 = times_s[end]
    domains = [t ∈ (t0, t1),
    x ∈ (rho_min_s, rho_max_s)]

    bcs = [c(0,x) ~ (V_init_s/(sqrt(2*pi)*sigma_s))*exp(-((x-mu_s)^2)/(2*sigma_s^2)), 
    y(0,x) ~ 0.0, Dx(c(t, rho_min_s)) ~ 0.0, Dx(c(t, rho_max_s)) ~ 0, 
    Dx(y(t, rho_min_s)) ~ 0.0, Dx(y(t, rho_max_s)) ~ 0]

    @named pdesys = PDESystem(eq, bcs, domains, [t,x], 
        [c(t,x), y(t,x), N(t), integrand(t,x), x̄(t)], [D_s, sigma_s, kd_s, s, kcl_s], 
        defaults=Dict(D_s=>p[1], sigma_s=>p[2],kd_s => p[3], s => p[4], kcl_s => p[5]))
    dx=0.1
    discretization = MOLFiniteDifference([x => 60], t)

    prob = discretize(pdesys,discretization)
    sol = solve(prob, Tsit5(), saveat=times_s)

    sol_total = sol[N(t)]
    sol_mean = sol[x̄(t)]
    discrete_t = sol[t]

    vol_int = []
    mean_int = []

    for i in eachindex(discrete_t)
        if discrete_t[i] in times_s
            push!(vol_int, sol_total[i])
            push!(mean_int, sol_mean[i])
        end
    end
    n_pred = vcat(vol_int, mean_int)
    return n_pred
end

mu = df_params[666,2]
K = df_log_params[idx, 2] # individual carrying capacity 

times_s = df_ctrl[(df_ctrl.replicate_id .== rep_id), :time_d].*mu

p0 = [0.0058, 0.2, 0.05, 0.8, 3]
# D, sigma, kd_s, s, kcl_s

mean_vec = fill(mu/mu, length(times_s))
vols_vec = df_ctrl[(df_ctrl.replicate_id .== rep_id), :tumor_volume_smoothed].*0.001./K
data = vcat(vols_vec, mean_vec)

fit = curve_fit(fit_trade_off_control, times_s, data, p0; lower=[0.0029, 0.05, 0.0, 1.0, 0.001/mu])
params_fit = fit.param
@save path*"results/"*"params_fit_tradeoff1_$(rep_id).jld2" params_fit
