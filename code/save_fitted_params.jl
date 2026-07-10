"""
This file collects and saves the fitted parameters for both trade-off terms. 
The fitting is run on the dimensionless version of the model. The parameters are 
converted to their respective units before saving. 
"""

path = "./"

file_data = "data/TV_control_long.csv"
df_ctrl = CSV.read(path*file_data, DataFrame)

file_params = "results/fit_K-rho_hierarchical_bayesian.csv"
df_params = CSV.read(path*file_params, DataFrame)

### Collect and export fitted parameters for monotonically increasing trade-off term

D_list = Float64[]
sigma_list = Float64[]
kd_list = Float64[]
s_list = Float64[]
kcl_list = Float64[]
for idx in 1:11
    rep_id = unique(df_ctrl.replicate_id)[idx]
    @load path*"results/"*"params_fit_tradeoff1_$(rep_id).jld2" params_fit
    
    mu = df_params[666,2]

    D_s = params_lsqfit[1]
    sigma_s = params_lsqfit[2]
    kd_s = params_lsqfit[3]
    s = params_lsqfit[4]
    kcl_s = params_lsqfit[5]

    kcl = kcl_s * mu
    kd = kd_s * mu^(1-s)
    D = D_s * mu^3
    sigma = sigma_s * mu

    push!(D_list, D)
    push!(sigma_list, sigma)
    push!(kd_list, kd)
    push!(s_list, s)
    push!(kcl_list, kcl)
end

df_fit_params = DataFrame(replicate_id = unique(df_ctrl.replicate_id),
kd = kd_list, s = s_list, kcl = kcl_list, D = D_list, sigma = sigma_list)

CSV.write(path*"results/fits_tradeoff1.csv", df_fit_params)

### Collect and export fitted parameters for u-shaped trade-off term

D_list = Float64[]
sigma_list = Float64[]
kd_list = Float64[]
xopt_list = Float64[]
kcl_list = Float64[]
for idx in 1:11
    rep_id = unique(df_ctrl.replicate_id)[idx]
    @load path*"results/"*"params_fit_tradeoff1_$(rep_id).jld2" params_fit

    mu = df_params[666,2]

    D_s = params_lsqfit[1]
    sigma_s = params_lsqfit[2]
    kd_s = params_lsqfit[3]
    x_opt_s = params_lsqfit[4]
    kcl_s = params_lsqfit[5]

    kcl = kcl_s * mu
    kd = kd_s * mu
    D = D_s * mu^3
    sigma = sigma_s * mu
    x_opt = x_opt_s * mu

    push!(D_list, D)
    push!(sigma_list, sigma)
    push!(kd_list, kd)
    push!(xopt_list, x_opt)
    push!(kcl_list, kcl)
end

df_fit_params = DataFrame(replicate_id = unique(df_ctrl.replicate_id),
kd = kd_list, x_opt = xopt_list, kcl = kcl_list, D = D_list, sigma = sigma_list)

CSV.write(path*"results/fits_tradeoff2.csv", df_fit_params)
