"""
## Treatment explorations 

This code was used to create the following figures:

- Figure S8 (a-f): 1, 2, 3, 4, 5, 6
- Figure S9 (a-b): 7, 8

"""

path = "./"
file_data = "data/TV_control_long.csv.csv"
file_params = "results/fit_K-rho_hierarchical_bayesian.csv"
file_log_params = "results/logistic_fit_params_control.csv"

df_ctrl = CSV.read(path*file_data, DataFrame)
df_params = CSV.read(path*file_params, DataFrame)
df_log_params = CSV.read(path*file_log_params, DataFrame)

delta_cont(t) = t >= 5 && t <= 15 ? 1.0 : 0.0
@register_symbolic delta_cont(t)

## Example for one tumor replicate (replicate 133)

rep_id = 133
@load path*"results/"*"params_fit_tradeoff2_$(rep_id).jld2" params_fit 

@parameters x
@variables c(..) y(..) N(..) C_v(..) C_d(..) integrand(..) x̄(..) x̄c(..) f(..)

Dxx = Differential(x)^2
Dtt = Differential(t)^2
Dx = Differential(x)

mu = df_params[666,2]
K = df_log_params[6, 2]

kd_s = params_fit[3]
x_opt_s = params_fit[4]
kcl_s = params_fit[5]

kcl = kcl_s * mu
kd = kd_s * mu

sigma_s = params_fit[2]
sigma = sigma_s * mu
x_opt = x_opt_s * mu

D_s = params_fit[1]
D = D_s * mu^3

V_init = df_ctrl[df_ctrl.replicate_id .== rep_id, :tumor_volume_smoothed][1]*0.001
times = df_ctrl[(df_ctrl.replicate_id .== rep_id), :time_d]

rho_min = 0.0
rho_max = 0.2 # 1/day

t_dose = 10.0
kt = 0.43
n = 1.5
sig = 0.2*mu

Ix = Integral(x in DomainSets.ClosedInterval(rho_min, rho_max))

### Untreated 

eq  = [Dt(c(t,x)) ~ D*Dxx(c(t, x)) + x*c(t,x)*(1 - N(t)/K) - kd * ((x-x_opt)^2/((x-x_opt)^2+x*(rho_max-x))) * c(t,x) - delta_cont(t) * f(t,x) * c(t,x) 
    Dt(y(t,x)) ~ kd * ((x-x_opt)^2/((x-x_opt)^2+x*(rho_max-x))) * c(t,x) - kcl * y(t,x) + delta_cont(t) * f(t,x) * c(t,x)
    N(t) ~ Ix(c(t,x)) + Ix(y(t,x))
    C_v(t) ~ Ix(c(t,x))
    C_d(t) ~ Ix(y(t,x))
    integrand(t,x) ~ x * c(t, x)
    x̄(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)) + Ix(y(t,x)))
    x̄c(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)))
    f(t,x) ~ 0]

t0 = times[1]
t1 = times[end]
domains = [t ∈ (t0, t1),
x ∈ (rho_min, rho_max)]

bcs = [c(0,x) ~ (V_init/(sqrt(2*pi)*sigma))*exp(-((x-mu)^2)/(2*sigma^2)), 
y(0,x) ~ 0.0, Dx(c(t, rho_min)) ~ 0.0, Dx(c(t, rho_max)) ~ 0, 
Dx(y(t, rho_min)) ~ 0.0, Dx(y(t, rho_max)) ~ 0]

@named pdesys = PDESystem(eq, bcs, domains, [t,x], 
    [c(t,x), y(t,x), N(t), C_v(t), C_d(t), integrand(t,x), x̄(t), x̄c(t), f(t,x)])
dx=0.1
discretization = MOLFiniteDifference([x => 240], t)

prob = discretize(pdesys,discretization)
sol = solve(prob, Tsit5(), saveat=t0:2.0:t1)

vols_vec = df_ctrl[(df_ctrl.replicate_id .== rep_id), :tumor_volume_smoothed].*0.001

discrete_x = sol[x]
discrete_t = sol[t]
solc_untreated = sol[c(t,x)]
soln_untreated = sol[N(t)]
sol_viable_untreated = sol[C_v(t)]
sol_doomed_untreated = sol[C_d(t)]
solx̄_untreated = sol[x̄(t)]
solx̄c_untreated = sol[x̄c(t)]

### Uniform treatment 

eq  = [Dt(c(t,x)) ~ D*Dxx(c(t, x)) + x*c(t,x)*(1 - N(t)/K) - kd * ((x-x_opt)^2/((x-x_opt)^2+x*(rho_max-x))) * c(t,x) - delta_cont(t) * f(t,x) * c(t,x) 
    Dt(y(t,x)) ~ kd * ((x-x_opt)^2/((x-x_opt)^2+x*(rho_max-x))) * c(t,x) - kcl * y(t,x) + delta_cont(t) * f(t,x) * c(t,x)
    N(t) ~ Ix(c(t,x)) + Ix(y(t,x))
    C_v(t) ~ Ix(c(t,x))
    C_d(t) ~ Ix(y(t,x))
    integrand(t,x) ~ x * c(t, x)
    x̄(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)) + Ix(y(t,x)))
    x̄c(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)))
    f(t,x) ~ 0.1]

t0 = times[1]
t1 = times[end]
domains = [t ∈ (t0, t1),
x ∈ (rho_min, rho_max)]

bcs = [c(0,x) ~ (V_init/(sqrt(2*pi)*sigma))*exp(-((x-mu)^2)/(2*sigma^2)), 
y(0,x) ~ 0.0, Dx(c(t, rho_min)) ~ 0.0, Dx(c(t, rho_max)) ~ 0, 
Dx(y(t, rho_min)) ~ 0.0, Dx(y(t, rho_max)) ~ 0]

@named pdesys = PDESystem(eq, bcs, domains, [t,x], 
    [c(t,x), y(t,x), N(t), C_v(t), C_d(t), integrand(t,x), x̄(t), x̄c(t), f(t,x)])
dx=0.1
discretization = MOLFiniteDifference([x => 240], t)

prob = discretize(pdesys,discretization)
sol = solve(prob, Tsit5(), saveat=t0:2.0:t1)

solc_uniform = sol[c(t,x)]
soln_uniform = sol[N(t)]
sol_viable_uniform = sol[C_v(t)]
sol_doomed_uniform = sol[C_d(t)]
solx̄_uniform = sol[x̄(t)]
solx̄c_uniform = sol[x̄c(t)]

### Low-rho targeting treatment 

eq  = [Dt(c(t,x)) ~ D*Dxx(c(t, x)) + x*c(t,x)*(1 - N(t)/K) - kd * ((x-x_opt)^2/((x-x_opt)^2+x*(rho_max-x))) * c(t,x) - delta_cont(t) * f(t,x) * c(t,x) 
    Dt(y(t,x)) ~ kd * ((x-x_opt)^2/((x-x_opt)^2+x*(rho_max-x))) * c(t,x) - kcl * y(t,x) + delta_cont(t) * f(t,x) * c(t,x)
    N(t) ~ Ix(c(t,x)) + Ix(y(t,x))
    C_v(t) ~ Ix(c(t,x))
    C_d(t) ~ Ix(y(t,x))
    integrand(t,x) ~ x * c(t, x)
    x̄(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)) + Ix(y(t,x)))
    x̄c(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)))
    f(t,x) ~ 0.26 * exp(-20*x)]

t0 = times[1]
t1 = times[end]
domains = [t ∈ (t0, t1),
x ∈ (rho_min, rho_max)]

bcs = [c(0,x) ~ (V_init/(sqrt(2*pi)*sigma))*exp(-((x-mu)^2)/(2*sigma^2)), 
y(0,x) ~ 0.0, Dx(c(t, rho_min)) ~ 0.0, Dx(c(t, rho_max)) ~ 0, 
Dx(y(t, rho_min)) ~ 0.0, Dx(y(t, rho_max)) ~ 0]

@named pdesys = PDESystem(eq, bcs, domains, [t,x], 
    [c(t,x), y(t,x), N(t), C_v(t), C_d(t), integrand(t,x), x̄(t), x̄c(t), f(t,x)])
dx=0.1
discretization = MOLFiniteDifference([x => 240], t)

prob = discretize(pdesys,discretization)
sol = solve(prob, Tsit5(), saveat=t0:2.0:t1)

solc_low = sol[c(t,x)]
soln_low = sol[N(t)]
sol_viable_low = sol[C_v(t)]
sol_doomed_low = sol[C_d(t)]
solx̄_low = sol[x̄(t)]
solx̄c_low = sol[x̄c(t)]

### Mid-rho targeting 

eq  = [Dt(c(t,x)) ~ D*Dxx(c(t, x)) + x*c(t,x)*(1 - N(t)/K) - kd * ((x-x_opt)^2/((x-x_opt)^2+x*(rho_max-x))) * c(t,x) - delta_cont(t) * f(t,x) * c(t,x) 
    Dt(y(t,x)) ~ kd * ((x-x_opt)^2/((x-x_opt)^2+x*(rho_max-x))) * c(t,x) - kcl * y(t,x) + delta_cont(t) * f(t,x) * c(t,x)
    N(t) ~ Ix(c(t,x)) + Ix(y(t,x))
    C_v(t) ~ Ix(c(t,x))
    C_d(t) ~ Ix(y(t,x))
    integrand(t,x) ~ x * c(t, x)
    x̄(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)) + Ix(y(t,x)))
    x̄c(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)))
    f(t,x) ~ kt * exp(-(x-mu)^2/(2*sig^2))]

t0 = times[1]
t1 = times[end]
domains = [t ∈ (t0, t1),
x ∈ (rho_min, rho_max)]

bcs = [c(0,x) ~ (V_init/(sqrt(2*pi)*sigma))*exp(-((x-mu)^2)/(2*sigma^2)), 
y(0,x) ~ 0.0, Dx(c(t, rho_min)) ~ 0.0, Dx(c(t, rho_max)) ~ 0, 
Dx(y(t, rho_min)) ~ 0.0, Dx(y(t, rho_max)) ~ 0]

@named pdesys = PDESystem(eq, bcs, domains, [t,x], 
    [c(t,x), y(t,x), N(t), C_v(t), C_d(t), integrand(t,x), x̄(t), x̄c(t), f(t,x)])
dx=0.1
discretization = MOLFiniteDifference([x => 240], t)

prob = discretize(pdesys,discretization)
sol = solve(prob, Tsit5(), saveat=t0:2.0:t1)

solc_mid = sol[c(t,x)]
soln_mid = sol[N(t)]
sol_viable_mid = sol[C_v(t)]
sol_doomed_mid = sol[C_d(t)]
solx̄_mid = sol[x̄(t)]
solx̄c_mid = sol[x̄c(t)]

### High-rho targeting 

eq  = [Dt(c(t,x)) ~ D*Dxx(c(t, x)) + x*c(t,x)*(1 - N(t)/K) - kd * ((x-x_opt)^2/((x-x_opt)^2+x*(rho_max-x))) * c(t,x) - delta_cont(t) * f(t,x) * c(t,x) 
    Dt(y(t,x)) ~ kd * ((x-x_opt)^2/((x-x_opt)^2+x*(rho_max-x))) * c(t,x) - kcl * y(t,x) + delta_cont(t) * f(t,x) * c(t,x)
    N(t) ~ Ix(c(t,x)) + Ix(y(t,x))
    C_v(t) ~ Ix(c(t,x))
    C_d(t) ~ Ix(y(t,x))
    integrand(t,x) ~ x * c(t, x)
    x̄(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)) + Ix(y(t,x)))
    x̄c(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)))
    f(t,x) ~ 0.01 * exp(30*x)]

t0 = times[1]
t1 = times[end]
domains = [t ∈ (t0, t1),
x ∈ (rho_min, rho_max)]

bcs = [c(0,x) ~ (V_init/(sqrt(2*pi)*sigma))*exp(-((x-mu)^2)/(2*sigma^2)), 
y(0,x) ~ 0.0, Dx(c(t, rho_min)) ~ 0.0, Dx(c(t, rho_max)) ~ 0, 
Dx(y(t, rho_min)) ~ 0.0, Dx(y(t, rho_max)) ~ 0]

@named pdesys = PDESystem(eq, bcs, domains, [t,x], 
    [c(t,x), y(t,x), N(t), C_v(t), C_d(t), integrand(t,x), x̄(t), x̄c(t), f(t,x)])
dx=0.1
discretization = MOLFiniteDifference([x => 240], t)

prob = discretize(pdesys,discretization)
sol = solve(prob, Tsit5(), saveat=t0:2.0:t1)

solc_high = sol[c(t,x)]
soln_high = sol[N(t)]
sol_viable_high = sol[C_v(t)]
sol_doomed_high = sol[C_d(t)]
solx̄_high = sol[x̄(t)]
solx̄c_high = sol[x̄c(t)]

## Plots 

### 1) Treatment profiles 

uniform_death(x) = 0.1
decreasing_death(x) = 0.26 * exp(-20*x)
increasing_death(x) = 0.01 * exp(30*x)
gaussian_death(x) = 0.43 * exp(-(x-mu)^2/(2*sig^2))

plot(uniform_death, 0, 0.12, label="uniform", ylims=(0,0.55), width=2)
plot!(decreasing_death, 0, 0.12, label="low "*L"\rho", width=2)
plot!(increasing_death, 0, 0.12, label="high "*L"\rho", width=2)
plot!(gaussian_death, 0, 0.12, label="mid "*L"\rho", width=2)
plot!(title="Treatment terms", xlabel="Proliferation rate "*L"\rho\,[\mathrm{day^{-1}}]", 
ylabel=L"[\mathrm{day^{-1}}]", size=(500,400), legend_title="Targeting:")
savefig(path*"results/treatment_profiles.svg")

### 2) tumor volumes after treatment 

plt = plot(500,400)
plot!(plt, discrete_t, soln_uniform, label="", xlabel="Time [days]", width=2, color=1,
ylabel="Tumor Volume [ml]")
plot!(plt, discrete_t, soln_low, label=" ", width=2, color=2)
plot!(plt, discrete_t, soln_high, label=" ", width=2, color=3)
plot!(plt, discrete_t, soln_mid, label=" ", width=2, color=4)
plot!(plt, discrete_t, soln_untreated, label="", color=:black, width=2)
savefig(plt,path*"results/TV_treatments_tradeoff2.svg")

### 3) Density C_v(t,rho) after uniform treatment 

pink = colorant"pink"
hawaii = get(ColorSchemes.hawaii, range(0.0, 1.0, length=length(discrete_t)-3))
cols = vcat(fill(pink,3), hawaii)

plt = plot(size=(500,400), xlims=(0,0.15))
for i in eachindex(discrete_t)
    Plots.plot!(plt, discrete_x, solc_uniform[i, :], label="", palette=cols, lw=1.5)
end
xlabel!(plt, L"\mathrm{\rho}"*" "*L"\mathrm{[day^{-1}}]")
ylabel!(plt, L"C_v(t,"*L"\rho"*")")
savefig(plt,path*"results/density_uniform_tradeoff2.svg")

### 4) Density C_v(t,rho) after low-rho targeting treatment 

plt = plot(size=(500,400), xlims=(0,0.15))
for i in eachindex(discrete_t)
    Plots.plot!(plt, discrete_x, solc_low[i, :], label="", palette=cols, lw=1.5)
end
xlabel!(plt, L"\mathrm{\rho}"*" "*L"\mathrm{[day^{-1}}]")
ylabel!(plt, L"C_v(t,"*L"\rho"*")")
savefig(plt,path*"results/density_low_tradeoff2.svg")

### 5) Density C_v(t,rho) after mid-rho targeting treatment 

plt = plot(size=(500,400), xlims=(0,0.15))
for i in eachindex(discrete_t)
    Plots.plot!(plt, discrete_x, solc_mid[i, :], label="", palette=cols, lw=1.5)
end
xlabel!(plt, L"\mathrm{\rho}"*" "*L"\mathrm{[day^{-1}}]")
ylabel!(plt, L"C_v(t,"*L"\rho"*")")
savefig(plt,path*"results/density_mid_tradeoff2.svg")

### 6) Density C_v(t,rho) after high-rho targeting treatment 

plt = plot(size=(500,400), xlims=(0,0.15))
for i in eachindex(discrete_t)
    Plots.plot!(plt, discrete_x, solc_high[i, :], label="", palette=cols, lw=1.5)
end
xlabel!(plt, L"\mathrm{\rho}"*" "*L"\mathrm{[day^{-1}}]")
ylabel!(plt, L"C_v(t,"*L"\rho"*")")
savefig(plt,path*"results/density_high_tradeoff2.svg")

### 7) Population-wide mean proliferation rate after treatment

plt = plot()
plot!(plt, discrete_t, solx̄_uniform, label="uniform", xlabel="Time [days]", st=:scatter,
ylabel=L"\mathrm{\bar{\rho}_{total} \, [day^{-1}]}", title="Mean proliferation rate "*L"\bar{\rho}_{total}(t)", msc=:white)
plot!(plt, discrete_t, solx̄_decreasing, label="low "*L"\rho", st=:scatter, msc=:white)
plot!(plt, discrete_t, solx̄_increasing, label="high "*L"\rho", st=:scatter, msc=:white)
plot!(plt, discrete_t, solx̄_gaussian, label="mid "*L"\rho", st=:scatter, msc=:white)
savefig(plt,path*"results/mean_total_trx_tradeoff2.svg")

### 8) Mean proliferation rate in viable subpopulation after treatment 

plt = plot()
plot!(plt, discrete_t, solx̄c_uniform, label="uniform", xlabel="Time [days]", st=:scatter,
ylabel=L"\mathrm{\bar{\rho}_{v} \, [day^{-1}]}", title="Mean proliferation rate "*L"\bar{\rho}_{v}(t)", msc=:white)
plot!(plt, discrete_t, solx̄c_low, label="low "*L"\rho", st=:scatter, msc=:white)
plot!(plt, discrete_t, solx̄c_high, label="high "*L"\rho", st=:scatter, msc=:white)
plot!(plt, discrete_t, solx̄c_mid, label="mid "*L"\rho", st=:scatter, msc=:white)
savefig(plt,path*"results/mean_viable_trx_tradeoff2.svg")
