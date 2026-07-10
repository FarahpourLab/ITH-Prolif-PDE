"""
## Simulation and visualization of example model

With this code, the following plots were created:

- Figure 1 (a-d): 1,2,3,4
- Figure 1 a (inset): 5
- Figure 1 d (inset): 6

"""

path = "./"

@parameters x
@variables c(..) y(..) N(..) C_v(..) C_d(..) integrand(..) x̄(..) x̄c(..)

Dxx = Differential(x)^2
Dtt = Differential(t)^2
Dx = Differential(x)

# Simplified parameter values 
mu = 0.05516 # 1/day
K = 1 # ml 
kcl = 0.001 # 1/day
s = 1.4
kd = 0.07
sigma = 0.007
D = 5e-7

# Initial conditions and boundary values 
V_init = 0.1 # ml 
times = [0, 2, 4, 7, 10, 12, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 48, 50, 52, 54, 56, 58, 60]

rho_min = 0.0
rho_max = 0.2 

Ix = Integral(x in DomainSets.ClosedInterval(rho_min, rho_max))

eq  = [Dt(c(t,x)) ~ D*Dxx(c(t, x)) + x*c(t,x)*(1 - N(t)/K) - kd * x^s * c(t,x)
    Dt(y(t,x)) ~ kd * x^s * c(t,x) - kcl * y(t,x)
    N(t) ~ Ix(c(t,x)) + Ix(y(t,x))
    C_v(t) ~ Ix(c(t,x))
    C_d(t) ~ Ix(y(t,x))
    integrand(t,x) ~ x * c(t, x)
    x̄(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)) + Ix(y(t,x)))
    x̄c(t) ~ Ix(integrand(t,x))/(Ix(c(t,x)))]

# Run model for long time (2000 days)

t0 = times[1]
t1 = 2000
domains = [t ∈ (t0, t1),
x ∈ (rho_min, rho_max)]

bcs = [c(0,x) ~ (V_init/(sqrt(2*pi)*sigma))*exp(-((x-mu)^2)/(2*sigma^2)), 
y(0,x) ~ 0.0, Dx(c(t, rho_min)) ~ 0.0, Dx(c(t, rho_max)) ~ 0, 
Dx(y(t, rho_min)) ~ 0.0, Dx(y(t, rho_max)) ~ 0]

@named pdesys = PDESystem(eq, bcs, domains, [t,x], 
    [c(t,x), y(t,x), N(t), C_v(t), C_d(t), integrand(t,x), x̄(t), x̄c(t)])
dx=0.1
discretization = MOLFiniteDifference([x => 120], t)

prob = discretize(pdesys,discretization)
sol = solve(prob, Tsit5(), saveat=t0:10.0:t1)

discrete_x = sol[x]
discrete_t = sol[t]
solc = sol[c(t,x)]
sol_total = sol[N(t)]
sol_viable = sol[C_v(t)]
sol_doomed = sol[C_d(t)]
solx̄ = sol[x̄(t)]
solx̄c = sol[x̄c(t)]

### 1) Plot total tumor volume over time
plt = plot(size=(600,400))
plot!(plt, discrete_t, sol_total, label="", xlabel="Time [days]",  zcolor=solx̄, color=palette(:vikO), msc=:white,
ylabel="Volume [ml]", title="", st=:scatter, colorbar_title=L"\bar{\rho}_{total}")
savefig(plt,path*"results/Ntotal_colored_long.svg")

### 2) Plot total tumor volume and volumes of viable and doomed subpopulations 
plt = plot(size=(600,400))
plot!(plt, discrete_t, sol_total, label=L"N_{total}", xlabel="Time [days]",  color=1,
ylabel="Volume [ml]", title="")
plot!(plt, discrete_t, sol_viable, label=L"N_v", xlabel="Time [days]",  color=2, 
ylabel="Volume [ml]", title="")
plot!(plt, discrete_t, sol_doomed, label=L"N_d", xlabel="Time [days]",  color=3, 
ylabel="Volume [ml]", title="")
savefig(plt,path*"results/Ntotal_Nv_Nd_long.svg")

### 3) Plot distribution of viable subpopulation over time 
cols = get(ColorSchemes.imola, range(0.0, 1.0, length=length(sol[t])))
plt = plot(size=(600,400))
plot!(plt, xlims=(0,0.15))
for i in eachindex(sol[t])
    Plots.plot!(plt, collect(sol[x]), sol[c(t,x)][i, :], label="", palette=cols, lw=1.5)
end
xlabel!(plt, L"\mathrm{\rho}"*" "*L"\mathrm{[day^{-1}}]")
ylabel!(plt, L"C_v(t,\rho)")
savefig(plt,path*"results/Cv_long.svg")

### 4) Plot mean proliferation rate rho_total and rho_v over time
plt = plot(size=(600,400))
plot!(plt, discrete_t, solx̄c, st=:scatter, msc=:white, xlabel="Time [days]", color=:orangered, 
markersize=3, ylabel=L"\mathrm{\bar{\rho} \, [day^{-1}]}")
plot!(plt, discrete_t, solx̄, st=:scatter, msc=:white, xlabel="Time [days]", color=:darkred,
markersize=3, ylabel=L"\mathrm{\bar{\rho} \, [day^{-1}]}")
savefig(plt,path*"results/rhototal_rhov_long.svg")

# Run model for short time (60 days)

t0 = times[1]
t1 = times[end]
domains = [t ∈ (t0, t1),
x ∈ (rho_min, rho_max)]

bcs = [c(0,x) ~ (V_init/(sqrt(2*pi)*sigma))*exp(-((x-mu)^2)/(2*sigma^2)), 
y(0,x) ~ 0.0, Dx(c(t, rho_min)) ~ 0.0, Dx(c(t, rho_max)) ~ 0, 
Dx(y(t, rho_min)) ~ 0.0, Dx(y(t, rho_max)) ~ 0]

@named pdesys = PDESystem(eq, bcs, domains, [t,x], 
    [c(t,x), y(t,x), N(t), C_v(t), C_d(t), integrand(t,x), x̄(t), x̄c(t)])
dx=0.1
discretization = MOLFiniteDifference([x => 120], t)

prob = discretize(pdesys,discretization)
sol = solve(prob, Tsit5(), saveat=times)

discrete_x = sol[x]
discrete_t = sol[t]
solc = sol[c(t,x)]
sol_total = sol[N(t)]
sol_viable = sol[C_v(t)]
sol_doomed = sol[C_d(t)]
solx̄ = sol[x̄(t)]
solx̄c = sol[x̄c(t)]

### 5) Plot total tumor volume over time (for inset)
plt = plot(size=(600,400))
plot!(plt, discrete_t, sol_total, label="", xlabel="Time [days]",  zcolor=solx̄, color=palette(:vikO), msc=:white,
ylabel="Volume [ml]", title="", st=:scatter, colorbar_title=L"\bar{\rho}_{total}")
savefig(plt,path*"results/Ntotal_colored_short.svg")

### 6) Plot mean proliferation rate rho_total and rho_v over time (for inset)
plt = plot(size=(600,400))
plot!(plt, discrete_t, solx̄c, st=:scatter, msc=:white, xlabel="Time [days]", color=:orangered, 
markersize=3, ylabel=L"\mathrm{\bar{\rho} \, [day^{-1}]}")
plot!(plt, discrete_t, solx̄, st=:scatter, msc=:white, xlabel="Time [days]", color=:darkred,
markersize=3, ylabel=L"\mathrm{\bar{\rho} \, [day^{-1}]}")
savefig(plt,path*"results/rhototal_rhov_short.svg")
