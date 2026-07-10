"""
## Compare ODE vs. PDE model solutions 

This code was used to create figure S2a. 

"""

### Create and solve ODE model for replicate 130 

path = "./"

df_ctrl = CSV.read(path*"data/tv_control.csv", DataFrame)
df_params = CSV.read(path*"results/fit_K-rho_hierarchical_bayesian.csv", DataFrame)
df_log_params = CSV.read(path*"results/logistic_fit_params_control.csv", DataFrame)

rep_id = 130

volumes = df_ctrl[df_ctrl.replicate_id .== rep_id, :tumor_volume_smoothed].*0.001
times = df_ctrl[df_ctrl.replicate_id .== rep_id, :time_d]

V0 = volumes[1]
rho = df_params[666,2] # 0.05561 day^-1, fitted using hierarchical Bayesian model 
K = df_params[42, 2] # 0.3353, fitted individually for replicate 130 using Bayesian approach

function ode_model!(du, u, p, t)
    du[1] = mu * (1 - (u[1]/K)) * u[1] 
end

prob_ode = ODEProblem(ode_model!, [V0], (times[1], times[end]))
sol_ode = solve(prob_ode, Tsit5(), saveat=1)

### Create and solve PDE model for replicate 130 with different initial standard deviations

@parameters x 
@variables c(..) N(..) integrand(..) x̄(..)
D = 1e-6 

Dxx = Differential(x)^2
Dtt = Differential(t)^2
Dx = Differential(x)

stds = [0.1, 0.15, 0.2, 0.25, 0.3]
pde_sols = []
for sigma in stds
    rep_id = 130
    times = df_ctrl[df_ctrl.replicate_id .== rep_id, :time_d]
    volumes = df_ctrl[df_ctrl.replicate_id .== rep_id, :tumor_volume_smoothed].*0.001
    V_0 = volumes[1]

    rho = df_params[666,2] 
    mu = df_params[42, 2] 
    std = mu * sigma

    Ix = Integral(x in DomainSets.ClosedInterval(0.0, 0.2))

    eq  = [Dt(c(t,x)) ~ D*Dxx(c(t,x)) + x*c(t,x)*(1-(N(t)/K))
        N(t) ~ Ix(c(t,x))
        integrand(t,x) ~ x * c(t, x)
        x̄(t) ~ Ix(integrand(t,x))/Ix(c(t,x))]

    domains = [t ∈ (times[1],times[end]),
    x ∈ (0.0,0.2)]

    bcs = [c(0,x) ~ 1/(sqrt(2*pi*std^2))*exp(-((x-mu)^2)/(2*std^2))*V_0,
    Dx(c(t, 0.0)) ~ 0.0, Dx(c(t, 0.2)) ~ 0]

    @named pdesys = PDESystem(eq, bcs, domains, [t,x], 
            [c(t,x), N(t), integrand(t,x), x̄(t)])
    dx=0.1
    discretization = MOLFiniteDifference([x => 120], t)

    prob = discretize(pdesys,discretization)
    sol = solve(prob, Tsit5(), saveat=1)

    push!(pde_sols, sol)
end

rep_id = 130
times = df_ctrl[df_ctrl.replicate_id .== rep_id, :time_d]
volume_data = df_ctrl[df_ctrl.replicate_id .== rep_id, :tumor_volume_smoothed].*0.001

sol_01 = pde_sols[1]
sol_015 = pde_sols[2]
sol_02 = pde_sols[3]
sol_025 = pde_sols[4]
sol_03 = pde_sols[5]

plt = plot(times, volume_data, seriestype=:scatter, marker=:xcross,
label="Data", color=:black, size=(600,600))
plot!(plt, sol_01[t], sol_01[N(t)], label=L"0.1 \cdot \bar{\rho}_0", color=RGBA(1, 128/255, 229/255, 1))
plot!(plt, sol_015[t], sol_015[N(t)], label=L"0.15 \cdot \bar{\rho}_0", color=RGBA(1, 42/255, 212/255, 1))
plot!(plt, sol_02[t], sol_02[N(t)], label=L"0.2 \cdot \bar{\rho}_0", color=RGBA(212/255, 0, 170/255, 1))
plot!(plt, sol_025[t], sol_025[N(t)], label=L"0.25 \cdot \bar{\rho}_0", color=RGBA(128/255, 0, 102/255, 1))
plot!(plt, sol_03[t], sol_03[N(t)], label=L"0.3 \cdot \bar{\rho}_0", color=RGBA(43/255, 0, 34/255, 1))
plot!(plt, sol_log.t, first.(sol_log.u), label="ODE Model", color=RGBA(0, 0, 1, 1))
xlabel!(plt, "Time [days]")
ylabel!(plt, L"\mathrm{[ml]}", yguidefontrotation=-90, left_margin=5mm)
savefig(path*"results/ode_vs_pde.svg")
display(plt)
