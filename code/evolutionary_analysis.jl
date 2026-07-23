"""
## Evolutionary Analysis

This code was used to create the following figures:

- Figure 2 (a-b): 1, 3
- Figure S5 (a-b): 2, 4

"""

# --- PARAMETERS (Consistent across panels) ---
k_d = 0.5
s_val = 2
K = 1000.0
N₀ = 10.0

### 1) Life-history trade-off 

ρ_range = range(0.0, 0.6, length=100)
pA = plot(legend=:topleft, xlabel=L"\rho", ylabel="Mortality Rate", 
          title="(a) Mortality-Production Trade-off", grid=true)

for s in [0.5, 1.0, 2.0]
    mortality = k_d .* ρ_range.^s
    plot!(ρ_range, mortality, lw=2, label="s = $s")
end
savefig("./results/life_history_tradeoff.svg")

### 2) Fitness landscape at fixed N/K 

function g(ρ, N_K; s=s_val, k_d=k_d)
    ρ * (1 - N_K) - k_d * ρ^s
end

N_K_fixed = 0.7
g_slice = [g(ρ, N_K_fixed) for ρ in ρ_range]
ρ_star_val = ((1 - N_K_fixed)/(k_d * s_val))^(1/(s_val-1))
g_star = g(ρ_star_val, N_K_fixed)

pB = plot(ρ_range, g_slice, lw=2, legend=false, 
          xlabel=L"\rho", ylabel=L"g(\rho)", 
          title="(b) Fitness Landscape at N/K = $N_K_fixed",
          grid=true)

vline!([ρ_star_val], ls=:dash, color=:black)
scatter!([ρ_star_val], [g_star], color=:red)
savefig("./results/fitness_landscape.svg")

### 3) Tumor growth and ESS dynamics 

ρmax=0.6
function rho_star(N, K, kd, s; ρmax=0.4)
    if s > 1
        return clamp(((1 - N/K)/(kd*s))^(1/(s-1)), 0.0, ρmax)
    else
        # Fitness at both boundaries
        g0 = 0.0
        gmax = ρmax*(1 - N/K) - kd*ρmax^s
        return gmax > g0 ? ρmax : 0.0
    end
end

function tumor_growth!(du, u, p, t)
    N, = u
    kd, s, K, ρmax = p
    ρ = rho_star(N, K, kd, s; ρmax=ρmax)
    du[1] = ρ*N*(1 - N/K) - kd*ρ^s*N
end

s_values = [0.5, 1.0, 2.0]
colors = [:red, :green, :blue]

pC = plot(
    xlabel="Time (days)",
    ylabel="Tumor burden",
    title="(c) Tumor Growth and ESS Adaptation",
    legend=:bottomright,
    grid=true
)

for (s,c) in zip(s_values, colors)
    params = (k_d, s, K, ρmax)
    sol = solve(
        ODEProblem(tumor_growth!, [N₀], (0.0,80.0), params),
        Tsit5(),
        saveat=0.1
    )
    plot!(pC, sol.t, sol[1,:],
          color=c,
          lw=2,
          label="N(t), s=$s")
end

ax2 = twinx()

for (s,c) in zip(s_values, colors)
    params = (k_d, s, K, ρmax)
    sol = solve(
        ODEProblem(
            tumor_growth!,
            [N₀],
            (0.0, 80.0),
            params
        ),
        Tsit5(),
        saveat = 0.1
    )
    ρ = [rho_star(n, K, k_d, s; ρmax=ρmax) for n in sol[1,:]]
    plot!(ax2,
          sol.t,
          ρ,
          ls=:dash,
          color=c,
          label="ρ*(t), s=$s")
end
savefig("./results/TV_vs_rho_opt.svg")

### 4) Fitness landscape dynamics 

ρ0 = ρ_star[1]
ρf = ρ_star[end]

ρ_targets = [
    ρ0,
    ρ0 - (ρ0-ρf)/3,
    ρ0 - 2*(ρ0-ρf)/3
]

idx = [argmin(abs.(ρ_star .- ρ)) for ρ in ρ_targets]

time_points = t[idx]
N_K_values = N[idx] ./ K
ρ_values = ρ_star[idx]

time_points = round.(time_points, digits=0)

# plot 
s_panel = 2.0
params = (k_d, s_panel, K, ρmax)

sol = solve(
    ODEProblem(tumor_growth!, [N₀], (0.0, 80.0), params),
    Tsit5(),
    saveat=0.1
)

t = sol.t
N = sol[1,:]

ρ_star = [rho_star(n, K, k_d, s_panel; ρmax=ρmax) for n in N]

ρ_range = range(0.0, ρmax, length=200)
colors = [:green, :orange, :purple]

pD = plot(
    xlabel=L"\rho",
    ylabel=L"g(\rho)",
    title="(d) Fitness Landscape Transformation",
    grid=true,
    legend=:topright
)

for (i, t_i) in enumerate(time_points)

    idx = argmin(abs.(t .- t_i))

    N_K = N[idx] / K
    ρ_opt = ρ_star[idx]
    println(ρ_opt)

    g_vals = [ρ*(1 - N_K) - k_d*ρ^s_panel for ρ in ρ_range]

    label = "Day $(Int(t_i)) (N/K=$(round(N_K, digits=2)))"

    plot!(
        pD,
        ρ_range,
        g_vals,
        lw=2,
        color=colors[i],
        label=label
    )

    scatter!(
        pD,
        [ρ_opt],
        [ρ_opt*(1 - N_K) - k_d*ρ_opt^s_panel],
        color=colors[i],
        marker=:star5,
        markersize=8,
        label=""
    )
end
ylims!(pD, 0, 0.45)
savefig("./results/fitness_landscape_dynamics.svg")
