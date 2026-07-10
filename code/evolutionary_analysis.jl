"""
## Evolutionary Analysis

This code was used to create the following figures:

- Figure 2 (a-b): 1, 3
- Figure S5 (a-b): 2, 4

"""

### 1) Life-history trade-off 

ρ_range = range(0.0, 0.2, length=100)
pA = plot(legend=:topleft, xlabel=L"\rho", ylabel="Death Rate", 
          title="(a) Life-history trade-off", grid=true)

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

display(pC)
savefig("./results/TV_vs_rho_opt.svg")

### 4) Fitness landscape dynamics 

ρ_range = range(0.0, 0.4, length=100)
time_points = [5.0, 25, 50.0]  # Early, mid, late
colors = [:green :orange :purple]
pD = plot(xlabel=L"\rho", ylabel=L"g(\rho)", 
          title="(d) Fitness Landscape Dynamics",
          grid=true, legend=:topright, ylims=(0,0.18))

for (i, t_i) in enumerate(time_points)
    idx = argmin(abs.(t .- t_i))
    N_K_i = N[idx]/K
    ρ_star_i = ρ_star[idx]
    g_vals = [g(ρ, N_K_i) for ρ in ρ_range]
    
    label_str = "Day $(Int(t_i))"*L"(\mathrm{\frac{N}{K}}"*"=$(round(N_K_i, digits=2)))"
    plot!(ρ_range, g_vals, lw=2, color=colors[i], label=label_str)
    scatter!([ρ_star_i], [g(ρ_star_i, N_K_i)], color=colors[i], marker=:star, markersize=6)
    println(ρ_star_i)
end
savefig("./results/fitness_landscape_dynamics.svg")
