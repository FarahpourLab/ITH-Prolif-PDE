"""
## Visualization of the fitted models

With this code, the following plots were created:

- Figure 3 (a-c): 1, 4, 3 (for replicate 174, idx = 11)
- Figure S3 (a-j): 1
- Figure S4 (a-j): 3
- Figure S6 (a-c): 5, 8, 7
- Figure S7 (a-d): 5, 6, 7, 8
- Figure S10 (a-j): 5
- Figure S11 (a-j): 7

"""

path = "./"

file_data = "data/TV_control_long.csv.csv"
df_ctrl = CSV.read(path*file_data, DataFrame)

function plots_tradeoff1(duration="short")
  for idx in 1:11
    @load path*"results/sol_fit_tradeoff1_$(duration)_$(rep_id).jld2" discrete_x discrete_t solc sol_total sol_Nv sol_Nd solx̄ solx̄c
  
    times = df_ctrl[(df_ctrl.replicate_id .== rep_id), :time_d]
    vols_vec = df_ctrl[(df_ctrl.replicate_id .== rep_id), :tumor_volume_smoothed].*0.001
  
    ### 1) plot tumor volume with experimental data, color represents population-wide mean rho
    plt = plot(size=(600,400))
    plot!(plt, times, vols_vec, seriestype=:scatter, marker=:xcross, color=:black,
    label="Data", xlabel="Time [days]", ylabel="Tumor Volume [ml]", title="")
    plot!(plt, discrete_t, sol_total, label="", xlabel="Time [days]",  zcolor=solx̄, color=palette(:vikO), msc=:white,
    ylabel="Tumor Volume[ml]", title="", st=:scatter, colorbar_title=L"\bar{\rho}_{total}")
    savefig(plt,path*"results/TV_tradeoff1_$(duration)_$(rep_id).svg")

    if duration != "short"
      ### 2) plot total tumor volume with volume of doomed and viable subpopulations 
      plt = plot(size=(600,400))
      plot!(plt, discrete_t, sol_total, label=L"N_{total}", xlabel="Time [days]",  color=1,
      ylabel="Volume [ml]", title="")
      plot!(plt, discrete_t, sol_viable, label=L"N_v", xlabel="Time [days]",  color=2, 
      ylabel="Volume [ml]", title="")
      plot!(plt, discrete_t, sol_doomed, label=L"N_d", xlabel="Time [days]",  color=3, 
      ylabel="Volume [ml]", title="")
      savefig(plt,path*"results/TV_compartments_tradeoff1_$(duration)_$(rep_id).svg")
    end
  
    # calculate density at population-wide mean proliferation rate 
    ymean = zeros(length(discrete_t))
    for i in eachindex(discrete_t)
        itp = LinearInterpolation(discrete_x, solc[i, :], extrapolation_bc=Interpolations.Flat())
        ymean[i] = itp(solx̄[i])
    end
  
    # calculate density at mean proliferation rate of viable subpopulation
    ymean_cv = zeros(length(discrete_t))
    for i in eachindex(discrete_t)
        itp = LinearInterpolation(discrete_x, solc[i, :], extrapolation_bc=Interpolations.Flat())
        ymean_cv[i] = itp(solx̄c[i])
    end
  
    ### 3) plot density of C(t,rho) with or without mean proliferation rate (total=red, viable=orange)
    cols = get(ColorSchemes.imola, range(0.0, 1.0, length=length(sol[t])))
    plt = plot(size=(400,400))
    plot!(plt, xlims=(0,0.15))
    for i in eachindex(sol[t])
        Plots.plot!(plt, collect(sol[x]), sol[c(t,x)][i, :], label="", palette=cols, lw=1.5)
    end
    xlabel!(plt, L"\mathrm{\rho}"*" "*L"\mathrm{[day^{-1}}]")
    ylabel!(plt, L"C_v(t,\rho)")
    plot!(plt, solx̄c, ymean_cv, st=:scatter, color=:orangered, label=L"\bar{\rho}_{v}", 
    msc=:orangered, markersize=2)
    plot!(plt, solx̄, ymean, st=:scatter, color=:darkred, label=L"\bar{\rho}_{total}", 
    msc=:darkred, markersize=2)
    savefig(plt,path*"results/density_mean_tradeoff1_$(duration)_$(rep_id).svg")
  
    ### 4) plot mean proliferation rates over time (red=total, orange=viable)
    plt = plot(size=(400,300))
    plot!(plt, discrete_t, solx̄c, st=:scatter, msc=:white, xlabel="Time [days]", color=:orangered, 
    markersize=3, ylabel=L"\mathrm{\bar{\rho} \, [day^{-1}]}")
    plot!(plt, discrete_t, solx̄, st=:scatter, msc=:white, xlabel="Time [days]", color=:darkred,
    markersize=3, ylabel=L"\mathrm{\bar{\rho} \, [day^{-1}]}")
    savefig(plt,path*"results/mean_rho_tradeoff1_$(duration)_$(rep_id).svg")
    
  end
end

function plots_tradeoff1(duration="short")
  for idx in 1:11
    @load path*"results/sol_fit_tradeoff2_$(duration)_$(rep_id).jld2" discrete_x discrete_t solc sol_total sol_Nv sol_Nd solx̄ solx̄c
  
    times = df_ctrl[(df_ctrl.replicate_id .== rep_id), :time_d]
    vols_vec = df_ctrl[(df_ctrl.replicate_id .== rep_id), :tumor_volume_smoothed].*0.001
  
    ### 5) plot tumor volume with experimental data, color represents population-wide mean rho
    plt = plot(size=(600,400))
    plot!(plt, times, vols_vec, seriestype=:scatter, marker=:xcross, color=:black,
    label="Data", xlabel="Time [days]", ylabel="Tumor Volume [ml]", title="")
    plot!(plt, discrete_t, sol_total, label="", xlabel="Time [days]",  zcolor=solx̄, color=palette(:vikO), msc=:white,
    ylabel="Tumor Volume[ml]", title="", st=:scatter, colorbar_title=L"\bar{\rho}_{total}")
    savefig(plt,path*"results/TV_tradeoff2_$(duration)_$(rep_id).svg")
  
    if duration != "short"
      ### 6) plot total tumor volume with volume of doomed and viable subpopulations 
      plt = plot(size=(600,400))
      plot!(plt, discrete_t, sol_total, label=L"N_{total}", xlabel="Time [days]",  color=1,
      ylabel="Volume [ml]", title="")
      plot!(plt, discrete_t, sol_viable, label=L"N_v", xlabel="Time [days]",  color=2, 
      ylabel="Volume [ml]", title="")
      plot!(plt, discrete_t, sol_doomed, label=L"N_d", xlabel="Time [days]",  color=3, 
      ylabel="Volume [ml]", title="")
      savefig(plt,path*"results/TV_compartments_tradeoff2_$(duration)_$(rep_id).svg")
    end
  
    # calculate density at population-wide mean proliferation rate 
    ymean = zeros(length(discrete_t))
    for i in eachindex(discrete_t)
        itp = LinearInterpolation(discrete_x, solc[i, :], extrapolation_bc=Interpolations.Flat())
        ymean[i] = itp(solx̄[i])
    end
  
    # calculate density at mean proliferation rate of viable subpopulation
    ymean_cv = zeros(length(discrete_t))
    for i in eachindex(discrete_t)
        itp = LinearInterpolation(discrete_x, solc[i, :], extrapolation_bc=Interpolations.Flat())
        ymean_cv[i] = itp(solx̄c[i])
    end
  
    ### 7) plot density of C(t,rho) with or without mean proliferation rate (total=red, viable=orange)
    cols = get(ColorSchemes.imola, range(0.0, 1.0, length=length(sol[t])))
    plt = plot(size=(400,400))
    plot!(plt, xlims=(0,0.15))
    for i in eachindex(sol[t])
        Plots.plot!(plt, collect(sol[x]), sol[c(t,x)][i, :], label="", palette=cols, lw=1.5)
    end
    xlabel!(plt, L"\mathrm{\rho}"*" "*L"\mathrm{[day^{-1}}]")
    ylabel!(plt, L"C_v(t,\rho)")
    plot!(plt, solx̄c, ymean_cv, st=:scatter, color=:orangered, label=L"\bar{\rho}_{v}", 
    msc=:orangered, markersize=2)
    plot!(plt, solx̄, ymean, st=:scatter, color=:darkred, label=L"\bar{\rho}_{total}", 
    msc=:darkred, markersize=2)
    savefig(plt,path*"results/density_mean_tradeoff2_$(duration)_$(rep_id).svg")
  
    ### 8) plot mean proliferation rates over time (red=total, orange=viable)
    plt = plot(size=(400,300))
    plot!(plt, discrete_t, solx̄c, st=:scatter, msc=:white, xlabel="Time [days]", color=:orangered, 
    markersize=3, ylabel=L"\mathrm{\bar{\rho} \, [day^{-1}]}")
    plot!(plt, discrete_t, solx̄, st=:scatter, msc=:white, xlabel="Time [days]", color=:darkred,
    markersize=3, ylabel=L"\mathrm{\bar{\rho} \, [day^{-1}]}")
    savefig(plt,path*"results/mean_rho_tradeoff2_$(duration)_$(rep_id).svg")
  
  end
end
