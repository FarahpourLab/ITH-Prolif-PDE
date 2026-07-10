"""
## Smooth tumor trajectories

This code was used to create figure S1. 

"""

path = "./"

smoothed_ctrl = []
s = [30000, 15000, 10000, 30000, 10000, 10000, 10000, 15000, 10000, 10000, 30000]
plt = plot(layout=grid(3,4), size=(1200,600), left_margin=10mm, bottom_margin=5mm,
plot_title="Tumor Volumes in Control Group", plot_titlefontsize=14)
for idx in eachindex(unique(df_ctrl_replicates_comp.replicate_id))
    df_tmp = df_ctrl_replicates_comp[df_ctrl_replicates_comp.replicate_id .== unique(df_ctrl_replicates_comp.replicate_id)[idx], :]

    spl = Spline1D(df_tmp.time_d, df_tmp.tumor_volume_mm3, s=s[idx])
    xs = df_tmp.time_d
    ys = spl(xs)

    plot!(plt, subplot=idx, xs, ys, label="Smoothed", color=:blue)
    plot!(plt, subplot=idx, df_tmp.time_d, df_tmp.tumor_volume_mm3, markersize=3,
    color=:red, st=:scatter, msc=:white, title="$(string(unique(df_ctrl_replicates_comp.replicate_id)[idx]))",
    label="Data", legend_font_pointsize=8, 
    titlefontsize=12)
    if idx in [8,9,10,11]
        xlabel!(plt, subplot=idx, "Time [days]", xguidefontsize=11)
    end
    if idx in [1,5,9]
        ylabel!(plt, subplot=idx, L"\mathrm{[mm^3]}", yguidefontrotation=-90, yguidefontsize=11)
    end
    push!(smoothed_ctrl, ys)
end
for i=1:7; plot!(plt[i],xformatter=_->""); end
plt

df_ctrl_replicates_comp[!, "tumor_volume_smoothed"] = vcat(smoothed_ctrl...)
CSV.write(path*"data/TV_control_long.csv", df_ctrl_replicates_comp)
