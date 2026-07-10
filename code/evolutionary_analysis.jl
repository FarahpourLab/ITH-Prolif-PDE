"""
Read in raw data xlsx file and write long dataframe
"""

path = "./data/"
fname_ctrl = "TVs_control_raw.xlsx"

xf_ctrl = XLSX.readxlsx(path*fname_ctrl)

data_ctrl = xf_ctrl[XLSX.sheetnames(xf_ctrl)[1]]

ids_ctrl = data_ctrl["A"][3:13, :]

times_ctrl = [0,2,4,7,10,12,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,48,50,52,54,56,58]

measurements_ctrl = data_ctrl["B:DI"][3:13, :]
cols_ctrl = collect(3:4:111) # select only column with tumor volumes
measurements_ctrl = measurements_ctrl[:, cols_ctrl]
measurements_ctrl = permutedims(measurements_ctrl)

df_ctrl = DataFrame(measurements_ctrl, :auto)
rename!(df_ctrl, string.(vec(ids_ctrl)))
df_ctrl[!, "time_d"] = times_ctrl

df_ctrl_replicates = stack(df_ctrl, Not("time_d"), variable_name = "replicate_id", value_name="tumor_volume_mm3")
df_ctrl_replicates = df_ctrl_replicates[completecases(df_ctrl_replicates), :]

stats_ctrl = combine(groupby(df_ctrl_replicates, :time_d),
                :tumor_volume_mm3 => mean => :vmean,
                :tumor_volume_mm3 => std  => :vstd)

df_ctrl_replicates_comp = leftjoin(df_ctrl_replicates, stats_ctrl, on = :time_d)

CSV.write(path*"TV_control_long.csv", df_ctrl_replicates_comp)
