## Import packages 

include("import_packages.jl")

## Include scripts

include("run_fitted_model_tradeoff1.jl")
include("run_fitted_model_tradeoff2.jl")
include("plots_fitted_model.jl")

## Run model with example parameters

include("run_representative_model.jl")

## Run fitted model for monotonically increasing trade-off function 

run_model_tradeoff1(t1=58) # experimental time frame 
run_model_tradeoff1(t1=2000) # long-term behavior

## Run fitted model for u-shaped trade-off function 

run_model_tradeoff2(t1=58) # experimental time frame 
run_model_tradeoff2(t1=2000) # long-term behavior

## Plot PDE model solutions 

plots_tradeoff1(duration="short") # monotonically increasing trade-off, experimental time frame 
plots_tradeoff1(duration="long") # monotonically increasing trade-off, long-term  
plots_tradeoff2(duration="short") # u-shaped trade-off, experimental time frame 
plots_tradeoff2(duration="long") # u-shaped trade-off, long-term  

## Treatment exploration for monotonically increasing trade-off term 

include("treatments_tradeoff1.jl")

## Treatment exploration for u-shaped trade-off term 

include("treatments_tradeoff2.jl")

## Compare ODE model to PDE model 

include("compare_ode_vs_pde.jl")
