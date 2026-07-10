# Model for _Eco-Evolutionary Dynamics of Proliferation Heterogeneity: A Phenotype-Structured Model for Tumor Growth and Treatment Response_

This repository contains the code used in: **Eco-Evolutionary Dynamics of Proliferation Heterogeneity: A Phenotype-Structured Model for Tumor Growth and Treatment Response** ([link to bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.13.711687v2))

# Overview

ITH-Prolif-PDE is a computational framework for simulating eco-evolutionary dynamics in tumors with heterogeneous proliferation rates. This phenotype-structured PDE model captures three fundamental biological principles: (1) continuous variation in proliferation phenotypes, (2) global resource competition, and (3) proliferation-mortality trade-offs. The framework enables quantitative prediction of tumor growth and therapy response through evolutionary dynamics.

The code reproduces the results shown in the main text and supplementary material of the manuscript.

# Repository Structure 

```bash
ith-prolif-pde
├── code
│   ├── compare_ode_vs_pde.jl
│   ├── data_preprocessing.jl
│   ├── evolutionary_analysis.jl
│   ├── fit_params_tradeoff1.jl
│   ├── fit_params_tradeoff2.jl
│   ├── main.jl # main script 
│   ├── plots_fitted_model.jl
│   ├── run_fitted_model_tradeoff1.jl
│   ├── run_fitted_model_tradeoff2.jl
│   ├── save_fitted_params.jl
│   ├── smooth_raw_TVs.jl
│   ├── treatments_tradeoff1.jl
│   └── treatments_tradeoff2.jl
├── results
│   ├── fit_K-rho_hierarchical_bayesian.csv
│   ├── fits_tradeoff1.csv
│   ├── fits_tradeoff2.csv
│   ├── logistic_fit_params_control.csv
│   ├── params_fit_tradeoff1_$(rep_id).jld2 # for each replicate 
│   └── params_fit_tradeoff2_$(rep_id).jld2 # for each replicate
├── .gitignore
├── Manifest.toml
├── Project.toml
└── README.md
```

# Reproduce the results

1. Install `Julia v1.12.0` from the [official website](https://julialang.org/downloads/)
2. Clone or download this repository
3. Set up the julia environment

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
4. Run the code

```julia
main.jl
```

**Note:** Some parts of the code (i.e., the PDE fits in `fit_params_tradeoff1.jl` and `fit_params_tradeoff2.jl`) are computationally intensive and were run on a server. We provide the results of these computations in the `results/` folder (i.e., `fits_tradeoff1.csv` and `fits_tradeoff2.csv`). You do not need to rerun this step to reproduce the results in the manuscript. 
