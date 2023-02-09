# import packages
using DifferentialEquations
using Plots

include("test_code.jl")


prob = ODEProblem(base_model_ODE!, base_model_ss_values, (0, 1e4), base_model_parameters)
sol = solve(prob, solver)

