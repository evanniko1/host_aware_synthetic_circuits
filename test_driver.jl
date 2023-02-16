#### TO-DO: 
# 1. ModelingToolkit.jl implementation -- DONE
# 2. solve and check steady state solution -- DONE
# 3. explore optimization with Jacobian -- DONE
#    1. several bits to check from documentation (?)
# 4. explore connected systems 
#    1. simplify growth model into blocks
#    2. if 4.1 == true, then examine host-circuit extensions
#    -- problems :: how to update shared resources ODEs ??
# 5. explore structural simplifications to build -if applicable- the leanest numerical representation of the system
# 6. is it possible to build a DAE version of the model?


# import packages
using ModelingToolkit
using DifferentialEquations
using Plots
plotly()

include("test_code.jl")

###################
# ModelingToolkit implementation
base_eqs, base_model_ss_values, base_model_parameters, ttrate_b = base_model();
ha_eqs = heter_model()
# base model with 14 ODEs
@named base_model_sys = ODESystem(base_eqs);

#tspan = (0.0,10000.0)
#prob = ODEProblem(base_model_sys, base_model_ss_values, tspan, base_model_parameters;jac=true,sparse=true)
#sol = solve(prob, Rodas4())
#plot(sol)

@named heterologous_sys = ODESystem(ha_eqs)
# try a structural simplification on the host-aware model
ha_eqn_simple = structural_simplify(heterologous_sys);

equations(ha_eqn_simple)


