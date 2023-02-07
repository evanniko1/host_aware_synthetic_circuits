# import packages
using DifferentialEquations
using Plots

# import custom stuff
include("values.jl")
include("helper.jl")

# init plotly
plotly();

# time span
tspan = (0.0, 1e5)
##########
# exposed relevant parameters
ns = 0.5
kurep = exp10(-2.6575)
kbrep = exp10(-1.3335)
wmaxrep = 150
kappa_ini = 0.1
# oganize in a vector
het_p = [ns, wmaxrep, kbrep, kurep, dmrep, dprep, kappa_ini]

# assemble into single parameter vector
append!(p, het_p)

ode_problem_dict = create_problem_dict!(init_values = u0, params_values = p, tspan = tspan)


##########
# one run of the model
#solve_once = solve_ode_problem!(model_def = ODE_model!, ode_problem_wrap = ode_problem_dict)

# specify and solve ODE problem
#plot(solve_once) # time trajectories for all species; the plot is interactive

##########
# reproduce the non-linear relationship from Cambray et al, 2018
phet, grate = trans_initiation!(ode_problem_dict = ode_problem_dict);
plot(phet, grate);


# TEST ZONE
# single perturbation -- WORKS
ll, kk, _ = perturb_one_param!(ode_problem_dict = ode_problem_dict, param_index = 27, range_bounds = (-0.65, 0), range_size = 10);

# double perturbation -- WORKS
ll, kk, _ = perturb_two_params!(ode_problem_dict = ode_problem_dict, param_index_inner = 27, param_index_outer = 26, range_bounds_inner = (-0.65, 0), range_bounds_outer = (0, 2), range_size = 10);
