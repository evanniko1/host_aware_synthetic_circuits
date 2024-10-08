# import packages
using DifferentialEquations
using Sundials
using Plots

# import custom code
include("values.jl")    # values for initial conditions 
include("helper.jl")    # code for parameter sweeps and to organize ODE solutions for better handling 
include("host_aware_models.jl") # ODE definitions for host-aware synthetic circuits

# path to save figures
sv_path = "./figures/"

# time span
tspan = (0.0, 1e9)

# exposed relevant parameters
ns = 0.5
dprep = log(2)/4 # default value: log(2)/4
kappa_ini = 0.8  # 

# choose a host-aware model
model_def = HETER_ODE_model!

# specify initial conditions and parameter values
# for different host-aware models
if model_def == HETER_ODE_model!
    kurep = exp10(-2.6575) # default value: exp10(-2.6575), ku_rng[end]
    kbrep = exp10(-1.3335) # default value: exp10(-1.3335), kb_rng[end]
    wmaxrep = 150

    # organize in a vector
    het_p = [ns, dmrep, dprep, kappa_ini, wmaxrep, kbrep, kurep]

elseif model_def == REPR_ODE_model!
    kurep_1 = kurep_2 = kurep_3 = exp10(-2.6575)
    kbrep_1 = kbrep_2 = kbrep_3 = exp10(-1.3335)
    wmaxrep_1 = wmaxrep_2 = wmaxrep_3 = 250

    Kq_rep_1, Kq_rep_2, Kq_rep_3 = 100, 100, 100
    nq_rep_1, nq_rep_2, nq_rep_3 = 2, 2, 2

    # oganize in a vector
    het_p = [ns, dmrep, dprep, kappa_ini, wmaxrep_1, kbrep_1, kurep_1,
                                          wmaxrep_2, kbrep_2, kurep_2,
                                          wmaxrep_3, kbrep_3, kurep_3,
                                          Kq_rep_1, nq_rep_1,
                                          Kq_rep_2, nq_rep_2,
                                          Kq_rep_3, nq_rep_3]

    # need to extend the initial conditions vector
    u0_extra = [0, 0, 1000, 0, 0, 0]
    append!(u0, u0_extra)

elseif model_def == NOT_gate_ODE_model! 
    kurep_1 = kurep_2 = exp10(-2)
    kbrep_1 = kbrep_2 = exp10(-2)
    wmaxrep_1 = wmaxrep_2 = 250

    Kq_rep_1 = 250
    nq_rep_1 = 2

    # oganize in a vector
    het_p = [ns, dmrep, dprep, kappa_ini, wmaxrep_1, kbrep_1, kurep_1,
                                          wmaxrep_2, kbrep_2, kurep_2,
                                          Kq_rep_1, nq_rep_1]

    # need to extend the initial conditions vector
    u0_extra = [0, 0, 0]
    append!(u0, u0_extra)

elseif model_def == AND_gate_ODE_model!
    kurep_1 = kurep_2 = exp10(-2)
    kbrep_1 = kbrep_2 = exp10(-2)

    kbrep_3 = exp10(-1.5)
    kurep_3 = exp10(-2.5)

    wmaxrep_1 = wmaxrep_2 = 150
    wmaxrep_3 = 375

    Kq_rep_1, Kq_rep_2 = 200, 3000
    nq_rep_1, nq_rep_2 = 2.381, 1.835

    # oganize in a vector
    het_p = [ns, dmrep, dprep, kappa_ini, wmaxrep_1, kbrep_1, kurep_1,
                                          wmaxrep_2, kbrep_2, kurep_2,
                                          wmaxrep_3, kbrep_3, kurep_3,
                                          Kq_rep_1, nq_rep_1,
                                          Kq_rep_2, nq_rep_2]

    # need to extend the initial conditions vector
    u0_extra = [1, 1, 0, 1, 1, 0]
    append!(u0, u0_extra)

elseif model_def == NAND_gate_ODE_model!
    kurep_1 = kurep_2 = kurep_3 = kurep_4 = exp10(-2)
    kbrep_1 = kbrep_2 = kbrep_3 = kbrep_4 = exp10(-2)

    wmaxrep_1 = wmaxrep_2 = 250
    wmaxrep_3 = 375
    wmaxrep_4 = 250

    Kq_rep_1, Kq_rep_2, Kq_rep_3 = 200, 3000, 250
    nq_rep_1, nq_rep_2, nq_rep_3 = 2.381, 1.835, 2

    # oganize in a vector
    het_p = [ns, dmrep, dprep, kappa_ini, wmaxrep_1, kbrep_1, kurep_1,
                                          wmaxrep_2, kbrep_2, kurep_2,
                                          wmaxrep_3, kbrep_3, kurep_3,
                                          wmaxrep_4, kbrep_4, kurep_4,
                                          Kq_rep_1, nq_rep_1,
                                          Kq_rep_2, nq_rep_2,
                                          Kq_rep_3, nq_rep_3
                                          ]

    # need to extend the initial conditions vector
    u0_extra = [1, 100, 0, 100, 1, 0, 1, 100, 0]
    append!(u0, u0_extra)

else
    print("Model definition does not exist ... Define in *host_aware_models.jl* ...")
end

# assemble into single parameter vector
append!(p, het_p)

# create a problem dictionary
ode_problem_dict = create_problem_dict!(model_choice = model_def, 
                                        init_values = u0, 
                                        params_values = p, 
                                        tspan = tspan,
                                        ode_solver = Rodas4(autodiff=false),
                                        abstol = 1e-8,
                                        reltol = 1e-8,
                                        maxiters = 1e7,
                                        show_progress=true)

#--------------#
# solve a specified system and plot the change of species over time
solve_once = solve_ode_problem!(ode_problem_wrap = ode_problem_dict);
plot(solve_once) 

# reproduce the non-linear relationship from Cambray et al, 2018
phet_content, grate, ribosomal_content = trans_initiation!(ode_problem_dict = ode_problem_dict, range_size = 50, kini_lower = -0.65, kini_upper = 0);

plot(phet_content["protein_1"], grate, fmt = :svg)
plot(grate)

#--------------#
# single perturbation
prot_exp_1D, grate_1D, = perturb_one_param!(ode_problem_dict = ode_problem_dict, 
                                            param_index = 25, 
                                            range_bounds = (0, 4), 
                                            range_size = 50);
plot(prot_exp_1D["protein_1"]) # Fig 2.3B
plot(grate_1D) # Fig 2.3B (inset)
plot(prot_exp_1D["protein_1"], grate_1D)

#--------------#
# double perturbation
prot_exp_2D, grate_2D = perturb_two_params!(ode_problem_dict = ode_problem_dict, 
                                            param_index_inner = 24, 
                                            param_index_outer = 25, 
                                            range_bounds_inner = (0, 1), 
                                            range_bounds_outer = (0, 4), 
                                            range_size = 50);

plot(prot_exp_2D["protein_1"], 
     grate_2D, 
     palette = palette([:purple, :green], 50), 
     legend = false
     )

#--------------#
prot_exp_RBS, grate_RBS = perturb_param_w_RBS!(ode_problem_dict = ode_problem_dict, 
                                               param_index = 24, 
                                               range_bounds = (0, 1), 
                                               RBS_bounds = (-4, -2, 0), 
                                               range_size = 50)
plot(prot_exp_RBS["protein_1"], 
     grate_RBS,
     palette = palette([:purple, :green], 50), 
     legend = false
     )

#--------------#
prot_exp_RBS_only, grate_RBS_only = perturb_RBS!(ode_problem_dict = ode_problem_dict, RBS_bounds = (-4, -2, 0), range_size = 50)
plot(prot_exp_RBS_only["protein_1"], grate_RBS_only)

# plot heatmap
nn = heatmap(vector_to_matrix(prot_exp_2D["protein_3"]))

# save figure
save_figure(img_to_sv = nn, 
            model_def = model_def, 
            custom_suffix = "");