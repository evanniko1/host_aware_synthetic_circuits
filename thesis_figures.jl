# This file contains a list of configurations to run the simulations that correspond to all the figures in the thesis.
# For each figure it is **imperative** to re-execute lines 13-23. This is a fast step.
# Following this, you need to run the code listed under each respective figure section.


# import packages
using DifferentialEquations
using Sundials
using Plots
using LaTeXStrings

# import custom stuff
include("values.jl")
include("helper.jl")
include("host_aware_models.jl")

# path to save figures
sv_path = "./figures/"

tspan = (0.0, 1e9) #the time span to simulate our system
ns = 0.5           #efficiency of nutrient quality
dprep = log(2)/4   #active degradation rate for heterologous proteins
kappa_ini = 0.8

#-------------------#
#     Chapter 2     #
#-------------------# 

# Figure 2.3
# define a host-aware model
model_def = HETER_ODE_model!

kbrep = kurep = exp10(-2)
wmaxrep = 0 #no induction strength for the heterologous construct is equivalent to simulating the wild strain

# assemble in a vector
het_p = [ns, dmrep, dprep, kappa_ini, wmaxrep, kbrep, kurep]
# append to the rest of the model parameters
append!(p, het_p)

# create problem dict
ode_problem_dict = create_problem_dict!(model_choice = model_def, 
                                        init_values = u0, 
                                        params_values = p, 
                                        tspan = tspan,
                                        ode_solver = Rodas4(autodiff=false),
                                        abstol = 1e-8,
                                        reltol = 1e-8,
                                        maxiters = 1e7,
                                        show_progress=true)

# Figure 2.3A :: CALL FUNCTION TO DRAW PIES
solve_once = solve_ode_problem!(ode_problem_wrap = ode_problem_dict);

# Figure 2.3B
# parameter sweep for induction strength
prot_exp_1D, grate_1D, = perturb_one_param!(ode_problem_dict = ode_problem_dict, 
                                            param_index = 25, 
                                            range_bounds = (0, 4), 
                                            range_size = 50);

# plot protein expresssion vs induction strength
prot_exp_het = plot(prot_exp_1D["parameter_values"],
                    prot_exp_1D["protein_1"],
                    label = false,
                    xlabel = "gene induction (mRNAs/min)",
                    ylabel = "expression (# of molecules)",
                    xaxis = :log10)

save_figure(img_to_sv = prot_exp_het, 
            model_def = model_def, 
            custom_suffix = "Figure2_3B"); 

# plot growth rate vs induction strength
grate_het = plot(prot_exp_1D["parameter_values"],
                 grate_1D,
                 label = false,
                 xlabel = "gene induction (mRNAs/min)",
                 ylabel = "growth rate (min\$^{-1}\$)",
                 xaxis = :log10)

save_figure(img_to_sv = grate_het,
            model_def = model_def,
            custom_suffix = "Figure2_3B(inset)")
#---------------------------------#

# Figure 2.5
# define a host-aware NOT gate
model_def = NOT_gate_ODE_model!

kurep_1 = kurep_2 = exp10(-2)
kbrep_1 = kbrep_2 = exp10(-2)
wmaxrep_1 = wmaxrep_2 = 250

Kq_rep_1 = 250
nq_rep_1 = 2

# oganize in a vector
het_p = [ns, dmrep, dprep, kappa_ini, wmaxrep_1, kbrep_1, kurep_1,
                                      wmaxrep_2, kbrep_2, kurep_2,
                                      Kq_rep_1, nq_rep_1]

# assemble into single parameter vector
append!(p, het_p)

# need to extend the initial conditions vector
u0_extra = [0, 0, 0]
append!(u0, u0_extra)

# create problem dict
ode_problem_dict = create_problem_dict!(model_choice = model_def, 
                                        init_values = u0, 
                                        params_values = p, 
                                        tspan = tspan,
                                        ode_solver = Rodas4(autodiff=false),
                                        abstol = 1e-8,
                                        reltol = 1e-8,
                                        maxiters = 1e7,
                                        show_progress=true)
                                  
# parameter sweep for gate input strength
prot_exp_1D, grate_1D, = perturb_one_param!(ode_problem_dict = ode_problem_dict, 
                                            param_index = 25, 
                                            range_bounds = (0, 4), 
                                            range_size = 50);

# plot protein 2 expresssion, i.e. output, vs induction strength
prot_exp_NOT = plot(prot_exp_1D["parameter_values"],
                    prot_exp_1D["protein_2"],
                    label = false,
                    xlabel = "input (mRNAs/min)",
                    ylabel = "output (# of molecules)",
                    xaxis = :log10)

save_figure(img_to_sv = prot_exp_NOT,
            model_def = model_def,
            custom_suffix = "Figure2_5B(left)")

# plot growth rate vs induction strength
grate_NOT = plot(prot_exp_1D["parameter_values"],
                 grate_1D,
                 label = false,
                 xlabel = "input (mRNAs/min)",
                 ylabel = "growth rate (min\$^{-1}\$)",
                 xaxis = :log10)

save_figure(img_to_sv = grate_NOT,
            model_def = model_def,
            custom_suffix = "Figure2_5B(right)")

#---------------------------------#

# Figure 2.6
# define a host-aware AND gate
model_def = AND_gate_ODE_model!

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

# assemble into single parameter vector
append!(p, het_p)

# need to extend the initial conditions vector
u0_extra = [1, 1, 0, 1, 1, 0]
append!(u0, u0_extra)

# create problem dictionary
ode_problem_dict = create_problem_dict!(model_choice = model_def, 
                                        init_values = u0, 
                                        params_values = p, 
                                        tspan = tspan,
                                        ode_solver = Rodas4(autodiff=false),
                                        abstol = 1e-8,
                                        reltol = 1e-8,
                                        maxiters = 1e7,
                                        show_progress=true)

# parameter sweep for the strength of the two gate inputs
prot_exp_2D, grate_2D = perturb_two_params!(ode_problem_dict = ode_problem_dict, 
                                            param_index_inner = 25, 
                                            param_index_outer = 28, 
                                            range_bounds_inner = (0, 4), 
                                            range_bounds_outer = (0, 4), 
                                            range_size = 50);

# plot protein 3 expresssion, i.e. output, vs induction strength
prot_exp_AND = heatmap(prot_exp_2D["parameter_values_inner"],
                       prot_exp_2D["parameter_values_outer"],
                       vector_to_matrix(prot_exp_2D["protein_3"]), 
                       xlabel = "input 2 (mRNAs/min)",
                       ylabel = "input 1 (mRNAs/min)",
                       xaxis = :log10,
                       yaxis = :log10
                       )

save_figure(img_to_sv = prot_exp_AND,
            model_def = model_def,
            custom_suffix = "Figure2_6B(left)")

# plot growth rate vs induction strength
# plot heatmap
grate_AND = heatmap(prot_exp_2D["parameter_values_inner"],
                    prot_exp_2D["parameter_values_outer"],
                    vector_to_matrix(grate_2D), 
                    xlabel = "input 2 (mRNAs/min)",
                    ylabel = "input 1 (mRNAs/min)",
                    xaxis = :log10,
                    yaxis = :log10
                    )
                    
save_figure(img_to_sv = grate_AND,
            model_def = model_def,
            custom_suffix = "Figure2_6B(right)")

#---------------------------------#

# Figure 2.7
# define a host-aware NAND gate
model_def = NAND_gate_ODE_model!

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

# assemble into single parameter vector
append!(p, het_p)

# need to extend the initial conditions vector
u0_extra = [1, 100, 0, 100, 1, 0, 1, 100, 0]
append!(u0, u0_extra)

# create problem dictionary
ode_problem_dict = create_problem_dict!(model_choice = model_def, 
                                        init_values = u0, 
                                        params_values = p, 
                                        tspan = tspan,
                                        ode_solver = Rodas4(autodiff=false),
                                        abstol = 1e-8,
                                        reltol = 1e-8,
                                        maxiters = 1e7,
                                        show_progress=true)


# parameter sweep for the strength of the two gate inputs
prot_exp_2D, grate_2D = perturb_two_params!(ode_problem_dict = ode_problem_dict, 
                                            param_index_inner = 25, 
                                            param_index_outer = 28, 
                                            range_bounds_inner = (0, 4), 
                                            range_bounds_outer = (0, 4), 
                                            range_size = 50);

# plot protein 4 expresssion, i.e. output, vs induction strength
prot_exp_NAND = heatmap(prot_exp_2D["parameter_values_inner"],
                       prot_exp_2D["parameter_values_outer"],
                       vector_to_matrix(prot_exp_2D["protein_4"]), 
                       xlabel = "input 2 (mRNAs/min)",
                       ylabel = "input 1 (mRNAs/min)",
                       xaxis = :log10,
                       yaxis = :log10
                       );

save_figure(img_to_sv = prot_exp_NAND,
            model_def = model_def,
            custom_suffix = "Figure2_7B(left)")

# plot growth rate vs induction strength
# plot heatmap
grate_NAND = heatmap(prot_exp_2D["parameter_values_inner"],
                     prot_exp_2D["parameter_values_outer"],
                     vector_to_matrix(grate_2D), 
                     xlabel = "input 2 (mRNAs/min)",
                     ylabel = "input 1 (mRNAs/min)",
                     xaxis = :log10,
                     yaxis = :log10
                     );
                    
save_figure(img_to_sv = grate_NAND,
            model_def = model_def,
            custom_suffix = "Figure2_7B(right)")   
#---------------------------------#

# Figure 2.8A :: change RBS values for AND gate output and execute same code
# Figure 2.8B :: change RBS values for NAND gate output and execute same code

# Figure 2.9 :: change nutrient quality values and execute induction strength
#               perturbations for AND and NAND gates with RBS = 1 for gate output



#-------------------#
#     Chapter 3     #
#-------------------# 

# Figure 3.5
model_def = HETER_ODE_model!

kurep = exp10(-2.6575)
kbrep = exp10(-1.3335)

wmaxrep = 150

# organize in a vector
het_p = [ns, dmrep, dprep, kappa_ini, wmaxrep, kbrep, kurep]

# assemble into single parameter vector
append!(p, het_p)

# create problem dictionary
ode_problem_dict = create_problem_dict!(model_choice = model_def, 
                                        init_values = u0, 
                                        params_values = p, 
                                        tspan = tspan,
                                        ode_solver = Rodas4(autodiff=false),
                                        abstol = 1e-8,
                                        reltol = 1e-8,
                                        maxiters = 1e7,
                                        show_progress=true)

# parameter sweep for initiation translation initiation
phet_content, grate, ribosomal_content = trans_initiation!(ode_problem_dict = ode_problem_dict, 
                                                           range_size = 50, 
                                                           kini_lower = -0.65, 
                                                           kini_upper = 0);

# Figure 3.5A :: protein expression versus growth rate
HETER_kappa_ini_prot_v_grate = plot(phet_content["protein_1"], 
                                    grate,
                                    legend = false,
                                    xlabel = "expression (# of molecules)",
                                    ylabel = "growth rate (min\$^{-1}\$)"
                                    )

save_figure(img_to_sv = HETER_kappa_ini_prot_v_grate,
            model_def = model_def,
            custom_suffix = "Figure3_5A")

# Figure 3.5B
plot(phet_content["protein_1"]) #top
plot(grate) #bottom

# Figure 3.5C
plot(phet_content["kappa_ini_values"], ribosomal_content["non_init_complexes_1"])
HETER_kappa_ini_complexes = plot!(phet_content["kappa_ini_values"],
                                  ribosomal_content["init_complexes_1"],
                                  legend = false,
                                  xlabel = "\$κ_{ini}\$",
                                  ylabel = "bound ribosomes (het.)")
save_figure(img_to_sv = HETER_kappa_ini_complexes,
            model_def = model_def,
            custom_suffix = "Figure3_5C")

# Figure 3.5D 
HETER_kappa_ini_endo_complexes = plot(phet_content["kappa_ini_values"],
                                      ribosomal_content["metab_ribo"] +
                                      ribosomal_content["ribo_ribo"] +
                                      ribosomal_content["housek_ribo"] +
                                      ribosomal_content["transf_ribo"],
                                      legend = false,
                                      xlabel = "\$κ_{ini}\$",
                                      ylabel = "bound ribosomes (endo.)")

save_figure(img_to_sv = HETER_kappa_ini_endo_complexes,
            model_def = model_def,
            custom_suffix = "Figure3_5D")
#---------------------------------#

# Figure 3.6

# Figure 3.6A :: incuction strength capacity curve
prot_exp_1D, grate_1D, = perturb_one_param!(ode_problem_dict = ode_problem_dict, 
                                            param_index = 25, 
                                            range_bounds = (0, 4), 
                                            range_size = 50);
                                
HETER_induction_capacity = plot(prot_exp_1D["protein_1"], 
                                grate_1D,
                                legend = false,
                                xlabel = "expression (# of molecules)",
                                ylabel = "growth rate (min\$^{-1}\$)")

save_figure(img_to_sv = HETER_induction_capacity,
            model_def = model_def,
            custom_suffix = "Figure3_6A")

# Figure 3.6B :: RBS strength capacity curve

prot_exp_RBS_only, grate_RBS_only = perturb_RBS!(ode_problem_dict = ode_problem_dict, 
                                                 RBS_bounds = (-4, -2, 0), 
                                                 range_size = 50)

HETER_RBS_capacity = plot(prot_exp_RBS_only["protein_1"], 
                          grate_RBS_only,
                          legend = false,
                          xlabel = "expression (# of molecules)",
                          ylabel = "growth rate (min\$^{-1}\$)")

save_figure(img_to_sv = HETER_RBS_capacity,
            model_def = model_def,
            custom_suffix = "Figure3_6B")

#---------------------------------#

# Figure 3.7

# Figure 3.7B :: induction strength and efficiency of translation initiation
prot_exp_2D, grate_2D = perturb_two_params!(ode_problem_dict = ode_problem_dict, 
                                            param_index_inner = 24, 
                                            param_index_outer = 25, 
                                            range_bounds_inner = (-0.65, 1), 
                                            range_bounds_outer = (0, 4), 
                                            range_size = 50);

# plot results                                            
HETER_kappa_ini_w_induction = plot(prot_exp_2D["protein_1"], 
                                   grate_2D, 
                                   palette = palette([:purple, :green], 50), 
                                   legend = false,
                                   xlabel = "expression (# of molecules)",
                                   ylabel = "growth rate (min\$^{-1}\$)")

save_figure(img_to_sv = HETER_kappa_ini_w_induction,
            model_def = model_def,
            custom_suffix = "Figure3_7B")

# Insets for Fig. 3.7C :: set new induction strength values for the HETER_ODE_model and execute the trans_initiation! function.

#---------------------------------#

# Figure 3.8

# Figure 3.8B :: RBS strength and efficiency of translation initiation
prot_exp_RBS, grate_RBS = perturb_param_w_RBS!(ode_problem_dict = ode_problem_dict, 
                                               param_index = 24, 
                                               range_bounds = (-0.65, 1), 
                                               RBS_bounds = (-4, -2, 0), 
                                               range_size = 50)
                                               
# plot results                                           
HETER_kappa_ini_w_RBS = plot(prot_exp_RBS["protein_1"], 
                             grate_RBS,
                             palette = palette([:purple, :green], 50), 
                             legend = false,
                             xlabel = "expression (# of molecules)",
                             ylabel = "growth rate (min\$^{-1}\$)"
                             )

save_figure(img_to_sv = HETER_kappa_ini_w_RBS,
            model_def = model_def,
            custom_suffix = "Figure3_8B")
        
# Insets for Fig. 3.8C :: set new RBS values for the HETER_ODE_model and execute the trans_initiation! function.