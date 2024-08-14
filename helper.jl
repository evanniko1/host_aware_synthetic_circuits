"""
	create_problem_dict!

Formats all vectors to be used/updated in downstream queries.
Returns a dictionary of:
	1. initial conditions
	2. parameters
	3. time span

Created by Evangelos-Marios Nikolados
"""
function create_problem_dict!(; model_choice, init_values, params_values, tspan, 
							    ode_solver=Rodas4(autodiff=false), 
								abstol=1e-10, reltol=1e-10, maxiters=1e5,
								show_progress=false)

	ode_problem_wrap = Dict("model_def" => model_choice,
							"initial_conditions" => init_values,
							"parameters" => params_values,
							"time_span"  => tspan,
							"ode_solver" => ode_solver,
							"solver_abstol" => abstol,
							"solver_reltol" => reltol,
							"max_iters"     => maxiters,
							"progress"      => show_progress)
	return ode_problem_wrap
end

"""
	calc_growth_rate!

Accepts as input a solved ODE system.
Calculates growth rate for different host-aware synthetic circuits.
Returns values for growth rate and total translation rate.

Created by Evangelos-Marios Nikolados
"""
function calc_growth_rate!(; sol, kappa_ini, gmax, Kgamma, M = 1e8, model_def = HETER_ODE_model!)
	if model_def == HETER_ODE_model!
		ttrate = (sol[end][1] + sol[end][3] + sol[end][4] + sol[end][6] + kappa_ini*sol[end][17])*(gmax*sol[end][14]/(Kgamma + sol[end][14]))
	elseif model_def == REPR_ODE_model!
		ttrate = (sol[end][1] + sol[end][3] + sol[end][4] + sol[end][6] + kappa_ini*sol[end][17] + kappa_ini*sol[end][20] + kappa_ini*sol[end][23]) * (gmax*sol[end][14]/(Kgamma + sol[end][14]))
	elseif model_def == NOT_gate_ODE_model!
		ttrate = (sol[end][1] + sol[end][3] + sol[end][4] + sol[end][6] + kappa_ini*sol[end][17] + kappa_ini*sol[end][20]) * (gmax*sol[end][14]/(Kgamma + sol[end][14]))
	elseif model_def == AND_gate_ODE_model!
		ttrate = (sol[end][1] + sol[end][3] + sol[end][4] + sol[end][6] + kappa_ini*sol[end][17] + kappa_ini*sol[end][20] + kappa_ini*sol[end][23]) * (gmax*sol[end][14]/(Kgamma + sol[end][14]))
	elseif model_def == NAND_gate_ODE_model!
		ttrate = (sol[end][1] + sol[end][3] + sol[end][4] + sol[end][6] + kappa_ini*sol[end][17] + kappa_ini*sol[end][20] + kappa_ini*sol[end][23] + kappa_ini*sol[end][26]) * (gmax*sol[end][14]/(Kgamma + sol[end][14]))
	else
		println("option not supported ...")
	end

	grate = ttrate / M
	
	return ttrate, grate
end

"""
	calc_het_expr!

Accepts as input a solved ODE system.
Returns dictionary of values for protein expression.

Created by Evangelos-Marios Nikolados
"""
function calc_het_expr!(; sol, model_def = HETER_ODE_model!)
	if model_def == HETER_ODE_model!
		protein_expression_dict = Dict("protein_1" => sol[end][15])
	elseif model_def == REPR_ODE_model!
		protein_expression_dict = Dict("protein_1" => sol[end][15],
									   "protein_2" => sol[end][18],
									   "protein_3" => sol[end][21])
	elseif model_def == NOT_gate_ODE_model!
		protein_expression_dict = Dict("protein_1" => sol[end][15],
									   "protein_2" => sol[end][18])
	elseif model_def == AND_gate_ODE_model!
		protein_expression_dict = Dict("protein_1" => sol[end][15],
									   "protein_2" => sol[end][18],
									   "protein_3" => sol[end][21])
	elseif model_def == NAND_gate_ODE_model!
		protein_expression_dict = Dict("protein_1" => sol[end][15],
									   "protein_2" => sol[end][18],
									   "protein_3" => sol[end][21],
									   "protein_4" => sol[end][24])
	else
		println("something went wrong ... check model definition again ...")
		return 1
	end

	return protein_expression_dict
end

"""
	calc_het_expr!

Accepts as input a solved ODE system.
Returns dictionary of values for ribosomal content.

Created by Evangelos-Marios Nikolados
"""
function calc_ribo_content!(; sol, kappa_ini, model_def = HETER_ODE_model!)

	# heterologous complexes
	if model_def == HETER_ODE_model!
		complexes_dict = Dict("init_complexes_1"     => kappa_ini * sol[end][17],
							  "non_init_complexes_1" => (1 - kappa_ini) * sol[end][17])
	elseif model_def == REPR_ODE_model!
		complexes_dict = Dict("init_complexes_1"     => kappa_ini * sol[end][17],
							  "non_init_complexes_1" => (1 - kappa_ini) * sol[end][17],
							  "init_complexes_2"     => kappa_ini * sol[end][20],
							  "non_init_complexes_2" => (1 - kappa_ini) * sol[end][20],
							  "init_complexes_3"     => kappa_ini * sol[end][23],
							  "non_init_complexes_3" => (1 - kappa_ini) * sol[end][23])
	elseif model_def == NOT_gate_ODE_model!
		complexes_dict = Dict("init_complexes_1"     => kappa_ini * sol[end][17],
							  "non_init_complexes_1" => (1 - kappa_ini) * sol[end][17],
							  "init_complexes_2"     => kappa_ini * sol[end][20],
							  "non_init_complexes_2" => (1 - kappa_ini) * sol[end][20])
	elseif model_def == AND_gate_ODE_model!
		complexes_dict = Dict("init_complexes_1"     => kappa_ini * sol[end][17],
							  "non_init_complexes_1" => (1 - kappa_ini) * sol[end][17],
							  "init_complexes_2"     => kappa_ini * sol[end][20],
							  "non_init_complexes_2" => (1 - kappa_ini) * sol[end][20],
							  "init_complexes_3"     => kappa_ini * sol[end][23],
							  "non_init_complexes_3" => (1 - kappa_ini) * sol[end][23])
	elseif model_def == NAND_gate_ODE_model!
		complexes_dict = Dict("init_complexes_1"     => kappa_ini * sol[end][17],
							  "non_init_complexes_1" => (1 - kappa_ini) * sol[end][17],
							  "init_complexes_2"     => kappa_ini * sol[end][20],
							  "non_init_complexes_2" => (1 - kappa_ini) * sol[end][20],
							  "init_complexes_3"     => kappa_ini * sol[end][23],
							  "non_init_complexes_3" => (1 - kappa_ini) * sol[end][23],
							  "init_complexes_4"	 => kappa_ini * sol[end][26],
							  "non_init_complexes_4" => (1 - kappa_ini) * sol[end][26])
	else
		println("something went wrong ... check model definition again ...")
		return 1
	end

	# endogenous complexes
	complexes_dict["ribo_complexes"]   = sol[end][1]
	complexes_dict["housek_complexes"] = sol[end][3]
	complexes_dict["metab_complexes"]  = sol[end][4]
	complexes_dict["trans_complexes"]  = sol[end][6]
	complexes_dict["free_ribo"]        = sol[end][13]

	return complexes_dict
end

"""
	solve_ode_problem!

Creates the ODEProblem object and passes it to the solver.
Returns a sol object containing the solved ODE system.

Created by Evangelos-Marios Nikolados
"""
function solve_ode_problem!(; ode_problem_wrap)
	prob = ODEProblem(ode_problem_wrap["model_def"], ode_problem_wrap["initial_conditions"], ode_problem_wrap["time_span"], ode_problem_wrap["parameters"])
	sol = solve(prob, ode_problem_wrap["ode_solver"], abstol=ode_problem_wrap["solver_abstol"],reltol=ode_problem_wrap["solver_reltol"], maxiters=ode_problem_wrap["max_iters"], progress=ode_problem_wrap["progress"], isoutofdomain = (m,p,t) -> any(x->x<0, m))

	return sol
end

"""
    trans_initiation!(; init_values, tspan, params_values, range_size=10, kini_lower=-0.65, kini_upper=0)

Helper function to replicate the non-linear expression-growth relationship from Cambray et al, 2018.
By default, the function will simulate 10 models for log-sampled κ_ini values in the range [10^-0.65, 1].

Returns two vectors:
     i.  heterologous protein expression values
     ii. growth rate values

Created by Evangelos-Marios Nikolados.
"""
function trans_initiation!(; ode_problem_dict, range_size = 10, kini_lower = -0.65, kini_upper = 0)

    grate_sols, param_vals = [], []
	# lists for ribosomal (bound & free) content
	# heterologous
	if ode_problem_dict["model_def"] == HETER_ODE_model!
		phet_sols_1, init_sols_1, non_init_sols_1 = [], [], []
	elseif ode_problem_dict["model_def"] == REPR_ODE_model!
		phet_sols_1, init_sols_1, non_init_sols_1 = [], [], []
		phet_sols_2, init_sols_2, non_init_sols_2 = [], [], []
		phet_sols_3, init_sols_3, non_init_sols_3 = [], [], []
	elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
		phet_sols_1, init_sols_1, non_init_sols_1 = [], [], []
		phet_sols_2, init_sols_2, non_init_sols_2 = [], [], []
	elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
		phet_sols_1, init_sols_1, non_init_sols_1 = [], [], []
		phet_sols_2, init_sols_2, non_init_sols_2 = [], [], []
		phet_sols_3, init_sols_3, non_init_sols_3 = [], [], []
	elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
		phet_sols_1, init_sols_1, non_init_sols_1 = [], [], []
		phet_sols_2, init_sols_2, non_init_sols_2 = [], [], []
		phet_sols_3, init_sols_3, non_init_sols_3 = [], [], []
		phet_sols_4, init_sols_4, non_init_sols_4 = [], [], []
	else
	end
	# endogenous bound
	ribo_ribo_sols, housek_ribo_sols, metab_ribo_sols, transf_ribo_sols = [], [], [], []
	# free
	free_ribo_sols = []

    for (_, kappa_ini) in enumerate(exp10.(range(kini_lower, kini_upper, length=range_size)))
        # update kappa_ini value
		ode_problem_dict["parameters"][24] = kappa_ini
		push!(param_vals, kappa_ini)
    
        # define & solve the new ODE problem
		sol = solve_ode_problem!(ode_problem_wrap = ode_problem_dict)

		# calculate relevant rates
		_, grate = calc_growth_rate!(sol = sol, kappa_ini = kappa_ini, gmax = gmax, Kgamma = Kgamma, model_def = ode_problem_dict["model_def"])
		het_expr_dict = calc_het_expr!(sol = sol, model_def = ode_problem_dict["model_def"])
		complexes_dict = calc_ribo_content!(sol = sol, kappa_ini = kappa_ini, model_def = ode_problem_dict["model_def"])

		# CONDITIONALS
        # push what we need
		if ode_problem_dict["model_def"] == HETER_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(init_sols_1, complexes_dict["init_complexes_1"])
			push!(non_init_sols_1, complexes_dict["non_init_complexes_1"])
		elseif ode_problem_dict["model_def"] == REPR_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(init_sols_1, complexes_dict["init_complexes_1"])
			push!(non_init_sols_1, complexes_dict["non_init_complexes_1"])

			push!(phet_sols_2, het_expr_dict["protein_2"])
			push!(init_sols_2, complexes_dict["init_complexes_2"])
			push!(non_init_sols_2, complexes_dict["non_init_complexes_2"])

			push!(phet_sols_3, het_expr_dict["protein_3"])
			push!(init_sols_3, complexes_dict["init_complexes_3"])
			push!(non_init_sols_3, complexes_dict["non_init_complexes_3"])
		elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(init_sols_1, complexes_dict["init_complexes_1"])
			push!(non_init_sols_1, complexes_dict["non_init_complexes_1"])

			push!(phet_sols_2, het_expr_dict["protein_2"])
			push!(init_sols_2, complexes_dict["init_complexes_2"])
			push!(non_init_sols_2, complexes_dict["non_init_complexes_2"])
		elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(init_sols_1, complexes_dict["init_complexes_1"])
			push!(non_init_sols_1, complexes_dict["non_init_complexes_1"])

			push!(phet_sols_2, het_expr_dict["protein_2"])
			push!(init_sols_2, complexes_dict["init_complexes_2"])
			push!(non_init_sols_2, complexes_dict["non_init_complexes_2"])

			push!(phet_sols_3, het_expr_dict["protein_3"])
			push!(init_sols_3, complexes_dict["init_complexes_3"])
			push!(non_init_sols_3, complexes_dict["non_init_complexes_3"])
		elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(init_sols_1, complexes_dict["init_complexes_1"])
			push!(non_init_sols_1, complexes_dict["non_init_complexes_1"])

			push!(phet_sols_2, het_expr_dict["protein_2"])
			push!(init_sols_2, complexes_dict["init_complexes_2"])
			push!(non_init_sols_2, complexes_dict["non_init_complexes_2"])

			push!(phet_sols_3, het_expr_dict["protein_3"])
			push!(init_sols_3, complexes_dict["init_complexes_3"])
			push!(non_init_sols_3, complexes_dict["non_init_complexes_3"])

			push!(phet_sols_4, het_expr_dict["protein_4"])
			push!(init_sols_4, complexes_dict["init_complexes_4"])
			push!(non_init_sols_4, complexes_dict["non_init_complexes_4"])			
		else
			println("something went wrong ...")
		end
        
        push!(grate_sols, grate)
		push!(ribo_ribo_sols, complexes_dict["ribo_complexes"])
		push!(housek_ribo_sols, complexes_dict["housek_complexes"])
		push!(metab_ribo_sols, complexes_dict["metab_complexes"])
		push!(transf_ribo_sols, complexes_dict["trans_complexes"])
		push!(free_ribo_sols, complexes_dict["free_ribo"])

    end
	# organize all lists of solutions for mRNA-ribosomal complexes in a dict

	if ode_problem_dict["model_def"] == HETER_ODE_model!
		ribosomal_content = Dict("init_complexes_1"     => init_sols_1,
								 "non_init_complexes_1" => non_init_sols_1)

		het_protein_content = Dict("protein_1" => phet_sols_1)

	elseif ode_problem_dict["model_def"] == REPR_ODE_model!
		ribosomal_content = Dict("init_complexes_1"     => init_sols_1,
								 "non_init_complexes_1" => non_init_sols_1,
								 "init_complexes_2"     => init_sols_2,
								 "non_init_complexes_2" => non_init_sols_2,
								 "init_complexes_3"     => init_sols_3,
								 "non_init_complexes_3" => non_init_sols_3)	
								 
		het_protein_content = Dict("protein_1" => phet_sols_1,
								"protein_2" => phet_sols_2,
								"protein_3" => phet_sols_3)	

	elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
				ribosomal_content = Dict("init_complexes_1"     => init_sols_1,
								 "non_init_complexes_1" => non_init_sols_1,
								 "init_complexes_2"     => init_sols_2,
								 "non_init_complexes_2" => non_init_sols_2)	

		het_protein_content = Dict("protein_1" => phet_sols_1,
								"protein_2" => phet_sols_2)	

	elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
				ribosomal_content = Dict("init_complexes_1"     => init_sols_1,
								 "non_init_complexes_1" => non_init_sols_1,
								 "init_complexes_2"     => init_sols_2,
								 "non_init_complexes_2" => non_init_sols_2,
								 "init_complexes_3"     => init_sols_3,
								 "non_init_complexes_3" => non_init_sols_3)	

		het_protein_content = Dict("protein_1" => phet_sols_1,
								"protein_2" => phet_sols_2,
								"protein_3" => phet_sols_3)	

	elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
				ribosomal_content = Dict("init_complexes_1"     => init_sols_1,
								 "non_init_complexes_1" => non_init_sols_1,
								 "init_complexes_2"     => init_sols_2,
								 "non_init_complexes_2" => non_init_sols_2,
								 "init_complexes_3"     => init_sols_3,
								 "non_init_complexes_3" => non_init_sols_3,
								 "init_complexes_4"     => init_sols_4,
								 "non_init_complexes_4" => non_init_sols_4)	

		het_protein_content = Dict("protein_1" => phet_sols_1,
								"protein_2" => phet_sols_2,
								"protein_3" => phet_sols_3,
								"protein_4" => phet_sols_4)	
	else
		println("something went wrong ...")
	end

	het_protein_content["kappa_ini_values"] = param_vals

	ribosomal_content["ribo_ribo"]   = ribo_ribo_sols
	ribosomal_content["housek_ribo"] = housek_ribo_sols
	ribosomal_content["transf_ribo"] = transf_ribo_sols
	ribosomal_content["metab_ribo"]  = metab_ribo_sols
	ribosomal_content["free_ribo"]   = free_ribo_sols

	return het_protein_content, grate_sols, ribosomal_content
end

"""
	perturb_one_param!(; ode_problem_dict, param_index, range_bounds, range_size)

Helper function to perturb a single parameter of the model. Parameter index specifies the parameter to perturb,
range bounds are then log10() transformed, and range_size specifies log-spaced sampling steps.

Returns three vectors:
     i.   heterologous protein expression
     ii.  growth rate
	 iii. total translation rate

Created by Evangelos-Marios Nikolados.
"""
function perturb_one_param!(; ode_problem_dict, param_index, range_bounds, range_size=10)
	# initialize phenotype and growth rate result vectors

    grate_sols, param_vals = [], []

	# heterologous
	if ode_problem_dict["model_def"] == HETER_ODE_model!
		phet_sols_1 = []
	elseif ode_problem_dict["model_def"] == REPR_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3 = [], [], []
	elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
		phet_sols_1, phet_sols_2 = [], [], []
	elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3 = [], [], []
	elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3, phet_sols_4 = [], [], [], []
	else
	end
    
    for (_, param_to_perturb) in enumerate(exp10.(range(range_bounds[1], range_bounds[2], length=range_size)))
        # update value for parameter
		ode_problem_dict["parameters"][param_index] = param_to_perturb
		# keep track of parameter values -- use later to plot axes
		push!(param_vals, param_to_perturb)
    
        # define & solve the new ODE problem
		sol = solve_ode_problem!(ode_problem_wrap = ode_problem_dict)

		# calculate relevant rates
		_, grate = calc_growth_rate!(sol = sol, kappa_ini = ode_problem_dict["parameters"][24], gmax = gmax, Kgamma = Kgamma, model_def = ode_problem_dict["model_def"])
		het_expr_dict = calc_het_expr!(sol = sol, model_def = ode_problem_dict["model_def"])

        # push what we need
		if ode_problem_dict["model_def"] == HETER_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])

		elseif ode_problem_dict["model_def"] == REPR_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(phet_sols_2, het_expr_dict["protein_2"])
			push!(phet_sols_3, het_expr_dict["protein_3"])

		elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(phet_sols_2, het_expr_dict["protein_2"])

		elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(phet_sols_2, het_expr_dict["protein_2"])
			push!(phet_sols_3, het_expr_dict["protein_3"])

		elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(phet_sols_2, het_expr_dict["protein_2"])
			push!(phet_sols_3, het_expr_dict["protein_3"])
			push!(phet_sols_4, het_expr_dict["protein_4"])
		
		else
			println("something went wrong ...")
		end
        
        push!(grate_sols, grate)

    end

	if ode_problem_dict["model_def"] == HETER_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1)

	elseif ode_problem_dict["model_def"] == REPR_ODE_model!
		
		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3)	

	elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2)	

	elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3)	

	elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3,
								   "protein_4" => phet_sols_4)	
	else
		println("something went wrong ...")
	end

	# add parameter values
	het_protein_content["parameter_values"] = param_vals

	return het_protein_content, grate_sols
end

"""
	perturb_two_params!(; ode_problem_dict, param_index_inner, param_index_outer, range_bounds_inner, range_bounds_outer, range_size)

Helper function to perturb two parameters of the model. Parameter index specifies the parameter to perturb, for inner and outer loop respectively,
range bounds are then log10() transformed, for inner and outer ranges respectively, and range_size specifies log-spaced sampling steps.

Returns three vectors of vectors:
     i.   heterologous protein expression
     ii.  growth rate
	 iii. total translation rate

Created by Evangelos-Marios Nikolados.
"""
function perturb_two_params!(; ode_problem_dict, param_index_inner, param_index_outer, range_bounds_inner, range_bounds_outer, range_size = 10)
	# initialize phenotype and growth rate result vectors
	grate_sols, param_vals_inner, param_vals_outer = [], [], []
	# heterologous
	if ode_problem_dict["model_def"] == HETER_ODE_model!
		phet_sols_1 = []
	elseif ode_problem_dict["model_def"] == REPR_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3 = [], [], []
	elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
		phet_sols_1, phet_sols_2 = [], [], []
	elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3 = [], [], []
	elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3, phet_sols_4 = [], [], [], []
	else
	end

	for (_, outer_param_to_perturb) in enumerate(exp10.(range(range_bounds_outer[1], range_bounds_outer[2], length = range_size)))
		# update value for outer parameter loop
		ode_problem_dict["parameters"][param_index_outer] = outer_param_to_perturb
		push!(param_vals_outer, outer_param_to_perturb)

		# INNER LOOP
		het_expr_dict, grate_pert_one = perturb_one_param!(ode_problem_dict = ode_problem_dict, param_index = param_index_inner, range_bounds = range_bounds_inner, range_size = range_size)
        
		# push what we need
		if ode_problem_dict["model_def"] == HETER_ODE_model!
			push!(phet_sols_1, values(het_expr_dict["protein_1"]))

		elseif ode_problem_dict["model_def"] == REPR_ODE_model!
			push!(phet_sols_1, values(het_expr_dict["protein_1"]))
			push!(phet_sols_2, values(het_expr_dict["protein_2"]))
			push!(phet_sols_3, values(het_expr_dict["protein_3"]))

		elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
			push!(phet_sols_1, values(het_expr_dict["protein_1"]))
			push!(phet_sols_2, values(het_expr_dict["protein_2"]))

		elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
			push!(phet_sols_1, values(het_expr_dict["protein_1"]))
			push!(phet_sols_2, values(het_expr_dict["protein_2"]))
			push!(phet_sols_3, values(het_expr_dict["protein_3"]))

		elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
			push!(phet_sols_1, values(het_expr_dict["protein_1"]))
			push!(phet_sols_2, values(het_expr_dict["protein_2"]))
			push!(phet_sols_3, values(het_expr_dict["protein_3"]))
			push!(phet_sols_4, values(het_expr_dict["protein_4"]))
		
		else
			println("something went wrong ...")
		end

        push!(grate_sols, grate_pert_one)
		push!(param_vals_inner, values(het_expr_dict["parameter_values"]))
	end

	# organize heterolous expression in a dictionary
	if ode_problem_dict["model_def"] == HETER_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1)

	elseif ode_problem_dict["model_def"] == REPR_ODE_model!
		
		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3)	

	elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2)	

	elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3)	

	elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3,
								   "protein_4" => phet_sols_4)	
	else
		println("something went wrong ...")
	end
	
	het_protein_content["parameter_values_inner"] = param_vals_inner[1]
	het_protein_content["parameter_values_outer"] = param_vals_outer
	return het_protein_content, grate_sols
end

"""
	perturb__RBS!(; ode_problem_dict, RBS_bounds, range_size)

Helper function to perturb the RBS of a given mRNA, defined as the ratio of binding and unbinding rates.
RBS bounds are used to determing the lower and upper ranges for the binding and unbinding rates.
Then, they are log10() transformed. Finally, range_size specifies log-spaced sampling steps.

Returns three vectors of vectors:
     i.   heterologous protein expression
     ii.  growth rate
	 iii. total translation rate

Created by Evangelos-Marios Nikolados.
"""
function perturb_RBS!(; ode_problem_dict, RBS_bounds, range_size = 10)
    # initialize phenotype and growth rate results vectors

	grate_sols = []
	# heterologous
	if ode_problem_dict["model_def"] == HETER_ODE_model!
		phet_sols_1 = []
	elseif ode_problem_dict["model_def"] == REPR_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3 = [], [], []
	elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
		phet_sols_1, phet_sols_2 = [], [], []
	elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3 = [], [], []
	elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3, phet_sols_4 = [], [], [], []
	else
	end
    # generate RBS bounds :: RBS is degined as the ratio kb/ku
    kb_rng = exp10.(range(RBS_bounds[2], RBS_bounds[3], length = range_size))
    ku_rng = reverse(exp10.(range(RBS_bounds[1], RBS_bounds[2], length = range_size)))

    for (kbrep_v, kurep_v) in zip(kb_rng, ku_rng)
        # assign kb and ku values for the heterologous reporter protein
		if ode_problem_dict["model_def"] == HETER_ODE_model!
			ode_problem_dict["parameters"][26], ode_problem_dict["parameters"][27] = kbrep_v, kurep_v
		elseif ode_problem_dict["model_def"] == REPR_ODE_model!
			ode_problem_dict["parameters"][26], ode_problem_dict["parameters"][27] = kbrep_v, kurep_v
			ode_problem_dict["parameters"][29], ode_problem_dict["parameters"][30] = kbrep_v, kurep_v
			ode_problem_dict["parameters"][32], ode_problem_dict["parameters"][33] = kbrep_v, kurep_v
		elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
			ode_problem_dict["parameters"][29], ode_problem_dict["parameters"][30] = kbrep_v, kurep_v
		elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
			ode_problem_dict["parameters"][32], ode_problem_dict["parameters"][33] = kbrep_v, kurep_v
		elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
			ode_problem_dict["parameters"][35], ode_problem_dict["parameters"][36] = kbrep_v, kurep_v
		else
			println("something went wrong ... check the model specified")
		end

        # define & solve the new ODE problem
		sol = solve_ode_problem!(ode_problem_wrap = ode_problem_dict)

		# calculate relevant rates
		_, grate = calc_growth_rate!(sol = sol, kappa_ini = ode_problem_dict["parameters"][24], gmax = gmax, Kgamma = Kgamma, model_def = ode_problem_dict["model_def"])
		het_expr_dict = calc_het_expr!(sol = sol, model_def = ode_problem_dict["model_def"])

		# CONDITIONALS
        # push what we need
		if ode_problem_dict["model_def"] == HETER_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
		elseif ode_problem_dict["model_def"] == REPR_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(phet_sols_2, het_expr_dict["protein_2"])
			push!(phet_sols_3, het_expr_dict["protein_3"])

		elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(phet_sols_2, het_expr_dict["protein_2"])

		elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(phet_sols_2, het_expr_dict["protein_2"])
			push!(phet_sols_3, het_expr_dict["protein_3"])

		elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
			push!(phet_sols_1, het_expr_dict["protein_1"])
			push!(phet_sols_2, het_expr_dict["protein_2"])
			push!(phet_sols_3, het_expr_dict["protein_3"])
			push!(phet_sols_4, het_expr_dict["protein_4"])
		else
			println("something went wrong ...")
		end
        
        push!(grate_sols, grate)
	end

	if ode_problem_dict["model_def"] == HETER_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1)

	elseif ode_problem_dict["model_def"] == REPR_ODE_model!
		
		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3)	

	elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2)	

	elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3)	

	elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3,
								   "protein_4" => phet_sols_4)	
	else
		println("something went wrong ...")
	end

	return het_protein_content, grate_sols
end

"""
	perturb_param_w_RBS!(; ode_problem_dict, param_index, range_bounds, RBS_bounds, range_size)

Helper function to perturb RBS and one user-specified parameter of the model. Parameter index specifies the parameter to perturb,
range bounds are then log10() transformed, and range_size specifies log-spaced sampling steps. The same apply to RBS.

Returns three vectors of vectors:
     i.   heterologous protein expression
     ii.  growth rate
	 iii. total translation rate

Created by Evangelos-Marios Nikolados.
"""
function perturb_param_w_RBS!(; ode_problem_dict, param_index, range_bounds, RBS_bounds, range_size = 10)
    # initialize phenotype and growth rate results vectors

	grate_sols = []
	# heterologous
	if ode_problem_dict["model_def"] == HETER_ODE_model!
		phet_sols_1 = []
	elseif ode_problem_dict["model_def"] == REPR_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3 = [], [], []
	elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
		phet_sols_1, phet_sols_2 = [], [], []
	elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3 = [], [], []
	elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
		phet_sols_1, phet_sols_2, phet_sols_3, phet_sols_4 = [], [], [], []
	else
	end

    # generate RBS bounds :: RBS is degined as the ratio kb/ku
    kb_rng = exp10.(range(RBS_bounds[2], RBS_bounds[3], length = range_size))
    ku_rng = reverse(exp10.(range(RBS_bounds[1], RBS_bounds[2], length = range_size)))

    for (kbrep_v, kurep_v) in zip(kb_rng, ku_rng)
        # assign kb and ku values for the heterologous reporter protein
		if ode_problem_dict["model_def"] == HETER_ODE_model!
			ode_problem_dict["parameters"][26], ode_problem_dict["parameters"][27] = kbrep_v, kurep_v
		elseif ode_problem_dict["model_def"] == REPR_ODE_model!
			ode_problem_dict["parameters"][26], ode_problem_dict["parameters"][27] = kbrep_v, kurep_v
			ode_problem_dict["parameters"][29], ode_problem_dict["parameters"][30] = kbrep_v, kurep_v
			ode_problem_dict["parameters"][32], ode_problem_dict["parameters"][33] = kbrep_v, kurep_v
		elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
			ode_problem_dict["parameters"][29], ode_problem_dict["parameters"][30] = kbrep_v, kurep_v
		elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
			ode_problem_dict["parameters"][32], ode_problem_dict["parameters"][33] = kbrep_v, kurep_v
		elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
			ode_problem_dict["parameters"][35], ode_problem_dict["parameters"][36] = kbrep_v, kurep_v
		else
			println("something went wrong ... check the model specified")
		end
        # INNER LOOP :: perturb user specified parameter
		het_expr_dict, grate_pert_one = perturb_one_param!(ode_problem_dict = ode_problem_dict, param_index = param_index, range_bounds = range_bounds, range_size = range_size)
        
		# push what we need
		if ode_problem_dict["model_def"] == HETER_ODE_model!
			push!(phet_sols_1, values(het_expr_dict["protein_1"]))

		elseif ode_problem_dict["model_def"] == REPR_ODE_model!
			push!(phet_sols_1, values(het_expr_dict["protein_1"]))
			push!(phet_sols_2, values(het_expr_dict["protein_2"]))
			push!(phet_sols_3, values(het_expr_dict["protein_3"]))

		elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!
			push!(phet_sols_1, values(het_expr_dict["protein_1"]))
			push!(phet_sols_2, values(het_expr_dict["protein_2"]))

		elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!
			push!(phet_sols_1, values(het_expr_dict["protein_1"]))
			push!(phet_sols_2, values(het_expr_dict["protein_2"]))
			push!(phet_sols_3, values(het_expr_dict["protein_3"]))

		elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!
			push!(phet_sols_1, values(het_expr_dict["protein_1"]))
			push!(phet_sols_2, values(het_expr_dict["protein_2"]))
			push!(phet_sols_3, values(het_expr_dict["protein_3"]))
			push!(phet_sols_4, values(het_expr_dict["protein_4"]))
		
		else
			println("something went wrong ...")
		end

        push!(grate_sols, grate_pert_one)

	end
	
	# organize heterolous expression in a dictionary
	if ode_problem_dict["model_def"] == HETER_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1)

	elseif ode_problem_dict["model_def"] == REPR_ODE_model!
		
		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3)	

	elseif ode_problem_dict["model_def"] == NOT_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2)	

	elseif ode_problem_dict["model_def"] == AND_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3)	

	elseif ode_problem_dict["model_def"] == NAND_gate_ODE_model!

		het_protein_content = Dict("protein_1" => phet_sols_1,
								   "protein_2" => phet_sols_2,
								   "protein_3" => phet_sols_3,
								   "protein_4" => phet_sols_4)	
	else
		println("something went wrong ...")
	end
	
	return het_protein_content, grate_sols
end

"""
	vector_to_matrix

Helper function to convert a vector to matrix.

Created by Evangelos-Marios Nikolados.
"""
function vector_to_matrix(vector)
    return reduce(hcat, vector)
end

"""
	save_figure

Helper function to save a generated figure at the specified path.

Created by Evangelos-Marios Nikolados.
"""
function save_figure(; img_to_sv, model_def, custom_suffix = "", path_to_sv = "./figures/")
    if custom_suffix == ""
        savefig(img_to_sv, sv_path * string(model_def)[1:end-11] * "_" * string(length(readdir(path_to_sv))+1))
    else
        if typeof(custom_suffix) == String
            savefig(img_to_sv, sv_path * string(model_def)[1:end-11] * "_" * custom_suffix)
        else
            savefig(img_to_sv, sv_path * string(model_def)[1:end-11] * "_" * string(custom_suffix))
        end
    end
end