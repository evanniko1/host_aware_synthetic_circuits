function create_problem_dict!(; init_values, params_values, tspan)
	ode_problem_wrap = Dict("initial_conditions" => init_values,
							"parameters" => params_values,
							"time_span"  => tspan)
	return ode_problem_wrap
end


"""
    ODE_model!(du,u,p,t)

Defines ODEs for extension of the growth model by Weisse et al, 2015, 
with effects of translation initiation on the expression of a heterologous protein.

Created by Evangelos-Marios Nikolados.
"""
function ODE_model!(du,u,p,t)

	thetar= p[1] # ribosomal threshold amount of energy at which transcription is half maximal 
	k_cm= p[2] # chloramphenicol-binding rate no idea why it is here
	s0= p[3] # external nutrient
	gmax= p[4] # maximal rate of translational elongation 
	cl= p[5] #
	thetax= p[6] # non-ribosomal threshold amount of energy at which transcription is half maximal
	Kt= p[7] # nutrient import threshold
	M= p[8] # fixed size of the cell i.e total cell mass
	we= p[9] # max. enzyme transcription rate we = wt = wm
	Km= p[10] # enzymatic threshold 
    f= cl*k_cm

	vm= p[11] # max. enzymatic rate
	nx= p[12] # length of non-ribosomal proteins
	Kq= p[13] # q-autoinhibition threshold (some constant for translation rate of q)
	Kp= p[14] # kgamma = gmax/Kp = 7
	vt= p[15] # max. nutrient import rate (something to do with nuimpt)
	wr= p[16] # maximal transcription rate of r
	wq= p[17] # maximal transcription rate of q 
	wp= p[18] # maximal transcription rate of p
	nq= p[19] # q-autoinhibition Hill coef (?)
	nr= p[20] # length of r (rybosome) protein
	ns= p[21] # one molecule of s yields ns molecules of a (energy) (nutrient efficiency)

    b= 0  # Hack - this variables also should be passed via p rather than defined here but for some reason it didn't work
    dm= 0.1
    kb= 0.0095
    ku= 1

    wmaxrep = p[22]
    kbrep= p[23]
    kurep= p[24]
    dmrep = p[25]
    dprep = p[26]
    kappa_ini = p[27]
    
	rmr= u[1] # complex between a ribosome and the mRNA for protein r
	em= u[2] # enzyme that metabolizes si (nutrinet) into a (energy)
	rmq= u[3] # complex between a ribosome and the mRNA for protein q
	rmt= u[4] # complex between a ribosome and the mRNA for protein et
	et= u[5] # enzyme that transports s (energy) into the cell
	rmm= u[6] # complex between a ribosome and the mRNA for protein em
	mt= u[7] # mRNA of et
	mm= u[8] # mRNA of em
	q= u[9] # house keeping proteins
	si= u[10] # nutrient (internalized)
	mq= u[11] # mRNA of q
	mr= u[12] # mRNA of r (rybosomes)
	r= u[13] # number of rybosomes
	a= u[14] # energy
    
	# rep is heterologous protein
    rep= u[15]
    mrep = u[16]
    rmrep = u[17]

	Kgamma= gmax/Kp # threshold amount of energy where elongation is halfmaximal
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation
	ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep)*gamma # term in expression of lambda
	lam= ttrate/M # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of si by em

    vrep = rmrep/nx*gamma*kappa_ini;
    wrep = wmaxrep*a/(thetax+a);
 
	du[1]= kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr # (6)
	du[2]= gamma/nx*rmm-lam*em # (4)
	du[3]= kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq # (6)
	du[4]= kb*r*mt+b-ku*rmt-gamma/nx*rmt-lam*rmt # (6)
	du[5]= gamma/nx*rmt-lam*et # (4)
	du[6]= kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm # (6)
	du[7]= (we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt # (5)
	du[8]= (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm # (5)
	du[9]= gamma/nx*rmq-lam*q # (4)
	du[10]= (et*vt*s0/(Kt + s0))-nucat-lam*si # (1)  import and metabolism of si
	du[11]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq # (5)
	du[12]= (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr # (5)
	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep*rmrep+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep*r*mrep-lam*r # (3) gamma/nr*rmr translation rate of r i.e. mur
	du[14]= ns*nucat-ttrate-lam*a # (2) energy is created from metabolism and lost from translation and dilution

	# heterologous equations, equations above also required modifications for heterologous case!
    du[15]= vrep - (lam + dprep)*rep
    du[16]= wrep - (lam + dmrep)*mrep + vrep - kbrep*r*mrep + (1-kappa_ini)*kurep*rmrep
    du[17]= -lam*rmrep + kbrep*r*mrep - (1-kappa_ini)*kurep*rmrep - vrep
end

function calc_growth_rate!(; sol, kappa_ini, gmax, Kgamma, M = 1e8)
	ttrate = ttrate = (sol[end][1] + sol[end][3] + sol[end][4] + sol[end][6] + kappa_ini*sol[end][17])*(gmax*sol[end][14]/(Kgamma + sol[end][14]))
	grate = ttrate / M
	
	return ttrate, grate
end

function calc_het_expr!(; sol)
	return sol[end][15]
end

function solve_ode_problem!(; model_def, ode_problem_wrap, solver = Rodas4())
	prob = ODEProblem(model_def, ode_problem_wrap["initial_conditions"], ode_problem_wrap["time_span"], ode_problem_wrap["parameters"])
	sol = solve(prob, solver)

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

    phet_sols, grate_sols = [], []
    
    for (idx, kappa_ini) in enumerate(exp10.(range(kini_lower, kini_upper, length=range_size)))
        # update kappa_ini value
		ode_problem_dict["parameters"][end] = kappa_ini
    
        # define & solve the new ODE problem
		sol = solve_ode_problem!(model_def = ODE_model!, ode_problem_wrap = ode_problem_dict)

		# calculate relevant rates
		_, grate = calc_growth_rate!(sol = sol, kappa_ini = kappa_ini, gmax = gmax, Kgamma = Kgamma)
		het_expr = calc_het_expr!(sol = sol)

        # push what we need
        push!(phet_sols, het_expr)
        push!(grate_sols, grate)
    end

    return phet_sols, grate_sols
end

########
# NEW #
########

# abstract parameter perturbations to avoid redundant code
function perturb_one_param!(; ode_problem_dict, param_index, range_bounds, range_size)
	# initialize phenotype, growth rate, and total translation rate result vectors
    phet_sols, grate_sols, ttrate_sols = [], [], []
    
    for (_, param_to_perturb) in enumerate(exp10.(range(range_bounds[1], range_bounds[2], length=range_size)))
        # update value for parameter
		ode_problem_dict["parameters"][param_index] = param_to_perturb
    
        # define & solve the new ODE problem
		sol = solve_ode_problem!(model_def = ODE_model!, ode_problem_wrap = ode_problem_dict)

		# calculate relevant rates
		ttrate, grate = calc_growth_rate!(sol = sol, kappa_ini = param_to_perturb, gmax = gmax, Kgamma = Kgamma)
		het_expr = calc_het_expr!(sol = sol)

        # push what we need
        push!(phet_sols, het_expr)
        push!(grate_sols, grate)
		push!(ttrate_sols, ttrate)
    end

    return phet_sols, grate_sols, ttrate_sols
end

function perturb_two_params!(; ode_problem_dict, param_index_inner, param_index_outer, range_bounds_inner, range_bounds_outer, range_size = 10)
	# initialize phenotype, growth rate, and total translation rate result vectors
	phet_sols, grate_sols, ttrate_sols = [], [], []

	for (_, outer_param_to_perturb) in enumerate(exp10.(range(range_bounds_outer[1], range_bounds_outer[2], length = range_size)))
		# update value for outer parameter loop
		ode_problem_dict["parameters"][param_index_outer] = outer_param_to_perturb

		# INNER LOOP
		phet_pert_one, grate_pert_one, ttrate_pert_one = perturb_one_param!(ode_problem_dict = ode_problem_dict, param_index = param_index_inner, range_bounds = range_bounds_inner, range_size = range_size)
        
		# push what we need
        push!(phet_sols, phet_pert_one)
        push!(grate_sols, grate_pert_one)
		push!(ttrate_sols, ttrate_pert_one)
	end

	return phet_sols, grate_sols, ttrate_sols
end

# utils functions zone 
# focus mostly on formatting for plotting
