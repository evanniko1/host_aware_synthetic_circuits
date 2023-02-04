# import packages
using DifferentialEquations
using Plots

# import custom stuff
include("helper.jl")

##########
# exposed relevant parameters
ns = 0.5
kurep= exp10(−2.6575) #default value is *1*
kbrep= exp10(−1.3335) #default value is *0.0095*
wmaxrep = 150
kappa_ini = 0.1

# assemble in a parameter vector
p= [thetar, k_cm, s0, gmax, cl, 
	thetax, Kt, M, we, Km, 
	vm, nx, Kq, Kp, vt, 
	wr, wq, wp, nq, nr, 
	ns, wmaxrep, kbrep, kurep, dmrep, dprep, kappa_ini]

##########
# one run of the model
# time span
tspan = (0.0, 1e5)

# specify ODE problem
prob = ODEProblem(ODE_model!,u0,tspan,p)

# solve ODE problem
sol = solve(prob,Rodas4())
plotly();
plot(sol) # time trajectories for all species; the plot is interactive

##########
# reproduce the non-linear relationship from Cambray et al, 2018
phet, grate = trans_initiation!(init_values = u0, tspan = (0,1e5), params_values = p);
plot(phet, grate)