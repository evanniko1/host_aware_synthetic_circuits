#### TO-DO: ModelingToolkit.jl implementation
# implementation
# solve and check steady state solution
# explore optimization with Jacobian
# explore connected systems
#   i.  simplify growth model into blocks
#   ii. if i. == true, then examine host-circuit extentions
#   problems :: how to update shared resources ODEs ??


# import packages
using DifferentialEquations
using Plots

include("test_code.jl")

prob = ODEProblem(base_model_ODE!, base_model_ss_values, (0, 1e4), base_model_parameters)
sol = solve(prob, Rodas4())
plot(sol)

#### KEEP THE FOLLOWING FOR A MODELINGTOOLKIT.JL IMPLEMENTATION

base_model_parameters_a = Dict(
    thetar => 426.8693338968694,
    k_cm => 0.005990373118888,
    s0 => 10000,
    gmax => 1260.0,
    thetax => 4.379733394834643,
    Kt => 1.0e3,
    M => 1.0e8,
    we => 4.139172187824451,
    Km => 1000,
    vm => 5800.0,
    nx => 300.0,
    Kq => 1.522190403737490e+05,
    Kp => 180.1378030928276,
    vt => 726.0,
    wr => 929.9678874564831,
    wq => 948.9349882947897,
    wp => 0.0,
    nq => 4,
    nr => 7549.0,
    dm => 0.1,
    kb => 0.0095,
    ku => 1,
    ns => 0.5
    );

base_model_ss_values_a = [
    rmr => 807.530517658162,
	em => 7066.62814403594,
	rmq => 2316.16746788774,
	rmt => 69.1538861165317,
	et => 7066.62758622631,
	rmm => 69.1538892849256,
	mt => 9.64096853393699,
	mm => 9.64096853393699,
	q => 236681.293193619,
	si => 128.404551112062,
	mq => 322.904581569518,
	mr => 5.94267303259607,
	r => 17.3282796528522,
	a => 9.17329183561663
];