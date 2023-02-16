# import packages
using ModelingToolkit
using DifferentialEquations
using Plots
plotly()

# variables of type Num
Kgamma= gmax/Kp
gamma= gmax*a/(Kgamma + a)
ttrate= (rmq + rmr + rmt + rmm)*gamma
lam= ttrate/M 
nucat= em*vm*si/(Km + si)


function gene(; name)
    @parameters t nx w_max theta gmax Kp kb ku dm
    @variables p(t) m(t) c(t) a(t) r(t)
    Kgamma = gmax/Kp
    gamma = gmax*a / (Kgamma + a)
    D = Differential(t)
    eqs = [
        D(p) ~ gamma/nx*c
        D(m) ~ w_max*a/(theta+a) + gamma/nx*c - kb*r*m + ku*c - dm*m
        D(c) ~ kb*r*m - ku*c - gamma/nx*c
    ];
    ODESystem(eqs; name = name)
end;

function alpha(; name)
    @parameters t vm Km vt Kt s ns
    @variables em(t) et(t) si(t) a(t)
    nucat= em*vm*si/(Km + si)
    D = Differential(t)
    eqs = [
        D(si) ~ (et*vt*s/(Kt + s))-nucat
        D(a)  ~ ns*nucat-ttrate
    ]
    ODESystem(eqs; name = name)
end

#########
@parameters t
D = Differential(t)

@parameters gmax Kp M
Kgamma  = gmax/Kp
@variables rmq(t) rmr(t) rmt(t) rmm(t)
ttrate= (rmq + rmr + rmt + rmm)*gamma
lam= ttrate/M 
gamma= gmax*a/(Kgamma + a)


@named seqn = ODESystem([D(S) ~ -β * S * I / N])
@named ieqn = ODESystem([D(I) ~ β * S * I / N - γ * I])
@named reqn = ODESystem([D(R) ~ γ * I])

sir = compose(ODESystem([
                            S ~ ieqn.S,
                            I ~ seqn.I,
                            R ~ ieqn.R,
                            ieqn.S ~ seqn.S,
                            seqn.I ~ ieqn.I,
                            seqn.R ~ reqn.R,
                            ieqn.R ~ reqn.R,
                            reqn.I ~ ieqn.I], t, [S, I, R], [β, γ],
                            
                        defaults = [seqn.β => β
                                    ieqn.β => β
                                    ieqn.γ => γ
                                    reqn.γ => γ], name = :sir), seqn, ieqn, reqn)


#################
include("test_code.jl")

base_eqs, base_params, base_vars, trate = base_model();

# method to map base_params and base_vars to dicts of parameter and initial values respectively


base_model_parameters = [
    thetar => 426.8693338968694
    k_cm => 0.005990373118888
    s0 => 10000
    gmax => 1260.0
    thetax => 4.379733394834643
    Kt => 1.0e3
    M => 1.0e8
    we => 4.139172187824451
    Km => 1000
    vm => 5800.0
    nx => 300.0
    Kq => 1.522190403737490e+05
    Kp => 180.1378030928276
    vt => 726.0
    wr => 929.9678874564831
    wq => 948.9349882947897
    wp => 0.0
    nq => 4
    nr => 7549.0
    dm => 0.1
    kb => 0.0095
    ku => 1
    ns => 0.5
    ];

base_model_ss_values = [
    rmr => 0
    em => 0
    rmq => 0
    rmt => 0
    et => 0
    rmm => 0
    mt => 0
    mm => 0
    q => 0
    si => 0
    mq => 0
    mr => 0
    r => 0
    a => 1e3
];

@named base_model_sys = ODESystem(base_eqs);
prob = ODEProblem(base_model_sys, base_model_ss_values, (0, 1e4), base_model_parameters; jac=true);
sol  = solve(prob, Rodas4());
plot(sol)

ha_eqs, _, _ = heter_model();

@variables t base_vars[1]

##################
include("GrowthModels.jl")
using Main.GrowthModels

base = Main.GrowthModels.base_model();
het = Main.GrowthModels.heter_model(; extend_from = base);

@named base_model = ODESystem(base)


dadteq = equations(base_model)[end]           
t      = independent_variable(base_model)   
@variables t a(t) c_h(t) 
@parameters M gmax Kp
Kgamma = gmax/Kp
gamma = gmax*a / (a + Kgamma)
dadteq = Equation(dadteq.lhs, dadteq.rhs - c_h*gamma/M)   
@named test_ext  = ODESystem([dadteq], t, states(base_model), parameters(base_model))
#oprob  = ODEProblem(test_ext, u0map, tspan, pmap)
#osol   = solve(oprob, Tsit5()