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
    @parameters t nx w_max theta gmax Kp kb ku
    @variables p(t) m(t) c(t) a(t) r(t)
    Kgamma = gmax/Kp
    gamma = gmax*a / (Kgamma + a)
    D = Differential(t)
    eqs = [
        D(p)  ~ gamma/nx*c
        D(m)  ~ w_max*a/(theta+a) + gamma/nx*c - kb*r*m + ku*c
        D(c)  ~ kb*r*m - ku*c - gamma/nx*c
    ];
    ODESystem(eqs; name = name)
end;



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