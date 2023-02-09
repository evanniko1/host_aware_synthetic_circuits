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
@parameters t thetar k_cm s0 gmax thetax Kt M we Km vm nx Kq Kp vt wr wq wp nq nr dm kb ku ns
@variables rmr(t) em(t) rmq(t) rmt(t) et(t) rmm(t) mt(t) mm(t) q(t) si(t) mq(t) mr(t) r(t) a(t) 
D = Differential(t)

# variables of type Num
Kgamma= gmax/Kp
gamma= gmax*a/(Kgamma + a)
ttrate= (rmq + rmr + rmt + rmm)*gamma
lam= ttrate/M 
nucat= em*vm*si/(Km + si)

eqs = [
    D(rmr) ~ kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr
    D(em) ~ gamma/nx*rmm-lam*em
    D(rmq) ~ kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq
    D(rmt) ~ kb*r*mt-ku*rmt-gamma/nx*rmt-lam*rmt
    D(et) ~ gamma/nx*rmt-lam*et
    D(rmm) ~ kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm
    D(mt) ~ (we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt
    D(mm) ~ (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm
    D(q) ~ gamma/nx*rmq-lam*q
    D(si) ~ (et*vt*s0/(Kt + s0))-nucat-lam*si
    D(mq) ~ (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq
    D(mr) ~ (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr
    D(r) ~ ku*rmr+ku*rmt+ku*rmm+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-lam*r
    D(a) ~ ns*nucat-ttrate-lam*a
]

@named base_model_sys = ODESystem(eqs)

tspan = (0.0,10000.0)
prob = ODEProblem(base_model_sys, base_model_ss_values, tspan, base_model_parameters;jac=true,sparse=true)
sol = solve(prob, Rodas4())
plot(sol)