using ModelingToolkit, DifferentialEquations

# global base model parameters and variables
@parameters thetar k_cm s0 gmax thetax Kt M we Km vm nx Kq Kp vt wr wq wp nq nr dm kb ku ns
@variables t rmr(t) em(t) rmq(t) rmt(t) et(t) rmm(t) mt(t) mm(t) q(t) si(t) mq(t) mr(t) r(t) a(t) 

#some global calculations
Kgamma = gmax/Kp
gamma  = gmax*a/(Kgamma + a)
ttrate = (rmq + rmr + rmt + rmm)*gamma
lam    = ttrate/M 
nucat  = em*vm*si/(Km + si)

function base_model()

    @variables t
    D = Differential(t)

    eqs_hk_pr = [
        D(mq)  ~ (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq,
        D(rmq) ~ kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq,
        D(q)   ~ gamma/nx*rmq-lam*q
    ];
    eqs_m_pr = [
        D(mm)  ~ (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm,
        D(rmm) ~ kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm.
        D(em)  ~ gamma/nx*rmm-lam*em,
    ];
    eqs_t_pr = [
        D(mt)  ~ (we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt,
        D(rmt) ~ kb*r*mt-ku*rmt-gamma/nx*rmt-lam*rmt,
        D(et)  ~ gamma/nx*rmt-lam*et
    ];
    eqs_r_pr = [
        D(mr)  ~ (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr,
        D(rmr) ~ kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr,
        D(r)   ~ ku*rmr+ku*rmt+ku*rmm+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-lam*r
    ];
    eqs_nutrient_energy = [
        D(si)  ~ (et*vt*s0/(Kt + s0))-nucat-lam*si,
        D(a)   ~ ns*nucat-ttrate-lam*a
    ];

    return eqs
end;

