using ModelingToolkit, DifferentialEquations

# global base model parameters and variables
@parameters thetar k_cm s0 gmax thetax Kt M we Km vm nx Kq Kp vt wr wq wp nq nr dm kb ku ns  
@variables t mq(t) rmq(t) q(t) mm(t) rmm(t) em(t) mt(t) rmt(t) et(t) mr(t) rmr(t) r(t) si(t) a(t)
#some global calculations
Kgamma = gmax/Kp
gamma  = gmax*a/(Kgamma + a)
ttrate = (rmq + rmr + rmt + rmm)*gamma
lam    = ttrate/M 
nucat  = em*vm*si/(Km + si)


function base_model_eqs()
    @variables t
    D = Differential(t)

    eqs_hk_pr = [
        D(mq)  ~ (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq,
        D(rmq) ~ kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq,
        D(q)   ~ gamma/nx*rmq-lam*q
    ];
    eqs_m_pr = [
        D(mm)  ~ (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm,
        D(rmm) ~ kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm,
        D(em)  ~ gamma/nx*rmm-lam*em
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

    eqs_dict = Dict(
        "house_keeping" => eqs_hk_pr,
        "metabolic" => eqs_m_pr,
        "transporter" => eqs_t_pr,
        "ribosomes" => eqs_r_pr,
        "nutr_energy" => eqs_nutrient_energy
    )
    return eqs_dict
end

function het_model_eqs()
    @parameters w_max dp kb_h ku_h
    @variables t m_h(t) c_h(t) p_h(t)
    D = Differential(t)
 
    # add contribution to growth
    ttrate_het = ttrate + c_h*gamma
    lam = ttrate_het / M

    eqs_het = [
        D(m_h) ~ w_max*a/(thetax+a) - (lam + dm)*m_h + gamma/nx*c_h - kb_h*r*m_h + ku_h*c_h
        D(c_h) ~ -lam*c_h + kb_h*r*m_h - ku_h*c_h - gamma/nx*c_h
        D(p_h) ~ gamma/nx*c_h - (lam + dp)*p_h
    ];
end


base_eqs_dict = base_model_eqs();
het_eqs = het_model_eqs();
#@named house_keeping = ODESystem(base_eqs_dict["house_keeping"])
#@named metabolic = ODESystem(base_eqs_dict["metabolic"])
#@named transporter = ODESystem(base_eqs_dict["transporter"])
#@named ribosomes = ODESystem(base_eqs_dict["ribosomes"])
#@named nutr_energy = ODESystem(base_eqs_dict["nutr_energy"])

# compose the base model from proteome fractions
D = Differential(t)
connected = compose(ODESystem([
    base_eqs_dict["house_keeping"]
    base_eqs_dict["metabolic"]
    base_eqs_dict["transporter"]
    base_eqs_dict["ribosomes"]
    base_eqs_dict["nutr_energy"]
    ], t; name = :connected));

het_iso = compose(ODESystem(het_eqs, t; name = :het_iso));
    
#equations(connected)
@variables t m_h(t) c_h(t)
@parameters ku_h kb_h
# extend base model with heterologous protein
connections = [
    connected.phldr1 ~ gamma*c_h
    connected.phldr2 ~ ku_h*c_h - kb_h*r*m_h + gamma/nx*c_h
    ]

D_h = Differential(t)
ha_connected = compose(
    ODESystem(
        connections, t; name = :ha_connected
    ), het_iso, connected);

simplified_ha = structural_simplify(ha_connected);

# try out extend()
nn = extend(het_iso, connected);