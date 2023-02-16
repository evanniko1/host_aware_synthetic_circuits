using ModelingToolkit, DifferentialEquations

# global base model parameters and variables
@parameters thetar k_cm s0 gmax thetax Kt M we Km vm nx Kq Kp vt wr wq wp nq nr dm kb ku ns  
@variables t mq(t) rmq(t) q(t) mm(t) rmm(t) em(t) mt(t) rmt(t) et(t) mr(t) rmr(t) r(t) si(t) a(t)
#some global calculations
Kgamma = gmax/Kp
gamma  = gmax*a/(Kgamma + a)
nucat  = em*vm*si/(Km + si)

function base_model_eqs()
    lam = ((rmq + rmr + rmt + rmm)*gamma)/M
    #ttrate = (rmq + rmr + rmt + rmm)*gamma
    #lam    = ttrate/M 

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
    eqs_r_mc = [
        D(mr)  ~ (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr,
        D(rmr) ~ kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr,
    ];
    eqs_r_pr = [
        D(r)   ~ ku*rmr+ku*rmt+ku*rmm+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-lam*r
    ];          
    eqs_s = [
        D(si)  ~ (et*vt*s0/(Kt + s0))-nucat-lam*si,
    ];
    eqs_a = [
        D(a)   ~ ns*nucat-lam*M-lam*a
    ];

    # organize all protein fractions into a nice dictionary
    eqs_dict = Dict(
        "house_keeping" => eqs_hk_pr,
        "metabolic"     => eqs_m_pr,
        "transporter"   => eqs_t_pr,
        "ribo_mc"       => eqs_r_mc,
        "ribosomes"     => eqs_r_pr,
        "nutrient"      => eqs_s,
        "energy"        => eqs_a
    )
    return eqs_dict, lam
end

function compose_base_model(; base_eqs_dict)
    D= Differential(t)
    connected = compose(ODESystem([
    base_eqs_dict["house_keeping"]
    base_eqs_dict["metabolic"]
    base_eqs_dict["transporter"]
    base_eqs_dict["ribo_mc"]
    base_eqs_dict["ribosomes"]
    base_eqs_dict["nutrient"]
    base_eqs_dict["energy"]
    ], t; name = :connected));

    return connected
end

function het_model_eqs(; base_eqs_dict, base_lam)
    @parameters w_max dp kb_h ku_h
    @variables t m_h(t) c_h(t) p_h(t) r(t) a(t)
    D = Differential(t)
 
    # add contribution to growth
    lam = ((rmq + rmr + rmt + rmm )*gamma)/M

    eqs_het = [
        D(m_h) ~ w_max*a/(thetax+a) - (lam + dm)*m_h + gamma/nx*c_h - kb_h*r*m_h + ku_h*c_h
        D(c_h) ~ -lam*c_h + kb_h*r*m_h - ku_h*c_h - gamma/nx*c_h
        D(p_h) ~ gamma/nx*c_h - (lam + dp)*p_h
    ];
    # update energy and ribosomal contributions
    eqs_r_pr_ha = [
        D(r) ~ base_eqs_dict["ribosomes"][1].rhs - kb_h*m_h*r + ku_h*c_h + gamma/nx*c_h
    ];
    eqs_a_ha = [
        D(a) ~ base_eqs_dict["energy"][1].rhs - gamma*c_h
    ];

    # organize host_aware model equations in a dictionary
    eqs_dict = Dict(
        "house_keeping" => base_eqs_dict["house_keeping"],
        "metabolic"     => base_eqs_dict["metabolic"],
        "transporter"   => base_eqs_dict["transporter"],
        "ribo_mc"       => base_eqs_dict["ribo_mc"],
        "ribosomes"     => eqs_r_pr_ha,
        "heterologous"  => eqs_het,
        "nutrient"      => base_eqs_dict["nutrient"],
        "energy"        => eqs_a_ha
    )

    return eqs_dict, ttrate_het
end


base_eqs_dict, base_ttrate = base_model_eqs();
het_eqs, _ = het_model_eqs(base_eqs_dict = base_eqs_dict, base_ttrate = base_ttrate);

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


# try out extend() -- WORKS 
nn = extend(het_iso, connected);

# TODO -- UPDATE SHARED RESOURCES EQUATIONS!!!