function base_model()
    @parameters thetar k_cm s0 gmax thetax Kt M we Km vm nx Kq Kp vt wr wq wp nq nr dm kb ku ns
    @variables t rmr(t) em(t) rmq(t) rmt(t) et(t) rmm(t) mt(t) mm(t) q(t) si(t) mq(t) mr(t) r(t) a(t) 
    D = Differential(t)

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
        r => 10
        a => 1e3
    ];

    # variables of type Num
    Kgamma= gmax/Kp
    gamma= gmax*a/(Kgamma + a)
    ttrate= (rmq + rmr + rmt + rmm)*gamma
    lam= ttrate/M 
    nucat= em*vm*si/(Km + si)

    # vector with 14 ODEs -- base model
    eqs = [
        D(rmr) ~ kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr,
        D(em)  ~ gamma/nx*rmm-lam*em,
        D(rmq) ~ kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq,
        D(rmt) ~ kb*r*mt-ku*rmt-gamma/nx*rmt-lam*rmt,
        D(et)  ~ gamma/nx*rmt-lam*et,
        D(rmm) ~ kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm,
        D(mt)  ~ (we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt,
        D(mm)  ~ (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm,
        D(q)   ~ gamma/nx*rmq-lam*q,
        D(si)  ~ (et*vt*s0/(Kt + s0))-nucat-lam*si,
        D(mq)  ~ (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq,
        D(mr)  ~ (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr,
        D(r)   ~ ku*rmr+ku*rmt+ku*rmm+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-lam*r,
        D(a)   ~ ns*nucat-ttrate-lam*a
    ];
    return eqs, base_model_parameters, base_model_ss_values, ttrate
end;

function heter_model()
    @parameters thetar k_cm s0 gmax thetax Kt M we Km vm nx Kq Kp vt wr wq wp nq nr dm kb ku n sw_max dp kb_h ku_h
    @variables t rmr(t) em(t) rmq(t) rmt(t) et(t) rmm(t) mt(t) mm(t) q(t) si(t) mq(t) mr(t) r(t) a(t) p_h(t) m_h(t) c_h(t)
    D = Differential(t)

    base_model_eqs, _, _, ttrate_base = base_model()
    Kgamma= gmax/Kp
    gamma= gmax*a/(Kgamma + a)
    ttrate = ttrate_base + c_h*gamma
    lam = ttrate / M
    nucat= em*vm*si/(Km + si)

    eqs_het = [
        D(p_h) ~ gamma/nx*c_h - (lam + dp)*p_h
        D(m_h) ~ wmaxhet*a/(thetax+a) - (lam + dm)*m_h + gamma/nx*c_h - kb_h*r*m_h + ku_h*c_h
        D(c_h) ~ -lam*c_h + kb_h*r*m_h - ku_h*c_h - gamma/nx*c_h
        D(r)   ~ ku*rmr+ku*rmt+ku*rmm+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-lam*r + ku_h*c_h - kb_h*r*m_h + gamma/nx*c_h
        D(a)   ~ ns*nucat-ttrate-lam*a
    ];

    host_aware_eqs = base_model_eqs[1:end-2]
    append!(host_aware_eqs, eqs_het)

    return host_aware_eqs
end