module GrowthModels
    
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
        return eqs
    end;
    
    function heter_model(; extend_from)
        @parameters w_max dp kb_h ku_h
        @variables t p_h(t) m_h(t) c_h(t)
        D = Differential(t)
    
        # get equations of base model
        base_model_eqs = extend_from
        
        # add contribution to growth
        ttrate_het = ttrate + c_h*gamma
        lam = ttrate_het / M
    
        eqs_het = [
            D(p_h) ~ gamma/nx*c_h - (lam + dp)*p_h
            D(m_h) ~ w_max*a/(thetax+a) - (lam + dm)*m_h + gamma/nx*c_h - kb_h*r*m_h + ku_h*c_h
            D(c_h) ~ -lam*c_h + kb_h*r*m_h - ku_h*c_h - gamma/nx*c_h
            D(r)   ~ ku*rmr+ku*rmt+ku*rmm+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-lam*r + ku_h*c_h - kb_h*r*m_h + gamma/nx*c_h
            D(a)   ~ ns*nucat-ttrate_het-lam*a
        ];
    
        #host_aware_eqs = base_model_eqs[1:end-2]
        #append!(host_aware_eqs, eqs_het)
        #host_aware_eqs = append!(base_model_eqs, eqs_het)
        #return host_aware_eqs
        
        return eqs_het
    end

end