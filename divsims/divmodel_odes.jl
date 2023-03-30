# divmodel_odes() is currently open -- NEEDS END

function divmodel_odes(t, y, data)
    rates = data.rates
    parameters = data.p
    
    b = rates[1]
    dm = rates[2]
    
    thetar = parameters[1]
    cl = parameters[2]
    we = parameters[3]
    kby = parameters[4]
    vt = parameters[5]
    s0 = parameters[6]
    nx = parameters[7]
    nq = parameters[8]
    nr = parameters[9]
    ns = parameters[10]
    wg = parameters[11]
    k_cm = parameters[12]
    gmax = parameters[13]
    thetax = parameters[14]
    Km = parameters[15]
    kbg = parameters[16]
    Kq = parameters[17]
    Kp = parameters[18]
    wr = parameters[19]
    wq = parameters[20]
    Kt = parameters[21]
    wy = parameters[22]
    kb = parameters[23]
    Mref = parameters[24]
    ku = parameters[25]
    vm = parameters[26]
    C = parameters[27]
    D = parameters[28]
    Vcrit = parameters[29]
    
    rmr = y[1]
    em = y[2]
    rmq = y[3]
    rmt = y[4]
    gfp = y[5]
    rmy = y[6]
    rmg = y[7]
    et = y[8]
    rmm = y[9]
    zmm = y[10]
    zmg = y[11]
    zmy = y[12]
    zmr = y[13]
    zmq = y[14]
    zmt = y[15]
    yfp = y[16]
    a = y[17]
    mg = y[18]
    mm = y[19]
    q = y[20]
    mt = y[21]
    si = y[22]
    r = y[23]
    mq = y[24]
    mr = y[25]
    my = y[26]

    if length(y) > 26
        fpos = y[27:end]
    end

    M = nr*(r+rmg+rmy+rmr+rmt+rmm+rmq+zmg+zmy+zmr+zmt+zmm+zmq)+nx*(gfp+yfp+q+et+em)
    Kgamma = gmax/Kp*M/Mref
    gamma = gmax*a/(Kgamma+a)
    ttrate = (rmg+rmy+rmq+rmr+rmt+rmm)*gamma
    lam = ttrate/M
    fr = nr*(r+rmg+rmy+rmr+rmt+rmm+rmq+zmg+zmy+zmr+zmt+zmm+zmq)/M
    nucat = em*vm*si/(Km*M/Mref+si)
    f = cl*k_cm*Mref/M

    dydt = zeros(26)
    dydt[1] = +kb*Mref/M*r*mr+b*zmr-ku*rmr-gamma/nr*rmr-f*rmr
    dydt[2] = +gamma/nx*rmm
    dydt[3] = +kb*Mref/M*r*mq+b*zmq-ku*rmq-gamma/nx*rmq-f*rmq
    dydt[4]= +kb*Mref/M*r*mt+b*zmt-ku*rmt-gamma/nx*rmt-f*rmt;
    dydt[5]= +gamma/nx*rmg;
    dydt[6]= +kby*Mref/M*r*my+b*zmy-ku*rmy-gamma/nx*rmy-f*rmy;
    dydt[7]= +kbg*Mref/M*r*mg+b*zmg-ku*rmg-gamma/nx*rmg-f*rmg;
    dydt[8]= +gamma/nx*rmt;
    dydt[9]= +kb*Mref/M*r*mm+b*zmm-ku*rmm-gamma/nx*rmm-f*rmm;
    dydt[10]= +f*rmm-b*zmm;
    dydt[11]= +f*rmg-b*zmg;
    dydt[12]= +f*rmy-b*zmy;
    dydt[13]= +f*rmr-b*zmr;
    dydt[14]= +f*rmq-b*zmq;
    dydt[15]= +f*rmt-b*zmt;
    dydt[16]= +gamma/nx*rmy;
    dydt[17]= +ns*nucat-ttrate;
    dydt[18]= +(M/Mref*wg*a/(thetax*M/Mref+a))+ku*rmg+gamma/nx*rmg-kbg*Mref/M*r*mg-dm*mg;
    dydt[19]= +(M/Mref*we*a/(thetax*M/Mref+a))+ku*rmm+gamma/nx*rmm-kb*Mref/M*r*mm-dm*mm;
    dydt[20]= +gamma/nx*rmq;
    dydt[21]= +(M/Mref*we*a/(thetax*M/Mref+a))+ku*rmt+gamma/nx*rmt-kb*Mref/M*r*mt-dm*mt;
    dydt[22]= +(et*vt*s0/(Kt+s0))-nucat;
    dydt[23]= +ku*rmr+ku*rmt+ku*rmm+ku*rmq+ku*rmg+ku*rmy+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+gamma/nx*rmg+gamma/nx*rmy-kb*Mref/M*r*mr-kb*Mref/M*r*mt-kb*Mref/M*r*mm-kb*Mref/M*r*mq-kbg*Mref/M*r*mg-kby*Mref/M*r*my;
    dydt[24]= +(M/Mref*wq*a/(thetax*M/Mref+a)/(1+(q/(Kq*M/Mref))^nq))+ku*rmq+gamma/nx*rmq-kb*Mref/M*r*mq-dm*mq;
    dydt[25]= +(M/Mref*wr*a/(thetar*M/Mref+a))+ku*rmr+gamma/nr*rmr-kb*Mref/M*r*mr-dm*mr;
    dydt[26]= +(M/Mref*wy*a/(thetax*M/Mref+a))+ku*rmy+gamma/nx*rmy-kby*Mref/M*r*my-dm*my;

if length(y) > 26
    dfposdt = ones(length(fpos));
    if haskey(data, "CD") && data.CD
        dtposdt = dfposdt * 1/(C+D) * a/(K_cd*M/Mref + a)
    end
    dydt = vcat(dydt, dfposdt)
end

flag = 0
new_data = []