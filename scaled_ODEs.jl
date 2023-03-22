using DifferentialEquations, Plots
plotly()

function scaled_ODEmodel!(du,u,p,t)
    # parameters

    dm= p[1];
    thetar= p[2];
    we= p[3];
    vt= p[4];
    s0= p[5];
    nx= p[6];
    nq= p[7];
    nr= p[8];
    ns= p[9];
    gmax= p[10];
    thetax= p[11];
    Km= p[12];
    Kq= p[13];
    Kp= p[14];
    wr= p[15];
    wq= p[16];
    Kt= p[17];
    kb= p[18];
    Mref= p[19];
    ku= p[20];
    vm= p[21];

    # initial conditions
    rmr= u[1];
    em= u[2];
    rmq= u[3];
    rmt= u[4];
    et= u[5];
    rmm= u[6];
    a= u[7];
    mm= u[8];
    q= u[9];
    mt= u[10];
    si= u[11];
    r= u[12];
    mq= u[13];
    mr= u[14];

    # algebraic relationships

    M= (nr*(r+rmr+rmt+rmm+rmq)+nx*(q+et+em));
    Kgamma= gmax/Kp*M/Mref;
    gamma= gmax*a/(Kgamma+a);
    ttrate= (rmq+rmr+rmt+rmm)*gamma;
    lam= ttrate/M;
    nucat= em*vm*si/(Km*M/Mref+si);
    
    # ODEs
    du[1]= +kb*Mref/M*r*mr+ku*rmr-gamma/nr*rmr;
    du[2]= +gamma/nx*rmm;
    du[3]= +kb*Mref/M*r*mq-ku*rmq-gamma/nx*rmq;
    du[4]= +kb*Mref/M*r*mt-ku*rmt-gamma/nx*rmt;
    du[5]= +gamma/nx*rmt;
    du[6]= +kb*Mref/M*r*mm-ku*rmm-gamma/nx*rmm;
    du[7]= +ns*nucat-ttrate;
    du[8]= +(M/Mref*we*a/(thetax*M/Mref+a))+ku*rmm+gamma/nx*rmm-kb*Mref/M*r*mm-dm*mm;
    du[9]= +gamma/nx*rmq;
    du[10]= +(M/Mref*we*a/(thetax*M/Mref+a))+ku*rmt+gamma/nx*rmt-kb*Mref/M*r*mt-dm*mt;
    du[11]= +(et*vt*s0/(Kt+s0))-nucat;
    du[12]= +ku*rmr+ku*rmt+ku*rmm+ku*rmq+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+-kb*Mref/M*r*mr-kb*Mref/M*r*mt-kb*Mref/M*r*mm-kb*Mref/M*r*mq;
    du[13]= +(M/Mref*wq*a/(thetax*M/Mref+a)/(1+(q/(Kq*M/Mref))^nq))+ku*rmq+gamma/nx*rmq-kb*Mref/M*r*mq-dm*mq;
    du[14]= +(M/Mref*wr*a/(thetar*M/Mref+a))+ku*rmr+gamma/nr*rmr-kb*Mref/M*r*mr-dm*mr;
 
end

dm= 0.1;
zeta = 2.8312e+06;

thetar= 426.8693338968694;
we= 4.139172187824451;
vt= 726.0;
s0= 1.0e4;
nx= 300.0;
nq= 4;
nr= 7549.0;
ns= 0.5;
gmax= 1260.0;
thetax= 4.379733394834643;
Km= 1.0e3;
Kq= 1.522190403737490e+05;
Kp= 180.1378030928276;
wr= 929.9678874564831;
wq= 948.9349882947897;
Kt= 1.0e3;
kb= 1.0;
Mref= 1.0e8;
ku= 1.0;
vm= 5800.0;

thetar = thetar*zeta;
thetax = thetax * zeta;
Kp = Kp/zeta;

# parameters vector
p = [dm, thetar, we, vt, s0, nx, nq, nr, ns, gmax, thetax, Km, Kq, Kp, wr, wq, Kt, kb, Mref, ku, vm]

#initial conditions vector
u0 = [1, 1, 1, 1, 1, 1, 1000, 1, 1, 1, 1, 10, 1, 1]

######
tspan = (0.0,1e5) 
prob = ODEProblem(scaled_ODEmodel!,u0,tspan,p)
sol = solve(prob, Rodas4(), callback=PositiveDomain())
plot(sol)