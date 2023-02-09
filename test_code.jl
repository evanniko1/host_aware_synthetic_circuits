# TEST FILE FOR DICTIONARY BASED MODEL DEFINITIONS
function ODE_model!(du,u,p,t)
    wmaxrep = p[22]
    kbrep= p[23]
    kurep= p[24]
    dmrep = p[25]
    dprep = p[26]
    kappa_ini = p[27]
    
	rmr= u[1] # complex between a ribosome and the mRNA for protein r
	em= u[2] # enzyme that metabolizes si (nutrinet) into a (energy)
	rmq= u[3] # complex between a ribosome and the mRNA for protein q
	rmt= u[4] # complex between a ribosome and the mRNA for protein et
	et= u[5] # enzyme that transports s (energy) into the cell
	rmm= u[6] # complex between a ribosome and the mRNA for protein em
	mt= u[7] # mRNA of et
	mm= u[8] # mRNA of em
	q= u[9] # house keeping proteins
	si= u[10] # nutrient (internalized)
	mq= u[11] # mRNA of q
	mr= u[12] # mRNA of r (rybosomes)
	r= u[13] # number of rybosomes
	a= u[14] # energy


	Kgamma= gmax/Kp # threshold amount of energy where elongation is halfmaximal
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation
	ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep)*gamma # term in expression of lambda
	lam= ttrate/M # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of si by e
 
	du[1]= kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr # (6)
	du[2]= gamma/nx*rmm-lam*em # (4)
	du[3]= kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq # (6)
	du[4]= kb*r*mt+b-ku*rmt-gamma/nx*rmt-lam*rmt # (6)
	du[5]= gamma/nx*rmt-lam*et # (4)
	du[6]= kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm # (6)
	du[7]= (we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt # (5)
	du[8]= (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm # (5)
	du[9]= gamma/nx*rmq-lam*q # (4)
	du[10]= (et*vt*s0/(Kt + s0))-nucat-lam*si # (1)  import and metabolism of si
	du[11]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq # (5)
	du[12]= (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr # (5)
	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep*rmrep+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep*r*mrep-lam*r # (3) gamma/nr*rmr translation rate of r i.e. mur
	du[14]= ns*nucat-ttrate-lam*a # (2) energy is created from metabolism and lost from translation and dilution
end

base_model_parameters = Dict(
    "thetar" => 426.8693338968694,
    "k_cm" => 0.005990373118888,
    "s0" => 10000,
    "gmax" => 1260.0,
    "thetax" => 4.379733394834643,
    "Kt" => 1.0e3,
    "M" => 1.0e8,
    "we" => 4.139172187824451,
    "Km" => 1000,
    "vm" => 5800.0,
    "nx" => 300.0,
    "Kq" => 1.522190403737490e+05,
    "Kp" => 180.1378030928276,
    "vt" => 726.0,
    "wr" => 929.9678874564831,
    "wq" => 948.9349882947897,
    "wp" => 0.0,
    "nq" => 4,
    "nr" => 7549.0,
    "dm" => 0.1,
    "kb" => 0.0095,
    "ku" => 1,
    "ns" => 0.5
    );

base_model_ss_values = Dict(
    "rmr" => 807.530517658162,
	"em" => 7066.62814403594,
	"rmq" => 2316.16746788774,
	"rmt" => 69.1538861165317,
	"et" => 7066.62758622631,
	"rmm" => 69.1538892849256,
	"mt" => 9.64096853393699,
	"mm" => 9.64096853393699,
	"q" => 236681.293193619,
	"si" => 128.404551112062,
	"mq" => 322.904581569518,
	"mr" => 5.94267303259607,
	"r" => 17.3282796528522,
	"a" => 9.17329183561663
);

function base_model_ODE!(base_model_parameters, u, du, t)

	rmr= u[1]
	em= u[2]
	rmq= u[3]
	rmt= u[4]
	et= u[5]
	rmm= u[6]
	mt= u[7]
	mm= u[8]
	q= u[9]
	si= u[10]
	mq= u[11]
	mr= u[12]
	r= u[13]
	a= u[14]

    # helpful definitions
    Kgamma= base_model_parameters["gmax"]/base_model_parameters["Kp"]
    gamma= base_model_parameters["gmax"]*a/(Kgamma + a)
    ttrate= (rmq + rmr + rmt + rmm)*gamma
    lam= ttrate/base_model_parameters["M"]
    nucat= em*vm*si/(base_model_parameters["Km"] + si)

    # species ODEs
    du[1]= base_model_parameters["kb"]*r*mr-base_model_parameters["ku"]*rmr-gamma/base_model_parameters["nr"]*rmr-lam*rmr
    du[2]= gamma/base_model_parameters["nx"]*rmm-lam*em
    du[3]= base_model_parameters["kb"]*r*mq-base_model_parameters["ku"]*rmq-gamma/base_model_parameters["nx"]*rmq-lam*rmq
    du[4]= base_model_parameters["kb"]*r*mt-base_model_parameters["ku"]*rmt-gamma/base_model_parameters["nx"]*rmt-lam*rmt
    du[5]= gamma/base_model_parameters["nx"]*rmt-lam*et
    du[6]= base_model_parameters["kb"]*r*mm-base_model_parameters["ku"]*rmm-gamma/base_model_parameters["nx"]*rmm-lam*rmm
    du[7]= (base_model_parameters["we"]*a/(base_model_parameters["thetax"] + a))+base_model_parameters["ku"]*rmt+gamma/base_model_parameters["nx"]*rmt-base_model_parameters["kb"]*r*mt-base_model_parameters["dm"]*mt-lam*mt
    du[8]= (base_model_parameters["we"]*a/(base_model_parameters["thetax"] + a))+base_model_parameters["ku"]*rmm+gamma/base_model_parameters["nx"]*rmm-base_model_parameters["kb"]*r*mm-base_model_parameters["dm"]*mm-lam*mm
    du[9]= gamma/base_model_parameters["nx"]*rmq-lam*q
    du[10]= (et*base_model_parameters["vt"]*base_model_parameters["s0"]/(base_model_parameters["Kt"] + base_model_parameters["s0"]))-nucat-lam*si
    du[11]= (base_model_parameters["wq"]*a/(base_model_parameters["thetax"] + a)/(1 + (q/base_model_parameters["Kq"])^base_model_parameters["nq"]))+base_model_parameters["ku"]*rmq+gamma/base_model_parameters["nx"]*rmq-base_model_parameters["kb"]*r*mq-base_model_parameters["dm"]*mq-lam*mq
    du[12]= (base_model_parameters["wr"]*a/(base_model_parameters["thetar"] + a))+base_model_parameters["ku"]*rmr+gamma/base_model_parameters["nr"]*rmr-base_model_parameters["kb"]*r*mr-base_model_parameters["dm"]*mr-lam*mr
    du[13]= base_model_parameters["ku"]*rmr+base_model_parameters["ku"]*rmt+base_model_parameters["ku"]*rmm+base_model_parameters["ku"]*rmq+gamma/nr*rmr+gamma/base_model_parameters["nr"]*rmr+gamma/base_model_parameters["nx"]*rmt+gamma/base_model_parameters["nx"]*rmm+gamma/base_model_parameters["nx"]*rmq-base_model_parameters["kb"]*r*mr-base_model_parameters["kb"]*r*mt-base_model_parameters["kb"]*r*mm-base_model_parameters["kb"]*r*mq-lam*r
    du[14]= base_model_parameters["ns"]*nucat-ttrate-lam*a
end