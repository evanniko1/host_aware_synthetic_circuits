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