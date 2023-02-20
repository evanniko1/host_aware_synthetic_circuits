module GrowthModels
    
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
    );

    ss_values = Dict([
        "rmr" => rmr => 0
        "em"  => em  => 0
        "rmq" => rmq => 0
        "rmt" => rmt => 0
        "et"  => et  => 0
        "rmm" => rmm => 0
        "mt"  => mt  => 0
        "mm"  => mm  => 0
        "q"   => q   => 0
        "si"  => si  => 0
        "mq"  => mq  => 0
        "mr"  => mr  => 0
        "r"   => r   => 10
        "a"   => a   => 1e3
    ]);

    param_values = Dict([
        "thetar" => thetar => 426.8693338968694
        "k_cm" => k_cm => 0.005990373118888
        "s0" => s0 => 10000
        "gmax" => gmax => 1260.0
        "thetax" => thetax => 4.379733394834643
        "Kt" => Kt => 1.0e3
        "M" => M => 1.0e8
        "we" => we => 4.139172187824451
        "Km" => Km => 1000
        "vm" => vm => 5800.0
        "nx" => nx => 300.0
        "Kq" => Kq => 1.522190403737490e+05
        "Kp" => Kp => 180.1378030928276
        "vt" => vt => 726.0
        "wr" => wr => 929.9678874564831
        "wq" => wq => 948.9349882947897
        "wp" => wp => 0.0
        "nq" => nq => 4
        "nr" => nr => 7549.0
        "dm" => dm => 0.1
        "kb" => kb => 0.0095
        "ku" => ku => 1
        "ns" => ns => 0.5
    ]);
    return eqs_dict, lam, param_values, ss_values
end;

function het_model_eqs(; input_eqs_dict::Dict, input_lam::Num, input_param_vals::Dict, input_ss_vals::Dict)
    @parameters w_max dp kb_h ku_h
    @variables t m_h(t) c_h(t) p_h(t) r(t) a(t)
    D = Differential(t)

    het_param_vals = Dict([
        "w_max" => w_max => 150
        "dp"   => dp    => 0.1
        "kb_h" => kb_h  => exp10(-1.3335)
        "ku_h" => ku_h  => exp10(-2.6575)
    ])
    het_ss_vals = Dict([
        "m_h" => m_h => 0.0
        "c_h" => c_h => 0.0
        "p_h" => p_h => 0.0
    ])

    # add contribution to growth
    lam = input_lam + c_h*gamma/M

    # update all growth rates for the input system of equations
    upd_eqs_dict = update_eqs!(eqs_dict= input_eqs_dict, term = input_lam - lam);

    # define equations for heterologous species
    eqs_het = [
        D(m_h) ~ w_max*a/(thetax+a) - (lam + dm)*m_h + gamma/nx*c_h - kb_h*r*m_h + ku_h*c_h
        D(c_h) ~ -lam*c_h + kb_h*r*m_h - ku_h*c_h - gamma/nx*c_h
        D(p_h) ~ gamma/nx*c_h - (lam + dp)*p_h
    ];
    # update energy and ribosomal contributions
    eqs_r_pr_ha = [
        D(r) ~ upd_eqs_dict["ribosomes"][1].rhs - kb_h*m_h*r + ku_h*c_h + gamma/nx*c_h
    ];
    eqs_a_ha = [
        D(a) ~ upd_eqs_dict["energy"][1].rhs - gamma*c_h
    ];

    # organize host_aware model equations in a dictionary
    eqs_dict = Dict(
        "house_keeping" => upd_eqs_dict["house_keeping"],
        "metabolic"     => upd_eqs_dict["metabolic"],
        "transporter"   => upd_eqs_dict["transporter"],
        "ribo_mc"       => upd_eqs_dict["ribo_mc"],
        "ribosomes"     => eqs_r_pr_ha,
        "heterologous"  => eqs_het,
        "nutrient"      => upd_eqs_dict["nutrient"],
        "energy"        => eqs_a_ha
    )

    param_vals = merge(input_param_vals, het_param_vals)
    ss_vals    = merge(input_ss_vals, het_ss_vals)

    return eqs_dict, lam, param_vals, ss_vals
end;

function update_eqs!(; eqs_dict::Dict, term::Num)
    # create a deep copy of the dictionary
    dict_to_updt = deepcopy(eqs_dict)

    for key in keys(eqs_dict)
        for idx in range(1, length(eqs_dict[key]))
            dict_to_updt[key][idx] = Equation(eqs_dict[key][idx].lhs, eqs_dict[key][idx].rhs + term)
        end
    end

    return dict_to_updt
end;

function compose_model(eqs_dict::Dict)
    D= Differential(t)
    connected = compose(ODESystem(

        [eqs_dict[key][idx] for key in keys(eqs_dict) for idx in range(1, length(eqs_dict[key]))], t; name = :connected));
    return connected
end;

# utility function for dictionary formating and updating
function _format(dict::Dict)
    return vec(collect(values(dict)))
end;

function _update!(dict::Dict, param::String, val::Float64)
    dict[param] = dict[param][1] => val
end;

#################
base_model_ss_values = Dict([
    "rmr" => rmr => 0
    "em"  => em  => 0
    "rmq" => rmq => 0
    "rmt" => rmt => 0
    "et"  => et  => 0
    "rmm" => rmm => 0
    "mt"  => mt  => 0
    "mm"  => mm  => 0
    "q"   => q   => 0
    "si"  => si  => 0
    "mq"  => mq  => 0
    "mr"  => mr  => 0
    "r"   => r   => 10
    "a"   => a   => 1e3
]);
base_model_parameters = Dict([
    "thetar" => thetar => 426.8693338968694
    "k_cm" => k_cm => 0.005990373118888
    "s0" => s0 => 10000
    "gmax" => gmax => 1260.0
    "thetax" => thetax => 4.379733394834643
    "Kt" => Kt => 1.0e3
    "M" => M => 1.0e8
    "we" => we => 4.139172187824451
    "Km" => Km => 1000
    "vm" => vm => 5800.0
    "nx" => nx => 300.0
    "Kq" => Kq => 1.522190403737490e+05
    "Kp" => Kp => 180.1378030928276
    "vt" => vt => 726.0
    "wr" => wr => 929.9678874564831
    "wq" => wq => 948.9349882947897
    "wp" => wp => 0.0
    "nq" => nq => 4
    "nr" => nr => 7549.0
    "dm" => dm => 0.1
    "kb" => kb => 0.0095
    "ku" => ku => 1
    "ns" => ns => 0.5
]);

end