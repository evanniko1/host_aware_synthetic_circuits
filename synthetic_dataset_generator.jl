#############################
# synthetic_dataset_generator.jl
#
# Generic sampler + dataset generator for host-aware synthetic circuits:
#   HETER_ODE_model!, REPR_ODE_model!, NOT_gate_ODE_model!,
#   AND_gate_ODE_model!, NAND_gate_ODE_model!
#
# Outputs:
#   X :: Array{Float64,3}  (N, T, D)  – all state variables over time
#   Z :: Array{Float64,2}  (N, P)     – design+environment features
#   (optionally) T :: Vector{Float64} – time grid
#
# Usage from CLI:
#   julia synthetic_dataset_generator.jl AND 500 and_500.jld2 123
#############################

using DifferentialEquations
using Sundials
using Random
using JLD2

# bring in your existing code
include("values.jl")            # defines u0, p, dmrep, dprep, etc.
include("helper.jl")            # create_problem_dict!, solve_ode_problem!, ...
include("host_aware_models.jl") # HETER_ODE_model!, REPR_ODE_model!, ...

# ----------------------------
# Global simulation defaults
# ----------------------------

const DEFAULT_TSPAN  = (0.0, 1e7)
const DEFAULT_N_TIME = 200

abstract type AbstractDesignEnv end

# ----------------------------
# Design + environment structs
# ----------------------------

struct HeterDesignEnv <: AbstractDesignEnv
    ns::Float64
    kappa_ini::Float64
    wmaxrep::Float64
    kbrep::Float64
    kurep::Float64
end

struct ReprDesignEnv  <: AbstractDesignEnv
    ns::Float64
    kappa_ini::Float64
    wmaxrep_1::Float64; kbrep_1::Float64; kurep_1::Float64
    wmaxrep_2::Float64; kbrep_2::Float64; kurep_2::Float64
    wmaxrep_3::Float64; kbrep_3::Float64; kurep_3::Float64
    Kq_rep_1::Float64; nq_rep_1::Float64
    Kq_rep_2::Float64; nq_rep_2::Float64
    Kq_rep_3::Float64; nq_rep_3::Float64
end

struct NotDesignEnv   <: AbstractDesignEnv
    ns::Float64
    kappa_ini::Float64
    wmaxrep_1::Float64; kbrep_1::Float64; kurep_1::Float64
    wmaxrep_2::Float64; kbrep_2::Float64; kurep_2::Float64
    Kq_rep_1::Float64; nq_rep_1::Float64
end

struct AndDesignEnv   <: AbstractDesignEnv
    ns::Float64
    kappa_ini::Float64
    wmaxrep_1::Float64; kbrep_1::Float64; kurep_1::Float64
    wmaxrep_2::Float64; kbrep_2::Float64; kurep_2::Float64
    wmaxrep_3::Float64; kbrep_3::Float64; kurep_3::Float64
    Kq_rep_1::Float64; nq_rep_1::Float64
    Kq_rep_2::Float64; nq_rep_2::Float64
end

struct NandDesignEnv  <: AbstractDesignEnv
    ns::Float64
    kappa_ini::Float64
    wmaxrep_1::Float64; kbrep_1::Float64; kurep_1::Float64
    wmaxrep_2::Float64; kbrep_2::Float64; kurep_2::Float64
    wmaxrep_3::Float64; kbrep_3::Float64; kurep_3::Float64
    wmaxrep_4::Float64; kbrep_4::Float64; kurep_4::Float64
    Kq_rep_1::Float64; nq_rep_1::Float64
    Kq_rep_2::Float64; nq_rep_2::Float64
    Kq_rep_3::Float64; nq_rep_3::Float64
end

# ----------------------------
# Sampling config types
# ----------------------------

struct ParamRange
    dist::Symbol      # :linear or :log10
    min::Float64
    max::Float64
end

struct ModelSamplingConfig
    params::Dict{Symbol,ParamRange}
end

# ----------------------------
# Helper samplers (log / linear)
# ----------------------------

log10_sample(a::Real, b::Real) = 10.0 ^ (a + (b - a) * rand())
linear_sample(a::Real, b::Real) = a + (b - a) * rand()

sample_from_range(r::ParamRange) =
    r.dist === :linear ? linear_sample(r.min, r.max) : log10_sample(r.min, r.max)

function get_param_range(
    cfg::Union{Nothing,ModelSamplingConfig},
    name::Symbol,
    default_dist::Symbol,
    default_min::Float64,
    default_max::Float64,
)
    if cfg === nothing || !haskey(cfg.params, name)
        return ParamRange(default_dist, default_min, default_max)
    else
        return cfg.params[name]
    end
end

# ----------------------------
# Model-specific samplers (config-aware)
# ----------------------------

function sample_design_env(
    ::typeof(HETER_ODE_model!),
    cfg::Union{Nothing,ModelSamplingConfig},
)::HeterDesignEnv
    ns_rng        = get_param_range(cfg, :ns,        :linear, 0.1,   1.0)
    kappa_ini_rng = get_param_range(cfg, :kappa_ini, :log10,  -0.65, 0.0)
    wmaxrep_rng   = get_param_range(cfg, :wmaxrep,   :log10,  0.0,   3.5)
    kbrep_rng     = get_param_range(cfg, :kbrep,     :log10, -3.0,  -1.0)
    kurep_rng     = get_param_range(cfg, :kurep,     :log10, -3.0,  -1.0)

    ns        = sample_from_range(ns_rng)
    kappa_ini = sample_from_range(kappa_ini_rng)
    wmaxrep   = sample_from_range(wmaxrep_rng)
    kbrep     = sample_from_range(kbrep_rng)
    kurep     = sample_from_range(kurep_rng)

    return HeterDesignEnv(ns, kappa_ini, wmaxrep, kbrep, kurep)
end

function sample_design_env(
    ::typeof(REPR_ODE_model!),
    cfg::Union{Nothing,ModelSamplingConfig},
)::ReprDesignEnv
    ns_rng        = get_param_range(cfg, :ns,        :linear, 0.1,   1.0)
    kappa_ini_rng = get_param_range(cfg, :kappa_ini, :log10,  -0.65, 0.0)

    wmaxrep_1_rng = get_param_range(cfg, :wmaxrep_1, :log10, 0.0, 3.0)
    wmaxrep_2_rng = get_param_range(cfg, :wmaxrep_2, :log10, 0.0, 3.0)
    wmaxrep_3_rng = get_param_range(cfg, :wmaxrep_3, :log10, 0.0, 3.5)

    kbrep_1_rng   = get_param_range(cfg, :kbrep_1, :log10, -3.0, -1.0)
    kbrep_2_rng   = get_param_range(cfg, :kbrep_2, :log10, -3.0, -1.0)
    kbrep_3_rng   = get_param_range(cfg, :kbrep_3, :log10, -3.0, -1.0)

    kurep_1_rng   = get_param_range(cfg, :kurep_1, :log10, -3.0, -1.0)
    kurep_2_rng   = get_param_range(cfg, :kurep_2, :log10, -3.0, -1.0)
    kurep_3_rng   = get_param_range(cfg, :kurep_3, :log10, -3.0, -1.0)

    Kq_rep_1_rng  = get_param_range(cfg, :Kq_rep_1, :log10, 1.7, 4.0)
    Kq_rep_2_rng  = get_param_range(cfg, :Kq_rep_2, :log10, 1.7, 4.0)
    Kq_rep_3_rng  = get_param_range(cfg, :Kq_rep_3, :log10, 1.7, 4.0)

    nq_rep_1_rng  = get_param_range(cfg, :nq_rep_1, :linear, 1.0, 4.0)
    nq_rep_2_rng  = get_param_range(cfg, :nq_rep_2, :linear, 1.0, 4.0)
    nq_rep_3_rng  = get_param_range(cfg, :nq_rep_3, :linear, 1.0, 4.0)

    ns        = sample_from_range(ns_rng)
    kappa_ini = sample_from_range(kappa_ini_rng)

    wmaxrep_1 = sample_from_range(wmaxrep_1_rng)
    wmaxrep_2 = sample_from_range(wmaxrep_2_rng)
    wmaxrep_3 = sample_from_range(wmaxrep_3_rng)

    kbrep_1   = sample_from_range(kbrep_1_rng)
    kbrep_2   = sample_from_range(kbrep_2_rng)
    kbrep_3   = sample_from_range(kbrep_3_rng)

    kurep_1   = sample_from_range(kurep_1_rng)
    kurep_2   = sample_from_range(kurep_2_rng)
    kurep_3   = sample_from_range(kurep_3_rng)

    Kq_rep_1  = sample_from_range(Kq_rep_1_rng)
    Kq_rep_2  = sample_from_range(Kq_rep_2_rng)
    Kq_rep_3  = sample_from_range(Kq_rep_3_rng)

    nq_rep_1  = sample_from_range(nq_rep_1_rng)
    nq_rep_2  = sample_from_range(nq_rep_2_rng)
    nq_rep_3  = sample_from_range(nq_rep_3_rng)

    return ReprDesignEnv(ns, kappa_ini,
                         wmaxrep_1, kbrep_1, kurep_1,
                         wmaxrep_2, kbrep_2, kurep_2,
                         wmaxrep_3, kbrep_3, kurep_3,
                         Kq_rep_1, nq_rep_1,
                         Kq_rep_2, nq_rep_2,
                         Kq_rep_3, nq_rep_3)
end

function sample_design_env(
    ::typeof(NOT_gate_ODE_model!),
    cfg::Union{Nothing,ModelSamplingConfig},
)::NotDesignEnv
    ns_rng        = get_param_range(cfg, :ns,        :linear, 0.1,   1.0)
    kappa_ini_rng = get_param_range(cfg, :kappa_ini, :log10,  -0.65, 0.0)

    wmaxrep_1_rng = get_param_range(cfg, :wmaxrep_1, :log10, 0.0, 3.0)
    wmaxrep_2_rng = get_param_range(cfg, :wmaxrep_2, :log10, 0.0, 3.5)

    kbrep_1_rng   = get_param_range(cfg, :kbrep_1, :log10, -3.0, -1.0)
    kbrep_2_rng   = get_param_range(cfg, :kbrep_2, :log10, -3.0, -1.0)

    kurep_1_rng   = get_param_range(cfg, :kurep_1, :log10, -3.0, -1.0)
    kurep_2_rng   = get_param_range(cfg, :kurep_2, :log10, -3.0, -1.0)

    Kq_rep_1_rng  = get_param_range(cfg, :Kq_rep_1, :log10, 1.7, 4.0)
    nq_rep_1_rng  = get_param_range(cfg, :nq_rep_1, :linear, 1.0, 4.0)

    ns        = sample_from_range(ns_rng)
    kappa_ini = sample_from_range(kappa_ini_rng)

    wmaxrep_1 = sample_from_range(wmaxrep_1_rng)
    wmaxrep_2 = sample_from_range(wmaxrep_2_rng)

    kbrep_1   = sample_from_range(kbrep_1_rng)
    kbrep_2   = sample_from_range(kbrep_2_rng)

    kurep_1   = sample_from_range(kurep_1_rng)
    kurep_2   = sample_from_range(kurep_2_rng)

    Kq_rep_1  = sample_from_range(Kq_rep_1_rng)
    nq_rep_1  = sample_from_range(nq_rep_1_rng)

    return NotDesignEnv(ns, kappa_ini,
                        wmaxrep_1, kbrep_1, kurep_1,
                        wmaxrep_2, kbrep_2, kurep_2,
                        Kq_rep_1, nq_rep_1)
end

function sample_design_env(
    ::typeof(AND_gate_ODE_model!),
    cfg::Union{Nothing,ModelSamplingConfig},
)::AndDesignEnv
    ns_rng        = get_param_range(cfg, :ns,        :linear, 0.1,   1.0)
    kappa_ini_rng = get_param_range(cfg, :kappa_ini, :log10,  -0.65, 0.0)

    wmaxrep_1_rng = get_param_range(cfg, :wmaxrep_1, :log10, 0.0, 3.0)
    wmaxrep_2_rng = get_param_range(cfg, :wmaxrep_2, :log10, 0.0, 3.0)
    wmaxrep_3_rng = get_param_range(cfg, :wmaxrep_3, :log10, 0.3, 3.5)

    kbrep_1_rng   = get_param_range(cfg, :kbrep_1, :log10, -3.0, -1.0)
    kbrep_2_rng   = get_param_range(cfg, :kbrep_2, :log10, -3.0, -1.0)
    kbrep_3_rng   = get_param_range(cfg, :kbrep_3, :log10, -3.0, -1.0)

    kurep_1_rng   = get_param_range(cfg, :kurep_1, :log10, -3.0, -1.0)
    kurep_2_rng   = get_param_range(cfg, :kurep_2, :log10, -3.0, -1.0)
    kurep_3_rng   = get_param_range(cfg, :kurep_3, :log10, -3.0, -1.0)

    Kq_rep_1_rng  = get_param_range(cfg, :Kq_rep_1, :log10, 1.7, 4.0)
    Kq_rep_2_rng  = get_param_range(cfg, :Kq_rep_2, :log10, 2.0, 4.2)

    nq_rep_1_rng  = get_param_range(cfg, :nq_rep_1, :linear, 1.0, 4.0)
    nq_rep_2_rng  = get_param_range(cfg, :nq_rep_2, :linear, 1.0, 4.0)

    ns        = sample_from_range(ns_rng)
    kappa_ini = sample_from_range(kappa_ini_rng)

    wmaxrep_1 = sample_from_range(wmaxrep_1_rng)
    wmaxrep_2 = sample_from_range(wmaxrep_2_rng)
    wmaxrep_3 = sample_from_range(wmaxrep_3_rng)

    kbrep_1   = sample_from_range(kbrep_1_rng)
    kbrep_2   = sample_from_range(kbrep_2_rng)
    kbrep_3   = sample_from_range(kbrep_3_rng)

    kurep_1   = sample_from_range(kurep_1_rng)
    kurep_2   = sample_from_range(kurep_2_rng)
    kurep_3   = sample_from_range(kurep_3_rng)

    Kq_rep_1  = sample_from_range(Kq_rep_1_rng)
    Kq_rep_2  = sample_from_range(Kq_rep_2_rng)

    nq_rep_1  = sample_from_range(nq_rep_1_rng)
    nq_rep_2  = sample_from_range(nq_rep_2_rng)

    return AndDesignEnv(ns, kappa_ini,
                        wmaxrep_1, kbrep_1, kurep_1,
                        wmaxrep_2, kbrep_2, kurep_2,
                        wmaxrep_3, kbrep_3, kurep_3,
                        Kq_rep_1, nq_rep_1,
                        Kq_rep_2, nq_rep_2)
end

function sample_design_env(
    ::typeof(NAND_gate_ODE_model!),
    cfg::Union{Nothing,ModelSamplingConfig},
)::NandDesignEnv
    ns_rng        = get_param_range(cfg, :ns,        :linear, 0.1,   1.0)
    kappa_ini_rng = get_param_range(cfg, :kappa_ini, :log10,  -0.65, 0.0)

    wmaxrep_1_rng = get_param_range(cfg, :wmaxrep_1, :log10, 0.0, 3.0)
    wmaxrep_2_rng = get_param_range(cfg, :wmaxrep_2, :log10, 0.0, 3.0)
    wmaxrep_3_rng = get_param_range(cfg, :wmaxrep_3, :log10, 0.3, 3.5)
    wmaxrep_4_rng = get_param_range(cfg, :wmaxrep_4, :log10, 0.0, 3.0)

    kbrep_1_rng   = get_param_range(cfg, :kbrep_1, :log10, -3.0, -1.0)
    kbrep_2_rng   = get_param_range(cfg, :kbrep_2, :log10, -3.0, -1.0)
    kbrep_3_rng   = get_param_range(cfg, :kbrep_3, :log10, -3.0, -1.0)
    kbrep_4_rng   = get_param_range(cfg, :kbrep_4, :log10, -3.0, -1.0)

    kurep_1_rng   = get_param_range(cfg, :kurep_1, :log10, -3.0, -1.0)
    kurep_2_rng   = get_param_range(cfg, :kurep_2, :log10, -3.0, -1.0)
    kurep_3_rng   = get_param_range(cfg, :kurep_3, :log10, -3.0, -1.0)
    kurep_4_rng   = get_param_range(cfg, :kurep_4, :log10, -3.0, -1.0)

    Kq_rep_1_rng  = get_param_range(cfg, :Kq_rep_1, :log10, 1.7, 4.0)
    Kq_rep_2_rng  = get_param_range(cfg, :Kq_rep_2, :log10, 2.0, 4.2)
    Kq_rep_3_rng  = get_param_range(cfg, :Kq_rep_3, :log10, 1.7, 4.0)

    nq_rep_1_rng  = get_param_range(cfg, :nq_rep_1, :linear, 1.0, 4.0)
    nq_rep_2_rng  = get_param_range(cfg, :nq_rep_2, :linear, 1.0, 4.0)
    nq_rep_3_rng  = get_param_range(cfg, :nq_rep_3, :linear, 1.0, 4.0)

    ns        = sample_from_range(ns_rng)
    kappa_ini = sample_from_range(kappa_ini_rng)

    wmaxrep_1 = sample_from_range(wmaxrep_1_rng)
    wmaxrep_2 = sample_from_range(wmaxrep_2_rng)
    wmaxrep_3 = sample_from_range(wmaxrep_3_rng)
    wmaxrep_4 = sample_from_range(wmaxrep_4_rng)

    kbrep_1   = sample_from_range(kbrep_1_rng)
    kbrep_2   = sample_from_range(kbrep_2_rng)
    kbrep_3   = sample_from_range(kbrep_3_rng)
    kbrep_4   = sample_from_range(kbrep_4_rng)

    kurep_1   = sample_from_range(kurep_1_rng)
    kurep_2   = sample_from_range(kurep_2_rng)
    kurep_3   = sample_from_range(kurep_3_rng)
    kurep_4   = sample_from_range(kurep_4_rng)

    Kq_rep_1  = sample_from_range(Kq_rep_1_rng)
    Kq_rep_2  = sample_from_range(Kq_rep_2_rng)
    Kq_rep_3  = sample_from_range(Kq_rep_3_rng)

    nq_rep_1  = sample_from_range(nq_rep_1_rng)
    nq_rep_2  = sample_from_range(nq_rep_2_rng)
    nq_rep_3  = sample_from_range(nq_rep_3_rng)

    return NandDesignEnv(ns, kappa_ini,
                         wmaxrep_1, kbrep_1, kurep_1,
                         wmaxrep_2, kbrep_2, kurep_2,
                         wmaxrep_3, kbrep_3, kurep_3,
                         wmaxrep_4, kbrep_4, kurep_4,
                         Kq_rep_1, nq_rep_1,
                         Kq_rep_2, nq_rep_2,
                         Kq_rep_3, nq_rep_3)
end

# Backwards-compatible convenience: call with no config
sample_design_env(model_def::Function) = sample_design_env(model_def, nothing)

# ----------------------------
# Problem construction per model (tspan-aware)
# ----------------------------

# HETER
function build_problem(::typeof(HETER_ODE_model!), des::HeterDesignEnv, tspan::Tuple{Float64,Float64})
    p_local = copy(p)
    het_p = [
        des.ns, dmrep, dprep, des.kappa_ini,
        des.wmaxrep, des.kbrep, des.kurep
    ]
    append!(p_local, het_p)

    u0_local = copy(u0) # already has one heterologous gene

    return create_problem_dict!(
        model_choice  = HETER_ODE_model!,
        init_values   = u0_local,
        params_values = p_local,
        tspan         = tspan,
        ode_solver    = Rodas4(autodiff = false),
        abstol        = 1e-8,
        reltol        = 1e-8,
        maxiters      = 1e7,
        show_progress = false
    )
end

# REPR
function build_problem(::typeof(REPR_ODE_model!), des::ReprDesignEnv, tspan::Tuple{Float64,Float64})
    p_local = copy(p)
    het_p = [
        des.ns, dmrep, dprep, des.kappa_ini,
        des.wmaxrep_1, des.kbrep_1, des.kurep_1,
        des.wmaxrep_2, des.kbrep_2, des.kurep_2,
        des.wmaxrep_3, des.kbrep_3, des.kurep_3,
        des.Kq_rep_1, des.nq_rep_1,
        des.Kq_rep_2, des.nq_rep_2,
        des.Kq_rep_3, des.nq_rep_3
    ]
    append!(p_local, het_p)

    u0_local = copy(u0)
    # add rep_2, mrep_2, rmrep_2, rep_3, mrep_3, rmrep_3
    u0_extra = [1.0, 1.0, 0.0, 1.0, 1.0, 0.0]
    append!(u0_local, u0_extra)

    return create_problem_dict!(
        model_choice  = REPR_ODE_model!,
        init_values   = u0_local,
        params_values = p_local,
        tspan         = tspan,
        ode_solver    = Rodas4(autodiff = false),
        abstol        = 1e-8,
        reltol        = 1e-8,
        maxiters      = 1e7,
        show_progress = false
    )
end

# NOT
function build_problem(::typeof(NOT_gate_ODE_model!), des::NotDesignEnv, tspan::Tuple{Float64,Float64})
    p_local = copy(p)
    het_p = [
        des.ns, dmrep, dprep, des.kappa_ini,
        des.wmaxrep_1, des.kbrep_1, des.kurep_1,
        des.wmaxrep_2, des.kbrep_2, des.kurep_2,
        des.Kq_rep_1,  des.nq_rep_1
    ]
    append!(p_local, het_p)

    u0_local = copy(u0)
    # add rep_2, mrep_2, rmrep_2
    u0_extra = [0.0, 0.0, 0.0]
    append!(u0_local, u0_extra)

    return create_problem_dict!(
        model_choice  = NOT_gate_ODE_model!,
        init_values   = u0_local,
        params_values = p_local,
        tspan         = tspan,
        ode_solver    = Rodas4(autodiff = false),
        abstol        = 1e-8,
        reltol        = 1e-8,
        maxiters      = 1e7,
        show_progress = false
    )
end

# AND
function build_problem(::typeof(AND_gate_ODE_model!), des::AndDesignEnv, tspan::Tuple{Float64,Float64})
    p_local = copy(p)
    het_p = [
        des.ns, dmrep, dprep, des.kappa_ini,
        des.wmaxrep_1, des.kbrep_1, des.kurep_1,
        des.wmaxrep_2, des.kbrep_2, des.kurep_2,
        des.wmaxrep_3, des.kbrep_3, des.kurep_3,
        des.Kq_rep_1,  des.nq_rep_1,
        des.Kq_rep_2,  des.nq_rep_2
    ]
    append!(p_local, het_p)

    u0_local = copy(u0)
    # add rep_2, mrep_2, rmrep_2, rep_3, mrep_3, rmrep_3
    u0_extra = [1.0, 1.0, 0.0, 1.0, 1.0, 0.0]
    append!(u0_local, u0_extra)

    return create_problem_dict!(
        model_choice  = AND_gate_ODE_model!,
        init_values   = u0_local,
        params_values = p_local,
        tspan         = tspan,
        ode_solver    = Rodas4(autodiff = false),
        abstol        = 1e-8,
        reltol        = 1e-8,
        maxiters      = 1e7,
        show_progress = false
    )
end

# NAND
function build_problem(::typeof(NAND_gate_ODE_model!), des::NandDesignEnv, tspan::Tuple{Float64,Float64})
    p_local = copy(p)
    het_p = [
        des.ns, dmrep, dprep, des.kappa_ini,
        des.wmaxrep_1, des.kbrep_1, des.kurep_1,
        des.wmaxrep_2, des.kbrep_2, des.kurep_2,
        des.wmaxrep_3, des.kbrep_3, des.kurep_3,
        des.wmaxrep_4, des.kbrep_4, des.kurep_4,
        des.Kq_rep_1, des.nq_rep_1,
        des.Kq_rep_2, des.nq_rep_2,
        des.Kq_rep_3, des.nq_rep_3
    ]
    append!(p_local, het_p)

    u0_local = copy(u0)
    # base u0 has gene 1; add genes 2–4
    u0_extra = [1.0, 100.0, 0.0, 100.0, 1.0, 0.0, 1.0, 100.0, 0.0]
    append!(u0_local, u0_extra)

    return create_problem_dict!(
        model_choice  = NAND_gate_ODE_model!,
        init_values   = u0_local,
        params_values = p_local,
        tspan         = tspan,
        ode_solver    = Rodas4(autodiff = false),
        abstol        = 1e-8,
        reltol        = 1e-8,
        maxiters      = 1e7,
        show_progress = false
    )
end

# ----------------------------
# Common simulation logic (tspan + grid)
# ----------------------------

function simulate_trajectory(
    model_def::Function,
    des::AbstractDesignEnv,
    tspan::Tuple{Float64,Float64},
    t_grid,
)
    ode_dict = build_problem(model_def, des, tspan)
    sol      = solve_ode_problem!(ode_problem_wrap = ode_dict)

    D = length(sol.u[end])
    X = Array{Float64}(undef, length(t_grid), D)

    for (i, t) in enumerate(t_grid)
        X[i, :] = sol(t)
    end
    return X
end

# ----------------------------
# Encode design/env to feature vector
# ----------------------------

function encode_design_env(::typeof(HETER_ODE_model!), des::HeterDesignEnv)
    [
        log10(des.ns),
        log10(des.kappa_ini),
        log10(des.wmaxrep),
        log10(des.kbrep),
        log10(des.kurep)
    ]
end

function encode_design_env(::typeof(REPR_ODE_model!), des::ReprDesignEnv)
    [
        log10(des.ns),
        log10(des.kappa_ini),

        log10(des.wmaxrep_1), log10(des.kbrep_1), log10(des.kurep_1),
        log10(des.wmaxrep_2), log10(des.kbrep_2), log10(des.kurep_2),
        log10(des.wmaxrep_3), log10(des.kbrep_3), log10(des.kurep_3),

        log10(des.Kq_rep_1), des.nq_rep_1,
        log10(des.Kq_rep_2), des.nq_rep_2,
        log10(des.Kq_rep_3), des.nq_rep_3
    ]
end

function encode_design_env(::typeof(NOT_gate_ODE_model!), des::NotDesignEnv)
    [
        log10(des.ns),
        log10(des.kappa_ini),

        log10(des.wmaxrep_1), log10(des.kbrep_1), log10(des.kurep_1),
        log10(des.wmaxrep_2), log10(des.kbrep_2), log10(des.kurep_2),

        log10(des.Kq_rep_1), des.nq_rep_1
    ]
end

function encode_design_env(::typeof(AND_gate_ODE_model!), des::AndDesignEnv)
    [
        log10(des.ns),
        log10(des.kappa_ini),

        log10(des.wmaxrep_1), log10(des.kbrep_1), log10(des.kurep_1),
        log10(des.wmaxrep_2), log10(des.kbrep_2), log10(des.kurep_2),
        log10(des.wmaxrep_3), log10(des.kbrep_3), log10(des.kurep_3),

        log10(des.Kq_rep_1), des.nq_rep_1,
        log10(des.Kq_rep_2), des.nq_rep_2
    ]
end

function encode_design_env(::typeof(NAND_gate_ODE_model!), des::NandDesignEnv)
    [
        log10(des.ns),
        log10(des.kappa_ini),

        log10(des.wmaxrep_1), log10(des.kbrep_1), log10(des.kurep_1),
        log10(des.wmaxrep_2), log10(des.kbrep_2), log10(des.kurep_2),
        log10(des.wmaxrep_3), log10(des.kbrep_3), log10(des.kurep_3),
        log10(des.wmaxrep_4), log10(des.kbrep_4), log10(des.kurep_4),

        log10(des.Kq_rep_1), des.nq_rep_1,
        log10(des.Kq_rep_2), des.nq_rep_2,
        log10(des.Kq_rep_3), des.nq_rep_3
    ]
end

# ----------------------------
# Z feature names per model (for metadata)
# ----------------------------

function z_feature_names(::typeof(HETER_ODE_model!))
    [
        "log10(ns)",
        "log10(kappa_ini)",
        "log10(wmaxrep)",
        "log10(kbrep)",
        "log10(kurep)",
    ]
end

function z_feature_names(::typeof(REPR_ODE_model!))
    [
        "log10(ns)",
        "log10(kappa_ini)",

        "log10(wmaxrep_1)", "log10(kbrep_1)", "log10(kurep_1)",
        "log10(wmaxrep_2)", "log10(kbrep_2)", "log10(kurep_2)",
        "log10(wmaxrep_3)", "log10(kbrep_3)", "log10(kurep_3)",

        "log10(Kq_rep_1)", "nq_rep_1",
        "log10(Kq_rep_2)", "nq_rep_2",
        "log10(Kq_rep_3)", "nq_rep_3",
    ]
end

function z_feature_names(::typeof(NOT_gate_ODE_model!))
    [
        "log10(ns)",
        "log10(kappa_ini)",

        "log10(wmaxrep_1)", "log10(kbrep_1)", "log10(kurep_1)",
        "log10(wmaxrep_2)", "log10(kbrep_2)", "log10(kurep_2)",

        "log10(Kq_rep_1)", "nq_rep_1",
    ]
end

function z_feature_names(::typeof(AND_gate_ODE_model!))
    [
        "log10(ns)",
        "log10(kappa_ini)",

        "log10(wmaxrep_1)", "log10(kbrep_1)", "log10(kurep_1)",
        "log10(wmaxrep_2)", "log10(kbrep_2)", "log10(kurep_2)",
        "log10(wmaxrep_3)", "log10(kbrep_3)", "log10(kurep_3)",

        "log10(Kq_rep_1)", "nq_rep_1",
        "log10(Kq_rep_2)", "nq_rep_2",
    ]
end

function z_feature_names(::typeof(NAND_gate_ODE_model!))
    [
        "log10(ns)",
        "log10(kappa_ini)",

        "log10(wmaxrep_1)", "log10(kbrep_1)", "log10(kurep_1)",
        "log10(wmaxrep_2)", "log10(kbrep_2)", "log10(kurep_2)",
        "log10(wmaxrep_3)", "log10(kbrep_3)", "log10(kurep_3)",
        "log10(wmaxrep_4)", "log10(kbrep_4)", "log10(kurep_4)",

        "log10(Kq_rep_1)", "nq_rep_1",
        "log10(Kq_rep_2)", "nq_rep_2",
        "log10(Kq_rep_3)", "nq_rep_3",
    ]
end

# Fallback (should never be used if MODEL_DICT is consistent)
function z_feature_names(model_def::Function)
    error("No z_feature_names defined for model $(model_def).")
end

# ----------------------------
# Internal dataset generator (returns X, Z, T)
# ----------------------------

function _generate_dataset(
    model_def::Function,
    N::Int;
    outfile::AbstractString = "synthetic_dataset.jld2",
    seed::Int = 42,
    sampling_cfg::Union{Nothing,ModelSamplingConfig} = nothing,
    tspan::Tuple{Float64,Float64} = DEFAULT_TSPAN,
    n_time::Int = DEFAULT_N_TIME,
)
    Random.seed!(seed)

    t_grid = collect(range(tspan[1], tspan[2], length = n_time))

    # First sample fixes shapes
    des1 = sample_design_env(model_def, sampling_cfg)
    X1   = simulate_trajectory(model_def, des1, tspan, t_grid)
    T, D = size(X1)
    P    = length(encode_design_env(model_def, des1))

    X = Array{Float64,3}(undef, N, T, D)
    Z = Array{Float64,2}(undef, N, P)

    X[1, :, :] .= X1
    Z[1, :]    .= encode_design_env(model_def, des1)

    n = 2
    while n <= N
        try
            des  = sample_design_env(model_def, sampling_cfg)
            traj = simulate_trajectory(model_def, des, tspan, t_grid)

            X[n, :, :] .= traj
            Z[n, :]    .= encode_design_env(model_def, des)
            n += 1
        catch err
            @warn "Sample $n failed with $err — resampling."
        end
    end

    # Store T_GRID as `T_GRID` in JLD2 for backward compatibility
    @save outfile X Z T_GRID=t_grid
    println("Saved dataset for $(model_def) with size: X=$(size(X)), Z=$(size(Z)) to '$outfile'")

    return X, Z, t_grid
end

# ----------------------------
# Public API
# ----------------------------

"""
    generate_dataset(model_def, N; kwargs...) -> X, Z

Backwards-compatible API: returns only X and Z, using any given kwargs
(outfile, seed, sampling_cfg, tspan, n_time).
"""
function generate_dataset(
    model_def::Function,
    N::Int;
    kwargs...
)
    X, Z, _ = _generate_dataset(model_def, N; kwargs...)
    return X, Z
end

"""
    generate_dataset_with_time(model_def, N; kwargs...) -> X, Z, T

Full API: returns X, Z, and the time grid T.
"""
function generate_dataset_with_time(
    model_def::Function,
    N::Int;
    kwargs...
)
    return _generate_dataset(model_def, N; kwargs...)
end

# ----------------------------
# Simple CLI interface
# ----------------------------

const MODEL_DICT = Dict{String,Function}(
    "HETER" => HETER_ODE_model!,
    "REPR"  => REPR_ODE_model!,
    "NOT"   => NOT_gate_ODE_model!,
    "AND"   => AND_gate_ODE_model!,
    "NAND"  => NAND_gate_ODE_model!,
)

if abspath(PROGRAM_FILE) == @__FILE__
    # ARGS: MODEL N [outfile] [seed]
    model_name = length(ARGS) >= 1 ? uppercase(ARGS[1]) : "AND"
    N          = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200
    outfile    = length(ARGS) >= 3 ? ARGS[3] : lowercase(model_name) * "_$(N).jld2"
    seed       = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 42

    if !haskey(MODEL_DICT, model_name)
        valid = join(collect(keys(MODEL_DICT)), ", ")
        error("Unknown model '$model_name'. Use one of: $valid.")
    end

    model_def = MODEL_DICT[model_name]

    # CLI uses default tspan / n_time and the backwards-compatible API
    generate_dataset(model_def, N; outfile = outfile, seed = seed)
end
