#############################
# run_from_config.jl
#
# Run synthetic_dataset_generator.jl using a YAML config.
#
# Usage:
#   julia run_from_config.jl config.yaml
#
# Expected keys in config.yaml:
#   dataset_name   :: String (optional, used only for logging)
#   model          :: String in ["HETER", "REPR", "NOT", "AND", "NAND"]
#   N              :: Int
#   seed           :: Int (optional, default = 42)
#   outfile        :: String (optional, JLD2 filename)
#   output_dir     :: String (optional, default = "generated_datasets")
#   python_format  :: String (optional, "npz" or "none"; default "npz")
#
#   time:          :: Dict (optional)
#     tspan: [t0, t1]     (optional, default = DEFAULT_TSPAN)
#     n_points: Int       (optional, default = DEFAULT_N_TIME)
#
#   params:        :: Dict (optional)
#     <MODEL>:
#       <param_name>:
#         dist: "linear" | "log10"
#         min:  Float64
#         max:  Float64
#############################

using YAML
using NPZ

# Bring in the generator, MODEL_DICT, config types, and helpers
include("synthetic_dataset_generator.jl")

# ----------------------------
# Build sampling config from YAML
# ----------------------------

"""
    build_sampling_config(cfg::Dict, model_name::String) -> Union{Nothing,ModelSamplingConfig}

Parse the `params` section of the YAML config for the given model name.

Expected structure in YAML:

params:
  AND:
    ns:
      dist: linear
      min: 0.1
      max: 1.0
    kappa_ini:
      dist: log10
      min: -0.65
      max: 0.0
    wmaxrep_1:
      dist: log10
      min: 0.0
      max: 3.0
    # ...

If the block or model subsection is missing, returns `nothing` and defaults
are used in the sampler.
"""
function build_sampling_config(cfg::Dict, model_name::String)
    # Top-level "params" may or may not exist
    params_by_model = get(cfg, "params", nothing)
    params_by_model === nothing && return nothing

    model_block_any = get(params_by_model, model_name, nothing)
    model_block_any === nothing && return nothing

    # YAML.jl returns Dict{Any,Any}
    param_dict = Dict{Symbol,ParamRange}()

    for (pname_any, spec_any) in model_block_any
        # pname_any is likely a String
        pname = Symbol(String(pname_any))

        # spec_any should be a Dict with "dist", "min", "max"
        dist_str = String(get(spec_any, "dist", "linear"))
        dist_sym = dist_str == "log10" ? :log10 : :linear

        min_val  = Float64(get(spec_any, "min", 0.0))
        max_val  = Float64(get(spec_any, "max", 1.0))

        param_dict[pname] = ParamRange(dist_sym, min_val, max_val)
    end

    return ModelSamplingConfig(param_dict)
end

# ----------------------------
# Main entrypoint
# ----------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        error("Usage: julia run_from_config.jl <config.yaml>")
    end

    cfg_path = ARGS[1]
    cfg = YAML.load_file(cfg_path)

    dataset_name = get(cfg, "dataset_name", "synthetic_dataset")
    model_name   = uppercase(get(cfg, "model", "AND"))
    N            = get(cfg, "N", 200)
    seed         = get(cfg, "seed", 42)

    # Directory where all results will be stored
    output_dir   = get(cfg, "output_dir", "generated_datasets")
    python_fmt   = get(cfg, "python_format", "npz")  # "npz" or "none"

    # Time configuration (optional)
    time_cfg = get(cfg, "time", nothing)
    tspan = DEFAULT_TSPAN
    n_time = DEFAULT_N_TIME

    if time_cfg !== nothing
        # Expect tspan: [t0, t1]
        tspan_vec = get(time_cfg, "tspan", [tspan[1], tspan[2]])
        if length(tspan_vec) != 2
            error("Config 'time.tspan' must be a 2-element array [t0, t1].")
        end
        tspan = (Float64(tspan_vec[1]), Float64(tspan_vec[2]))
        n_time = Int(get(time_cfg, "n_points", n_time))
    end

    # derive default outfile (JLD2) if not provided
    outfile_cfg  = get(cfg, "outfile", "$(lowercase(dataset_name)).jld2")
    outfile_base = String(outfile_cfg)
    if !endswith(outfile_base, ".jld2")
        outfile_base *= ".jld2"
    end

    # check model
    if !haskey(MODEL_DICT, model_name)
        valid = join(collect(keys(MODEL_DICT)), ", ")
        error("Unknown model '$model_name' in config. Use one of: $valid.")
    end

    model_def = MODEL_DICT[model_name]

    # Build sampling config from YAML, if provided
    sampling_cfg = build_sampling_config(cfg, model_name)

    # Ensure output directory exists
    mkpath(output_dir)

    # Final paths
    jld_path = joinpath(output_dir, outfile_base)
    npz_path = replace(jld_path, r"\.jld2$" => ".npz")

    println("▶ Generating dataset '$dataset_name'")
    println("  model        = $model_name")
    println("  N            = $N")
    println("  seed         = $seed")
    println("  output_dir   = $output_dir")
    println("  jld2_outfile = $jld_path")
    println("  python_fmt   = $python_fmt")
    println("  tspan        = $tspan")
    println("  n_time       = $n_time")

    # Generate data with explicit time grid
    X, Z, T = generate_dataset_with_time(
        model_def,
        N;
        outfile      = jld_path,
        seed         = seed,
        sampling_cfg = sampling_cfg,
        tspan        = tspan,
        n_time       = n_time,
    )

    # Names for Z (design/env features); currently only used on the Julia side
    z_names = z_feature_names(model_def)

    if lowercase(python_fmt) == "npz"
        # Export Python-friendly NPZ (numeric-only; NPZ.jl does not support String arrays)
        # Python usage:
        #   import numpy as np
        #   data  = np.load("path/to/file.npz")
        #   X     = data["X"]      # (N, T, D)
        #   Z     = data["Z"]      # (N, P)
        #   T     = data["T"]      # (T,)
        #   N     = int(data["N"])
        #   tspan = data["tspan"]  # array([t0, t1])
        NPZ.npzwrite(
            npz_path,
            Dict(
                "X"     => X,
                "Z"     => Z,
                "T"     => T,
                "N"     => N,
                "tspan" => [tspan[1], tspan[2]],
            ),
        )
        println("Also wrote Python-compatible NPZ to '$npz_path'")
    else
        println("Skipping Python NPZ export (python_format='$python_fmt').")
    end
end
