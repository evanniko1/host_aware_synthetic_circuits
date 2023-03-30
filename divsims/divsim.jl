using Random

# Clear all variables
for name in names(Main)
    if name != :Random && name != :julia && name != :Base
        Base.eval(Main, Expr(:quote, :($name = nothing)))
    end
end

# Set clearfigs to 0, to not clear the figures
clearfigs = 0

# Set accuracy of event location
evacc = 1e-3

# Set number of generations to simulate
Ngens = 1e2

# Set time frame in minutes at which variables should be 'measured'
tframe = 1e5

# Set probability of symmetric division
pdiv = 0.5

# Set cell partitioning at division ('det', 'stoch1', or 'stoch2')
part = "stoch2"

# Set printout to 0, to not print output for each birth & initiation event
printout = 0

# Set rescale to 1, to rescale distributions when plotting
rescale = 1

# Set partition energy (a) deterministically
adet = 0

# Create vector of ns values
ns_vec = range(0.1, 0.5, length=5)

# Set set_wg to 0
set_wg = 0

# Set fileprefix
fileprefix = "sims/test_$(part)_adet$(adet)"

# Run the model for several nutrient conditions
for nsi in eachindex(ns_vec)
    # Load parameters and initial values
    include("initcond.jl")
    M = (nr * (r + rmg + rmy + rmr + rmt + rmm + rmq + zmg + zmy + zmr + zmt + zmm + zmq) +
         nx * (gfp + yfp + q + et + em))
    ns = ns_vec[nsi]
    include("run_divmodel.jl")
    # Save the results
    save("$(fileprefix)_ns$(ns).mat", Dict(
        "parameters" => parameters,
        "results" => (n, cells, bcells, ecells, qcells, pcells, ccells, tcells, M_cells)
    ))
end