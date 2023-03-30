# vectors for saving thing
tvec = []
tauvec = []
Tvec = []
birthidx = []
tinit = []
t0 = 0
tint = 1e4 # max time for cell cycle
tf = tint
forks = []

data[:p] = parameters
data[:tframe] = tframe
data[:rates] = rates

IE = 0
tbirth = 0
t = 0
tic()

for k = 1:length(varnames)
eval(Symbol(varnames[k])) = []
end

# simulate
for k = 1:Ngens
    tcycle = []
    if 1/Vcrit - 2^(length(init)-26)/M > 0
        # at division it can happen that volume is randomly larger such
        # that the rep-initiation criterium is already fulfilled
        # => start new round of replication right away
        push!(init, 0)
        tinit = vcat(tinit, t[end])
        ftmp = zeros(10, 1)
        ftmp[1:length(init)-26] = init[27:end]
        forks = hcat(forks, ftmp)
        if printout
            println("Initiating right after division... t=$(t[end]), forks: $(length(init)-1)")
        end
    end
end


while isempty(intersect(IE,2))   # while not divided
    options = CVodeSetOptions("UserData", data,
                              "RelTol", 1e-6,
                              "AbsTol", 1e-12 .* ones(size(init)),
                              "LinearSolver", "Dense",
                              "InitialStep", 1e-4,
                              "MaxStep", 0.1,
                              "MaxNumSteps", 1e5,
                              "RootsFn", divmodel_ev, "NumRoots", 3)
    CVodeInit(divmodel_odes, "BDF", "Newton", t0, init, options)
    status, t, y = CVode([tf], "Normal")

    for j in 1:length(varnames)
        eval(parse("$(varnames[j]) = cat(2, $(varnames[j]), y[j, :]))"))
    end

    evstate = divmodel_ev(t, y, data)
    IE = find(abs.(evstate) ./ [1e-8 C+D 1] .< evacc)
    init = y[:, end]

    if !isempty(intersect(IE, 1))
        init = vcat(init, 0)
        tinit = cat(1, tinit, t[end])
        ftmp = zeros(10, 1)
        ftmp[1:length(init)-26] = init[27:end]
        forks = cat(2, forks, ftmp)
        if printout
            println("Initiating... t=$(t[end]), forks: $(length(init)-26)")
        end
    end

    tvec = cat(2, tvec, t)
    tcycle = cat(2, tcycle, t)
    t0 = t[end]
    if !isempty(intersect(IE, 2))
        if printout
            println("---------------Birth...t=$(t[end]), gen = $(k), forks: $(length(init)-26)")
        end

        birthidx = vcat(birthidx, length(tvec)+1)
        if !isempty(tcycle)
            tauvec = vcat(tauvec, tvec[end]-tbirth)
            tbirth = tcycle[end]
            Tvec = vcat(Tvec, tvec[end])
        else
            tauvec = vcat(tauvec, t)
            tbirth = t[end]
            Tvec = vcat(Tvec, t[end])
        end
        tf = tbirth + tint
        init[27] = []       # remove oldest fork
        ftmp = zeros(10, 1)
        ftmp[1:length(init)-26] = init[27:end]
        forks = cat(2, forks, ftmp)
        ppart = 0
        alpha = 200
        beta = alpha / pdiv - alpha
        # HERE :: change syntax to be Julia-specific
        # pattern matching :: https://thautwarm.github.io/MLStyle.jl/latest/syntax/pattern.html#guards
        # also check Match.jl :: https://github.com/kmsquire/Match.jl
        switch part
            case "stoch1"
                # just partition binomially with pdiv
                ppart = pdiv
            case "stoch2"
                ppart = betarnd(alpha, beta)
        end
        switch part
            case "det"
                # partition deterministically
                init[1:26] = y[1:26, end] .* pdiv
            otherwise
                neg = true
                while neg
                    binidx = sqrt.(round.(y[1:26, end])) .< 3 * sqrt(ppart*(1-ppart)) / ppart
                    normidx = .!binidx
                    init[normidx] = randn(sum(normidx), 1)


%%%
# Post calculations
# Get rid of last stored birth index
birthidx = birthidx[1:end-1]

# Division indices are those just before birth
dividx = birthidx .- 1

# Compute protein mass and growth rate
M = nr * (r + rmg + rmy + rmr + rmt + rmm + rmq + zmg + zmy + zmr + zmt + zmm + zmq) +
    nx * (gfp + yfp + q + et + em)
Kgamma = gmax / Kp * M / Mref
gamma = gmax .* a ./ (Kgamma .+ a)
ttrate = (rmg .+ rmy .+ rmq .+ rmr .+ rmt .+ rmm) .* gamma
lam = ttrate ./ M

# Birth, division, and added size
vb = M[birthidx]        # Birth size
vd = M[dividx]          # Division size
addv = vd[2:end] .- vb[1:end-1]