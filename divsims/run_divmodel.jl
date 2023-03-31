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
                println("---------------Birth...t=$(t[end]), gen = $k, forks: $(length(init)-26)")
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

            tf = tbirth+tint
            init = init[1:26] # remove oldest fork            
            ftmp = zeros(10)
            ftmp[1:length(init)-26] = init[27:end]
            forks = hcat(forks, ftmp)

            if part == "stoch1"
                # just partition binomially with pdiv
                ppart = pdiv
            elseif part == "stoch2"
                # first sample size partitioning and then partition
                # species accordingly
                alpha = 200 # chosen to give CV of ~5% as reported in Koppes, L., C. Woldringh, and N. Nanninga. 1978. Size variations and correlation of different cell cycle events in slow growing Escherichia coli. J. Bacteriol. 134:423?433.
                beta = alpha/pdiv - alpha
                ppart = rand(Beta(alpha, beta))
            end

            if part == "det"
                # partition deterministically
                init[1:26] = y[1:26,end]*pdiv
            else
                neg = true
                while neg
                    # sample until no negative initial values there anymore (shouldnt happen too often)
                    binidx = sqrt.(round.(y[1:26,end])) .< 3 * sqrt(ppart*(1-ppart)) / ppart
                    normidx = .!binidx
                    init[normidx] = randn(sum(normidx)) .* sqrt(y[normidx,end]*ppart*(1-ppart)) + y[normidx,end]*ppart
                    if adet
                        # partition a deterministically
                        binidx[varidx["a"]] .= false
                        init[varidx["a"]] = y[varidx["a"],end] * ppart
                    end
                    init[binidx] = binomial(round.(y[binidx,end]), ppart)
                    neg = sum(init[1:26] .< 0) > 0
                end
            end

            # remember values after division
            for j in 1:length(varnames)
                eval(Meta.parse("$(varnames[j]) = hcat($varnames[j], init[j])"))
            end
            
            M = nr*(r[end] + rmg[end] + rmy[end] + rmr[end] + rmt[end] + rmm[end] + rmq[end] + zmg[end] + zmy[end] + zmr[end] + zmt[end] + zmm[end] + zmq[end]) + nx * (gfp[end] + yfp[end] + q[end] + et[end] + em[end])
        
            Tvec = vcat(Tvec, t[end])
            tvec = hcat(tvec, t)
        
            if any(isnan.(init))
                error("NaNs detected after partitioning!")
            end
        end
    end
    IE = 0
end
toc()

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