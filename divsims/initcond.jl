# using Symbolics #not sure if needed

# function to convert varnames to workspace variables
function string_as_varname(s::AbstractString,v::Any)
    s=Symbol(s)
    @eval (($s) = ($v))
end

varnames = ["rmr", "em", "rmq", "rmt", "gfp", "rmy", "rmg", "et", "rmm",
            "zmm", "zmg", "zmy", "zmr", "zmq", "zmt", "yfp", "a", "mg", "mm",
            "q", "mt", "si", "r", "mq", "mr", "my"]

varidx = Dict{String, Int}()
for (idx, k) in enumerate(varnames)
    varidx[k] = idx
end

# parameters
# gsize = 4.6e6/2 # bp BID: 100269 (division by 2 -> each fork)
# reprate = gsize/40 #600*60 # nt/s -> nt/min BID: 109251

zeta = 1; 2.8312e+06; 1e0

thetar = 426.8693338968694
cl = 0
we = 4.139172187824451
kby = 1.0
vt = 726.0
s0 = 1.0e4
nx = 300.0
nq = 4
nr = 7549.0
ns = 0.5
wg = 0.0 #set_wg() 
k_cm = 0.005990373118888
gmax = 1260.0
thetax = 4.379733394834643
Km = 1.0e3
kbg = 1.0
Kq = 1.522190403737490e+05
Kp = 180.1378030928276
wr = 929.9678874564831
wq = 948.9349882947897
Kt = 1.0e3
wy = 0.0
kb = 1.0
Mref = 1.0e8
ku = 1.0
vm = 5800.0

thetar *= zeta
thetax *= zeta
Kp /= zeta

K_cd = 2.8 # fitted from Wallden et al., arXiv 2015

C = 40
D = 20
Vcrit = 1.2e8

# define rate constants
b = 0
dm = 0.1
rates = [b, dm]

# define initial conditions
fpos = []  # fork positions, here initially none

using MAT #https://github.com/JuliaIO/MAT.jl
init4expl = matread("init4expl.mat") #reads .mat file into a dict
gfp = 0; yfp = 0; mg = 0; my = 0; rmg = 0; rmy = 0
init = zeros(length(varnames))

for k in varnames
    string_as_varname(k, init4expl[k])

end

for (idx, k) in enumerate(varnames)
    eval(Meta.parse("init[$idx] = $k"))
end

# init= [rmr em rmq rmt gfp rmr rmg et rmm zmm zmg zmy zmr zmq zmt yfp a mg...
#     mm q mt si r mq mr my fpos]
