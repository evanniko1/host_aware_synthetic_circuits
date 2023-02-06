# parameter values
thetar= 426.8693338968694
k_cm= 0.005990373118888
s0= 10000
gmax= 1260.0
cl= 0
thetax= 4.379733394834643
Kt= 1.0e3
M= 1.0e8
we= 4.139172187824451
Km= 1000
vm= 5800.0
nx= 300.0
Kq= 1.522190403737490e+05
Kp= 180.1378030928276
vt= 726.0
wr= 929.9678874564831
wq= 948.9349882947897
wp= 0.0
nq= 4
nr= 7549.0
dmrep= log(2)/2
dprep= log(2)/4
Kgamma= gmax/Kp

# base model steady state values
# vector's three last values are protein, mRNA, and ribosomal complex ss values respectively
u0 = [807.530517658162,7066.62814403594,2316.16746788774,69.1538861165317,7066.62758622631,69.1538892849256,9.64096853393699,9.64097064704617,236681.293193619,128.404551112062,322.904581569518,5.94267303259607,17.3282796528522,9.17329183561663,0,0,0]

# organize parameters in a vector
# TO-DO: change to dict-based implementation

p= [thetar, k_cm, s0, gmax, cl, 
	thetax, Kt, M, we, Km, 
	vm, nx, Kq, Kp, vt, 
	wr, wq, wp, nq, nr]