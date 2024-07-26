# model definition file

# host coupled with a single heterologous gene
"""
    ODE_model!(du,u,p,t)

Defines ODEs for extension of the growth model by Weisse et al, 2015, 
with effects of translation initiation on the expression of a heterologous protein.

Created by Evangelos-Marios Nikolados.
"""
function HETER_ODE_model!(du,u,p,t)

	thetar= p[1] # ribosomal threshold amount of energy at which transcription is half maximal 
	k_cm= p[2] # chloramphenicol-binding rate no idea why it is here
	s0= p[3] # external nutrient
	gmax= p[4] # maximal rate of translational elongation 
	cl= p[5] #
	thetax= p[6] # non-ribosomal threshold amount of energy at which transcription is half maximal
	Kt= p[7] # nutrient import threshold
	M= p[8] # fixed size of the cell i.e total cell mass
	we= p[9] # max. enzyme transcription rate we = wt = wm
	Km= p[10] # enzymatic threshold 
    f= cl*k_cm

	vm= p[11] # max. enzymatic rate
	nx= p[12] # length of non-ribosomal proteins
	Kq= p[13] # q-autoinhibition threshold (some constant for translation rate of q)
	Kp= p[14] # kgamma = gmax/Kp = 7
	vt= p[15] # max. nutrient import rate (something to do with nuimpt)
	wr= p[16] # maximal transcription rate of r
	wq= p[17] # maximal transcription rate of q 
	wp= p[18] # maximal transcription rate of p
	nq= p[19] # q-autoinhibition Hill coef (?)
	nr= p[20] # length of r (rybosome) protein
	ns= p[21] # one molecule of s yields ns molecules of a (energy) (nutrient efficiency)

    b= 0  # Hack - this variables also should be passed via p rather than defined here but for some reason it didn't work
    dm= 0.1
    kb= 0.0095
    ku= 1
	
    dmrep = p[22]
    dprep = p[23]
    kappa_ini = p[24]
	wmaxrep = p[25]
	kbrep= p[26]
	kurep= p[27]
	
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
    
    # rep is heterologous protein
	rep= u[15]
	mrep = u[16]
	rmrep = u[17]

	Kgamma= gmax/Kp # threshold amount of energy where elongation is halfmaximal
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation

	ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep)*gamma # term in expression of lambda
	vrep = gamma/nx*(kappa_ini*rmrep);
	wrep = wmaxrep*a/(thetax+a);

	lam= ttrate/M # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of si by em
 
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

	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep*rmrep+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep*r*mrep-lam*r 
	
	du[14]= ns*nucat-ttrate-lam*a # (2) energy is created from metabolism and lost from translation and dilution

	# heterologous equations, equations above also required modifications for heterologous case!
	du[15]= vrep - (lam + dprep)*rep
	du[16]= wrep - (lam + dmrep)*mrep + vrep - kbrep*r*mrep + (1-kappa_ini)*kurep*rmrep
	du[17]= -lam*rmrep + kbrep*r*mrep - (1-kappa_ini)*kurep*rmrep - vrep
	
end

# repressilator
function REPR_ODE_model!(du,u,p,t)

	thetar= p[1] # ribosomal threshold amount of energy at which transcription is half maximal 
	k_cm= p[2] # chloramphenicol-binding rate no idea why it is here
	s0= p[3] # external nutrient
	gmax= p[4] # maximal rate of translational elongation 
	cl= p[5] #
	thetax= p[6] # non-ribosomal threshold amount of energy at which transcription is half maximal
	Kt= p[7] # nutrient import threshold
	M= p[8] # fixed size of the cell i.e total cell mass
	we= p[9] # max. enzyme transcription rate we = wt = wm
	Km= p[10] # enzymatic threshold 
    f= cl*k_cm

	vm= p[11] # max. enzymatic rate
	nx= p[12] # length of non-ribosomal proteins
	Kq= p[13] # q-autoinhibition threshold (some constant for translation rate of q)
	Kp= p[14] # kgamma = gmax/Kp = 7
	vt= p[15] # max. nutrient import rate (something to do with nuimpt)
	wr= p[16] # maximal transcription rate of r
	wq= p[17] # maximal transcription rate of q 
	wp= p[18] # maximal transcription rate of p
	nq= p[19] # q-autoinhibition Hill coef (?)
	nr= p[20] # length of r (rybosome) protein
	ns= p[21] # one molecule of s yields ns molecules of a (energy) (nutrient efficiency)

    b= 0  # Hack - this variables also should be passed via p rather than defined here but for some reason it didn't work
    dm= 0.1
    kb= 0.0095
    ku= 1
	
    dmrep = p[22]
    dprep = p[23]
    kappa_ini = p[24]

	wmaxrep_1 = p[25]
	kbrep_1= p[26]
	kurep_1= p[27]	

	wmaxrep_2 = p[28]
	kbrep_2= p[29]
	kurep_2= p[30]	

	wmaxrep_3 = p[31]
	kbrep_3= p[32]
	kurep_3= p[33]	

    Kq_rep_1 = p[34]
    nq_rep_1 = p[35]

    Kq_rep_2 = p[36]
    nq_rep_2 = p[37]   

    Kq_rep_3 = p[38]
    nq_rep_3 = p[39]

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
	mr= u[12] # mRNA of r (ribosomes)
	r= u[13] # number of ribosomes
	a= u[14] # energy
    
	rep_1= u[15]
	mrep_1 = u[16]
	rmrep_1 = u[17]

	rep_2= u[18]
	mrep_2 = u[19]
	rmrep_2 = u[20]
		
	rep_3= u[21]
	mrep_3 = u[22]
	rmrep_3 = u[23]		

	Kgamma= gmax/Kp # threshold amount of energy where elongation is halfmaximal
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation

	ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep_1 + kappa_ini*rmrep_2 + kappa_ini*rmrep_3)*gamma
				
	vrep_1 = gamma/nx*(kappa_ini*rmrep_1);
	wrep_1 = wmaxrep_1*a/(thetax+a);

	vrep_2 = gamma/nx*(kappa_ini*rmrep_2);
	wrep_2 = wmaxrep_2*a/(thetax+a);

	vrep_3 = gamma/nx*(kappa_ini*rmrep_3);
	wrep_3 = wmaxrep_3*a/(thetax+a);	

	lam= ttrate/M # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of si by em
 
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

	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep_1*rmrep_1 + (1-kappa_ini)*kurep_2*rmrep_2 + (1-kappa_ini)*kurep_3*rmrep_3 +gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep_1 +vrep_2 +vrep_3 -kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep_1*r*mrep_1 -kbrep_2*r*mrep_2 -kbrep_3*r*mrep_3 -lam*r

	du[14]= ns*nucat-ttrate-lam*a # (2) energy is created from metabolism and lost from translation and dilution

	du[15]= vrep_1 - (lam + dprep)*rep_1
	du[16]= wrep_1 * (1 / (1 + (rep_3/Kq_rep_3)^nq_rep_3)) - (lam + dmrep)*mrep_1 + vrep_1 - kbrep_1*r*mrep_1 + (1-kappa_ini)*kurep_1*rmrep_1
	du[17]= -lam*rmrep_1 + kbrep_1*r*mrep_1 - (1-kappa_ini)*kurep_1*rmrep_1 - vrep_1

	du[18]= vrep_2 - (lam + dprep)*rep_2
	du[19]= wrep_2 * (1 / (1 + (rep_1/Kq_rep_1)^nq_rep_1)) - (lam + dmrep)*mrep_2 + vrep_2 - kbrep_2*r*mrep_2 + (1-kappa_ini)*kurep_2*rmrep_2
	du[20]= -lam*rmrep_2 + kbrep_2*r*mrep_2 - (1-kappa_ini)*kurep_2*rmrep_2 - vrep_2

	du[21]= vrep_3 - (lam + dprep)*rep_3
	du[22]= wrep_3 * (1 / (1 + (rep_2/Kq_rep_2)^nq_rep_2)) - (lam + dmrep)*mrep_3 + vrep_3 - kbrep_3*r*mrep_3 + (1-kappa_ini)*kurep_3*rmrep_3
	du[23]= -lam*rmrep_3 + kbrep_3*r*mrep_3 - (1-kappa_ini)*kurep_3*rmrep_3 - vrep_3

end

# NOT gate
function NOT_gate_ODE_model!(du,u,p,t)

	thetar= p[1] # ribosomal threshold amount of energy at which transcription is half maximal 
	k_cm= p[2] # chloramphenicol-binding rate no idea why it is here
	s0= p[3] # external nutrient
	gmax= p[4] # maximal rate of translational elongation 
	cl= p[5] #
	thetax= p[6] # non-ribosomal threshold amount of energy at which transcription is half maximal
	Kt= p[7] # nutrient import threshold
	M= p[8] # fixed size of the cell i.e total cell mass
	we= p[9] # max. enzyme transcription rate we = wt = wm
	Km= p[10] # enzymatic threshold 
    f= cl*k_cm

	vm= p[11] # max. enzymatic rate
	nx= p[12] # length of non-ribosomal proteins
	Kq= p[13] # q-autoinhibition threshold (some constant for translation rate of q)
	Kp= p[14] # kgamma = gmax/Kp = 7
	vt= p[15] # max. nutrient import rate (something to do with nuimpt)
	wr= p[16] # maximal transcription rate of r
	wq= p[17] # maximal transcription rate of q 
	wp= p[18] # maximal transcription rate of p
	nq= p[19] # q-autoinhibition Hill coef (?)
	nr= p[20] # length of r (rybosome) protein
	ns= p[21] # one molecule of s yields ns molecules of a (energy) (nutrient efficiency)

    b= 0  # Hack - this variables also should be passed via p rather than defined here but for some reason it didn't work
    dm= 0.1
    kb= 0.0095
    ku= 1
	
    dmrep = p[22]
    dprep = p[23]
    kappa_ini = p[24]

	wmaxrep_1 = p[25]
	kbrep_1= p[26]
	kurep_1= p[27]	

	wmaxrep_2 = p[28]
	kbrep_2= p[29]
	kurep_2= p[30]
    
    Kq_rep_1 = p[31]
    nq_rep_1 = p[32]

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

	# NOT gate
	rep_1= u[15]
	mrep_1 = u[16]
	rmrep_1 = u[17]

	rep_2= u[18]
	mrep_2 = u[19]
	rmrep_2 = u[20]

	Kgamma= gmax/Kp # threshold amount of energy where elongation is halfmaximal
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation

    ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep_1 + kappa_ini*rmrep_2)*gamma

	vrep_1 = gamma/nx*(kappa_ini*rmrep_1);
	wrep_1 = wmaxrep_1*a/(thetax+a);

	vrep_2 = gamma/nx*(kappa_ini*rmrep_2);
	wrep_2 = wmaxrep_2*a/(thetax+a);

	lam= ttrate/M # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of si by em

 
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

	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep_1*rmrep_1 + (1-kappa_ini)*kurep_2*rmrep_2+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep_1 +vrep_2-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep_1*r*mrep_1 -kbrep_2*r*mrep_2-lam*r
	
	du[14]= ns*nucat-ttrate-lam*a # (2) energy is created from metabolism and lost from translation and dilution

	# ODEs to implement a NOT gate :-> contains two genes
	du[15]= vrep_1 - (lam + dprep)*rep_1
	du[16]= wrep_1 - (lam + dmrep)*mrep_1 + vrep_1 - kbrep_1*r*mrep_1 + (1-kappa_ini)*kurep_1*rmrep_1
	du[17]= -lam*rmrep_1 + kbrep_1*r*mrep_1 - (1-kappa_ini)*kurep_1*rmrep_1 - vrep_1

	du[18]= vrep_2 - (lam + dprep)*rep_2
	du[19]= wrep_2 * (1 / (1 + (rep_1/Kq_rep_1)^nq_rep_1)) - (lam + dmrep)*mrep_2 + vrep_2 - kbrep_2*r*mrep_2 + (1-kappa_ini)*kurep_2*rmrep_2
	du[20]= -lam*rmrep_2 + kbrep_2*r*mrep_2 - (1-kappa_ini)*kurep_2*rmrep_2 - vrep_2

end

# AND gate
function AND_gate_ODE_model!(du,u,p,t)

	thetar= p[1] # ribosomal threshold amount of energy at which transcription is half maximal 
	k_cm= p[2] # chloramphenicol-binding rate no idea why it is here
	s0= p[3] # external nutrient
	gmax= p[4] # maximal rate of translational elongation 
	cl= p[5] #
	thetax= p[6] # non-ribosomal threshold amount of energy at which transcription is half maximal
	Kt= p[7] # nutrient import threshold
	M= p[8] # fixed size of the cell i.e total cell mass
	we= p[9] # max. enzyme transcription rate we = wt = wm
	Km= p[10] # enzymatic threshold 
    f= cl*k_cm

	vm= p[11] # max. enzymatic rate
	nx= p[12] # length of non-ribosomal proteins
	Kq= p[13] # q-autoinhibition threshold (some constant for translation rate of q)
	Kp= p[14] # kgamma = gmax/Kp = 7
	vt= p[15] # max. nutrient import rate (something to do with nuimpt)
	wr= p[16] # maximal transcription rate of r
	wq= p[17] # maximal transcription rate of q 
	wp= p[18] # maximal transcription rate of p
	nq= p[19] # q-autoinhibition Hill coef (?)
	nr= p[20] # length of r (rybosome) protein
	ns= p[21] # one molecule of s yields ns molecules of a (energy) (nutrient efficiency)

    b= 0  # Hack - this variables also should be passed via p rather than defined here but for some reason it didn't work
    dm= 0.1
    kb= 0.0095
    ku= 1
	
    dmrep = p[22]
    dprep = p[23]
    kappa_ini = p[24]

	wmaxrep_1 = p[25]
	kbrep_1= p[26]
	kurep_1= p[27]	

	wmaxrep_2 = p[28]
	kbrep_2= p[29]
	kurep_2= p[30]	

	wmaxrep_3 = p[31]
	kbrep_3= p[32]
	kurep_3= p[33]	

    Kq_rep_1 = p[34]
    nq_rep_1 = p[35]

    Kq_rep_2 = p[36]
    nq_rep_2 = p[37]   

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
	mr= u[12] # mRNA of r (ribosomes)
	r= u[13] # number of ribosomes
	a= u[14] # energy
    
	rep_1= u[15]
	mrep_1 = u[16]
	rmrep_1 = u[17]

	rep_2= u[18]
	mrep_2 = u[19]
	rmrep_2 = u[20]
		
	rep_3= u[21]
	mrep_3 = u[22]
	rmrep_3 = u[23]		

	Kgamma= gmax/Kp # threshold amount of energy where elongation is halfmaximal
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation

	ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep_1 + kappa_ini*rmrep_2 + kappa_ini*rmrep_3)*gamma
				
	vrep_1 = gamma/nx*(kappa_ini*rmrep_1);
	wrep_1 = wmaxrep_1*a/(thetax+a);

	vrep_2 = gamma/nx*(kappa_ini*rmrep_2);
	wrep_2 = wmaxrep_2*a/(thetax+a);

	vrep_3 = gamma/nx*(kappa_ini*rmrep_3);
	wrep_3 = wmaxrep_3*a/(thetax+a);	

	lam= ttrate/M # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of si by em
 
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

	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep_1*rmrep_1 + (1-kappa_ini)*kurep_2*rmrep_2 + (1-kappa_ini)*kurep_3*rmrep_3 +gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep_1 +vrep_2 +vrep_3 -kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep_1*r*mrep_1 -kbrep_2*r*mrep_2 -kbrep_3*r*mrep_3 -lam*r

	du[14]= ns*nucat-ttrate-lam*a # (2) energy is created from metabolism and lost from translation and dilution

	du[15]= vrep_1 - (lam + dprep)*rep_1
	du[16]= wrep_1 - (lam + dmrep)*mrep_1 + vrep_1 - kbrep_1*r*mrep_1 + (1-kappa_ini)*kurep_1*rmrep_1
	du[17]= -lam*rmrep_1 + kbrep_1*r*mrep_1 - (1-kappa_ini)*kurep_1*rmrep_1 - vrep_1

	du[18]= vrep_2 - (lam + dprep)*rep_2
	du[19]= wrep_2 - (lam + dmrep)*mrep_2 + vrep_2 - kbrep_2*r*mrep_2 + (1-kappa_ini)*kurep_2*rmrep_2
	du[20]= -lam*rmrep_2 + kbrep_2*r*mrep_2 - (1-kappa_ini)*kurep_2*rmrep_2 - vrep_2

	du[21]= vrep_3 - (lam + dprep)*rep_3
	du[22]= wrep_3 * ((((rep_1/Kq_rep_1)^nq_rep_1) / (1 + (rep_1/Kq_rep_1)^nq_rep_1)) * (((rep_2/Kq_rep_2)^nq_rep_2) / (1 + (rep_2/Kq_rep_2)^nq_rep_2))) - (lam + dmrep)*mrep_3 + vrep_3 - kbrep_3*r*mrep_3 + (1-kappa_ini)*kurep_3*rmrep_3
	du[23]= -lam*rmrep_3 + kbrep_3*r*mrep_3 - (1-kappa_ini)*kurep_3*rmrep_3 - vrep_3

end

# NAND gate
function NAND_gate_ODE_model!(du,u,p,t)

	thetar= p[1] # ribosomal threshold amount of energy at which transcription is half maximal 
	k_cm= p[2] # chloramphenicol-binding rate no idea why it is here
	s0= p[3] # external nutrient
	gmax= p[4] # maximal rate of translational elongation 
	cl= p[5] #
	thetax= p[6] # non-ribosomal threshold amount of energy at which transcription is half maximal
	Kt= p[7] # nutrient import threshold
	M= p[8] # fixed size of the cell i.e total cell mass
	we= p[9] # max. enzyme transcription rate we = wt = wm
	Km= p[10] # enzymatic threshold 
    f= cl*k_cm

	vm= p[11] # max. enzymatic rate
	nx= p[12] # length of non-ribosomal proteins
	Kq= p[13] # q-autoinhibition threshold (some constant for translation rate of q)
	Kp= p[14] # kgamma = gmax/Kp = 7
	vt= p[15] # max. nutrient import rate (something to do with nuimpt)
	wr= p[16] # maximal transcription rate of r
	wq= p[17] # maximal transcription rate of q 
	wp= p[18] # maximal transcription rate of p
	nq= p[19] # q-autoinhibition Hill coef (?)
	nr= p[20] # length of r (rybosome) protein
	ns= p[21] # one molecule of s yields ns molecules of a (energy) (nutrient efficiency)

    b= 0  # Hack - this variables also should be passed via p rather than defined here but for some reason it didn't work
    dm= 0.1
    kb= 0.0095
    ku= 1
	
    dmrep = p[22]
    dprep = p[23]
    kappa_ini = p[24]

	wmaxrep_1 = p[25]
	kbrep_1= p[26]
	kurep_1= p[27]	

	wmaxrep_2 = p[28]
	kbrep_2= p[29]
	kurep_2= p[30]	

	wmaxrep_3 = p[31]
	kbrep_3= p[32]
	kurep_3= p[33]	

	wmaxrep_4 = p[34]
	kbrep_4= p[35]
	kurep_4= p[36]	

    Kq_rep_1 = p[37]
    nq_rep_1 = p[38]

    Kq_rep_2 = p[39]
    nq_rep_2 = p[40]   

    Kq_rep_3 = p[41]
    nq_rep_3 = p[42]

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
	mr= u[12] # mRNA of r (ribosomes)
	r= u[13] # number of ribosomes
	a= u[14] # energy
    
	rep_1= u[15]
	mrep_1 = u[16]
	rmrep_1 = u[17]

	rep_2= u[18]
	mrep_2 = u[19]
	rmrep_2 = u[20]
		
	rep_3= u[21]
	mrep_3 = u[22]
	rmrep_3 = u[23]	
		
	rep_4= u[24]
	mrep_4 = u[25]
	rmrep_4 = u[26]	


	Kgamma= gmax/Kp # threshold amount of energy where elongation is halfmaximal
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation

	ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep_1 + kappa_ini*rmrep_2 + kappa_ini*rmrep_3 + kappa_ini*rmrep_4)*gamma
				
	vrep_1 = gamma/nx*(kappa_ini*rmrep_1);
	wrep_1 = wmaxrep_1*a/(thetax+a);

	vrep_2 = gamma/nx*(kappa_ini*rmrep_2);
	wrep_2 = wmaxrep_2*a/(thetax+a);

	vrep_3 = gamma/nx*(kappa_ini*rmrep_3);
	wrep_3 = wmaxrep_3*a/(thetax+a);

	vrep_4 = gamma/nx*(kappa_ini*rmrep_4);
	wrep_4 = wmaxrep_4*a/(thetax+a);	

	lam= ttrate/M # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of si by em
 
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

	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep_1*rmrep_1 + (1-kappa_ini)*kurep_2*rmrep_2 + (1-kappa_ini)*kurep_3*rmrep_3 + (1-kappa_ini)*kurep_4*rmrep_4+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep_1 +vrep_2 +vrep_3 +vrep_4-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep_1*r*mrep_1 -kbrep_2*r*mrep_2 -kbrep_3*r*mrep_3 -kbrep_4*r*mrep_4-lam*r

	du[14]= ns*nucat-ttrate-lam*a # (2) energy is created from metabolism and lost from translation and dilution

	du[15]= vrep_1 - (lam + dprep)*rep_1
	du[16]= wrep_1 - (lam + dmrep)*mrep_1 + vrep_1 - kbrep_1*r*mrep_1 + (1-kappa_ini)*kurep_1*rmrep_1
	du[17]= -lam*rmrep_1 + kbrep_1*r*mrep_1 - (1-kappa_ini)*kurep_1*rmrep_1 - vrep_1

	du[18]= vrep_2 - (lam + dprep)*rep_2
	du[19]= wrep_2 - (lam + dmrep)*mrep_2 + vrep_2 - kbrep_2*r*mrep_2 + (1-kappa_ini)*kurep_2*rmrep_2
	du[20]= -lam*rmrep_2 + kbrep_2*r*mrep_2 - (1-kappa_ini)*kurep_2*rmrep_2 - vrep_2

	du[21]= vrep_3 - (lam + dprep)*rep_3
	du[22]= wrep_3 * ((((rep_1/Kq_rep_1)^nq_rep_1) / (1 + (rep_1/Kq_rep_1)^nq_rep_1)) * (((rep_2/Kq_rep_2)^nq_rep_2) / (1 + (rep_2/Kq_rep_2)^nq_rep_2))) - (lam + dmrep)*mrep_3 + vrep_3 - kbrep_3*r*mrep_3 + (1-kappa_ini)*kurep_3*rmrep_3
	du[23]= -lam*rmrep_3 + kbrep_3*r*mrep_3 - (1-kappa_ini)*kurep_3*rmrep_3 - vrep_3

	du[24]= vrep_4 - (lam + dprep)*rep_4
	du[25]= wrep_4 * (1 / (1 + (rep_3/Kq_rep_3)^nq_rep_3)) - (lam + dmrep)*mrep_4 + vrep_4 - kbrep_4*r*mrep_4 + (1-kappa_ini)*kurep_4*rmrep_4
	du[26]= -lam*rmrep_4 + kbrep_4*r*mrep_4 - (1-kappa_ini)*kurep_4*rmrep_4 - vrep_4
end