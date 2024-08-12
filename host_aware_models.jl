# model definition file

###############
# host coupled with a single heterologous gene
"""
    HETER_ODE_model!(du,u,p,t)

Defines ODEs for extension of the growth model by Weisse et al, 2015, 
with effects of translation initiation on the expression of a heterologous protein.

Created by Evangelos-Marios Nikolados.
"""
function HETER_ODE_model!(du,u,p,t)

	thetar= p[1]  # ribosome transcription threshold
	k_cm= p[2]    # chloramphenicol-binding rate
	s0= p[3]      # external nutrient
	gmax= p[4]    # maximal translation elongation rate
	cl= p[5]      # 
	thetax= p[6]  # non-ribosomal transcription threshold
	Kt= p[7]      # nutrient import threshold
	M= p[8]       # total cell mass
	we= p[9]      # max. enzyme transcription rate we = wt = wm
	Km= p[10]     # enzymatic threshold 
    #f= cl*k_cm
	vm= p[11]     # max. enzymatic rate
	nx= p[12]     # length of non-ribosomal proteins
	Kq= p[13]     # q-autoinhibition threshold
	Kp= p[14]     # 
	vt= p[15]     # max. nutrient import rate (something to do with nuimpt)
	wr= p[16]     # maximal transcription rate of r
	wq= p[17]     # maximal transcription rate of q 
	wp= p[18]     # maximal transcription rate of p
	nq= p[19]     # q-autoinhibition Hill coefficient
	nr= p[20]     # ribosome length
	ns= p[21]     # nutrient efficiency

	dm= 0.1        # endogenous mRNA degradation rate
    kb= 0.0095     # endogenous mRNA-ribosome binding rate
    ku= 1          # endogenous mRNA-ribosome unbinding rate
	
    dmrep = p[22]  # heterologous mRNA degradation rate
    dprep = p[23]  # active protein degradation rate
    kappa_ini = p[24] # translation initiation efficiency
	wmaxrep = p[25]   # maximal transcription rate for heterologous gene
	kbrep= p[26]	  # heterologous mRNA-ribosome binding rate
	kurep= p[27]      # heterologous mRNA-ribosome unbinding rate
	
	rmr= u[1]   # ribosomal mRNA :: ribosome complex
	em= u[2]    # metabolic enzyme
	rmq= u[3]   # housekeeping protein mRNA :: ribosome complex
	rmt= u[4]   # transporter protein mRNA :: ribosome complex
	et= u[5]    # transporter protein
	rmm= u[6]   # metabolic enzyme mRNA :: ribosome complex
	mt= u[7]    # transporter protein mRNA
	mm= u[8]    # metabolic enzyme mRNA
	q= u[9]     # house keeping proteins
	si= u[10]   # internalized nutrient
	mq= u[11]   # housekeeping protein mRNA
	mr= u[12]   # ribosomal mRNA
	r= u[13]    # free ribosomes
	a= u[14]    # energy
    
    # rep is heterologous protein
	rep= u[15]    # reporter heterologous protein
	mrep = u[16]  # reporter heterologous protein mRNA
	rmrep = u[17] # reporter heterologous protein mRNA :: ribosome complex

	Kgamma= gmax/Kp            # translation elongation threshold
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation

	ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep)*gamma # total translation rate
	vrep = gamma/nx*(kappa_ini*rmrep);                      # heterologous protein translation rate
	wrep = wmaxrep*a/(thetax+a);                            # heterologous transcription rate

	lam= ttrate/M              # growth rate
	nucat= em*vm*si/(Km + si)  # rate of metabolism of internalized nutrient
 
	# Ordinary Differential Equations
	du[1]= kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr
	du[2]= gamma/nx*rmm-lam*em
	du[3]= kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq
	du[4]= kb*r*mt+-ku*rmt-gamma/nx*rmt-lam*rmt
	du[5]= gamma/nx*rmt-lam*et
	du[6]= kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm
	du[7]= (we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt
	du[8]= (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm
	du[9]= gamma/nx*rmq-lam*q
	du[10]= (et*vt*s0/(Kt + s0))-nucat-lam*si
	du[11]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq
	du[12]= (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr
	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep*rmrep+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep*r*mrep-lam*r 
	du[14]= ns*nucat-ttrate-lam*a
	du[15]= vrep - (lam + dprep)*rep
	du[16]= wrep - (lam + dmrep)*mrep + vrep - kbrep*r*mrep + (1-kappa_ini)*kurep*rmrep
	du[17]= -lam*rmrep + kbrep*r*mrep - (1-kappa_ini)*kurep*rmrep - vrep
	
end

###############
# host coupled with a repressilator
"""
    REPR_ODE_model!(du,u,p,t)


Defines ODEs for extension of the growth model by Weisse et al, 2015, 
with effects of translation initiation on the expression of a repressilator synthetic circuit.

Created by Evangelos-Marios Nikolados.
"""
function REPR_ODE_model!(du,u,p,t)

	thetar= p[1]  # ribosome transcription threshold
	k_cm= p[2]    # chloramphenicol-binding rate
	s0= p[3]      # external nutrient
	gmax= p[4]    # maximal translation elongation rate
	cl= p[5]      # 
	thetax= p[6]  # non-ribosomal transcription threshold
	Kt= p[7]      # nutrient import threshold
	M= p[8]       # total cell mass
	we= p[9]      # max. enzyme transcription rate we = wt = wm
	Km= p[10]     # enzymatic threshold 
    #f= cl*k_cm
	vm= p[11]     # max. enzymatic rate
	nx= p[12]     # length of non-ribosomal proteins
	Kq= p[13]     # q-autoinhibition threshold
	Kp= p[14]     # 
	vt= p[15]     # max. nutrient import rate (something to do with nuimpt)
	wr= p[16]     # maximal transcription rate of r
	wq= p[17]     # maximal transcription rate of q 
	wp= p[18]     # maximal transcription rate of p
	nq= p[19]     # q-autoinhibition Hill coefficient
	nr= p[20]     # ribosome length
	ns= p[21]     # nutrient efficiency

	dm= 0.1        # endogenous mRNA degradation rate
    kb= 0.0095     # endogenous mRNA-ribosome binding rate
    ku= 1          # endogenous mRNA-ribosome unbinding rate
	
    dmrep = p[22]     # heterologous mRNA degradation rate
    dprep = p[23]     # active protein degradation rate
    kappa_ini = p[24] # translation initiation efficiency

	wmaxrep_1 = p[25] # maximal transcription rate for heterologous gene 1
	kbrep_1= p[26]    # heterologous mRNA-ribosome 1 binding rate
	kurep_1= p[27]	  # heterologous mRNA-ribosome 1 unbinding rate

	wmaxrep_2 = p[28] # maximal transcription rate for heterologous gene 2
	kbrep_2= p[29]    # heterologous mRNA-ribosome 2 binding rate
	kurep_2= p[30]	  # heterologous mRNA-ribosome 2 unbinding rate

	wmaxrep_3 = p[31] # maximal transcription rate for heterologous gene 3
	kbrep_3= p[32]    # heterologous mRNA-ribosome 3 binding rate
	kurep_3= p[33]	  # heterologous mRNA-ribosome 3 unbinding rate

    Kq_rep_1 = p[34]  # heterologous protein 1 inhibition threshold
    nq_rep_1 = p[35]  # heterologous protein 1 inhibition Hill coefficient

    Kq_rep_2 = p[36]  # heterologous protein 2 inhibition threshold
    nq_rep_2 = p[37]  # heterologous protein 2 inhibition Hill coefficient 

    Kq_rep_3 = p[38]  # heterologous protein 3 inhibition threshold
    nq_rep_3 = p[39]  # heterologous protein 3 inhibition Hill coefficient

	rmr= u[1]   # ribosomal mRNA :: ribosome complex
	em= u[2]    # metabolic enzyme
	rmq= u[3]   # housekeeping protein mRNA :: ribosome complex
	rmt= u[4]   # transporter protein mRNA :: ribosome complex
	et= u[5]    # transporter protein
	rmm= u[6]   # metabolic enzyme mRNA :: ribosome complex
	mt= u[7]    # transporter protein mRNA
	mm= u[8]    # metabolic enzyme mRNA
	q= u[9]     # house keeping proteins
	si= u[10]   # internalized nutrient
	mq= u[11]   # housekeeping protein mRNA
	mr= u[12]   # ribosomal mRNA
	r= u[13]    # free ribosomes
	a= u[14]    # energy

	# Repressilator
	rep_1= u[15]    # reporter heterologous protein 1
	mrep_1 = u[16]  # reporter heterologous protein 1 mRNA 
	rmrep_1 = u[17] # reporter heterologous protein 1 mRNA :: ribosome complex

	rep_2= u[18]    # reporter heterologous protein 2
	mrep_2 = u[19]  # reporter heterologous protein 2 mRNA 
	rmrep_2 = u[20] # reporter heterologous protein 2 mRNA :: ribosome complex
		
	rep_3= u[21]    # reporter heterologous protein 3
	mrep_3 = u[22]  # reporter heterologous protein 3 mRNA 
	rmrep_3 = u[23] # reporter heterologous protein 3 mRNA :: ribosome complex

	Kgamma= gmax/Kp            # translation elongation threshold
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation

	ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep_1 + kappa_ini*rmrep_2 + kappa_ini*rmrep_3)*gamma # total translation rate
				
	vrep_1 = gamma/nx*(kappa_ini*rmrep_1); # heterologous protein 1 translation rate
	wrep_1 = wmaxrep_1*a/(thetax+a);       # heterologous transcription rate 1

	vrep_2 = gamma/nx*(kappa_ini*rmrep_2); # heterologous protein 2 translation rate
	wrep_2 = wmaxrep_2*a/(thetax+a);       # heterologous transcription rate 2

	vrep_3 = gamma/nx*(kappa_ini*rmrep_3); # heterologous protein 3 translation rate
	wrep_3 = wmaxrep_3*a/(thetax+a);       # heterologous transcription rate 3

	lam= ttrate/M             # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of internalized nutrient
 
	# Ordinary Differential Equations
	du[1]= kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr
	du[2]= gamma/nx*rmm-lam*em
	du[3]= kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq
	du[4]= kb*r*mt-ku*rmt-gamma/nx*rmt-lam*rmt
	du[5]= gamma/nx*rmt-lam*et
	du[6]= kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm
	du[7]= (we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt
	du[8]= (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm
	du[9]= gamma/nx*rmq-lam*q
	du[10]= (et*vt*s0/(Kt + s0))-nucat-lam*si
	du[11]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq
	du[12]= (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr
	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep_1*rmrep_1 + (1-kappa_ini)*kurep_2*rmrep_2 + (1-kappa_ini)*kurep_3*rmrep_3 +gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep_1 +vrep_2 +vrep_3 -kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep_1*r*mrep_1 -kbrep_2*r*mrep_2 -kbrep_3*r*mrep_3 -lam*r
	du[14]= ns*nucat-ttrate-lam*a
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

###############
# host coupled with a NOT gate
"""
    NOT_gate_ODE_model!(du,u,p,t)

Defines ODEs for extension of the growth model by Weisse et al, 2015, 
with effects of translation initiation on the expression of a NOT gate synthetic circuit.

Created by Evangelos-Marios Nikolados.
"""
function NOT_gate_ODE_model!(du,u,p,t)

	thetar= p[1]  # ribosome transcription threshold
	k_cm= p[2]    # chloramphenicol-binding rate
	s0= p[3]      # external nutrient
	gmax= p[4]    # maximal translation elongation rate
	cl= p[5]      # 
	thetax= p[6]  # non-ribosomal transcription threshold
	Kt= p[7]      # nutrient import threshold
	M= p[8]       # total cell mass
	we= p[9]      # max. enzyme transcription rate we = wt = wm
	Km= p[10]     # enzymatic threshold 
    #f= cl*k_cm
	vm= p[11]     # max. enzymatic rate
	nx= p[12]     # length of non-ribosomal proteins
	Kq= p[13]     # q-autoinhibition threshold
	Kp= p[14]     # 
	vt= p[15]     # max. nutrient import rate (something to do with nuimpt)
	wr= p[16]     # maximal transcription rate of r
	wq= p[17]     # maximal transcription rate of q 
	wp= p[18]     # maximal transcription rate of p
	nq= p[19]     # q-autoinhibition Hill coefficient
	nr= p[20]     # ribosome length
	ns= p[21]     # nutrient efficiency

	dm= 0.1        # endogenous mRNA degradation rate
    kb= 0.0095     # endogenous mRNA-ribosome binding rate
    ku= 1          # endogenous mRNA-ribosome unbinding rate
	
    dmrep = p[22]     # heterologous mRNA degradation rate
    dprep = p[23]     # active protein degradation rate
    kappa_ini = p[24] # translation initiation efficiency

	wmaxrep_1 = p[25] # maximal transcription rate for heterologous gene 1
	kbrep_1= p[26]    # heterologous mRNA-ribosome 1 binding rate
	kurep_1= p[27]	  # heterologous mRNA-ribosome 1 unbinding rate

	wmaxrep_2 = p[28] # maximal transcription rate for heterologous gene 2
	kbrep_2= p[29]    # heterologous mRNA-ribosome 2 binding rate
	kurep_2= p[30]	  # heterologous mRNA-ribosome 2 unbinding rate
    
    Kq_rep_1 = p[31]  # heterologous protein 1 inhibition threshold
    nq_rep_1 = p[32]  # heterologous protein 1 inhibition Hill coefficient

	rmr= u[1]   # ribosomal mRNA :: ribosome complex
	em= u[2]    # metabolic enzyme
	rmq= u[3]   # housekeeping protein mRNA :: ribosome complex
	rmt= u[4]   # transporter protein mRNA :: ribosome complex
	et= u[5]    # transporter protein
	rmm= u[6]   # metabolic enzyme mRNA :: ribosome complex
	mt= u[7]    # transporter protein mRNA
	mm= u[8]    # metabolic enzyme mRNA
	q= u[9]     # house keeping proteins
	si= u[10]   # internalized nutrient
	mq= u[11]   # housekeeping protein mRNA
	mr= u[12]   # ribosomal mRNA
	r= u[13]    # free ribosomes
	a= u[14]    # energy

	# NOT gate
	rep_1= u[15]    # reporter heterologous protein 1
	mrep_1 = u[16]  # reporter heterologous protein 1 mRNA 
	rmrep_1 = u[17] # reporter heterologous protein 1 mRNA :: ribosome complex

	rep_2= u[18]    # reporter heterologous protein 2
	mrep_2 = u[19]  # reporter heterologous protein 2 mRNA 
	rmrep_2 = u[20] # reporter heterologous protein 2 mRNA :: ribosome complex

	Kgamma= gmax/Kp            # translation elongation threshold
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation

    ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep_1 + kappa_ini*rmrep_2)*gamma # total translation rate

	vrep_1 = gamma/nx*(kappa_ini*rmrep_1); # heterologous protein 1 translation rate
	wrep_1 = wmaxrep_1*a/(thetax+a);       # heterologous transcription rate 1

	vrep_2 = gamma/nx*(kappa_ini*rmrep_2); # heterologous protein 2 translation rate
	wrep_2 = wmaxrep_2*a/(thetax+a);       # heterologous transcription rate 2

	lam= ttrate/M             # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of internalized nutrient
 
	# Ordinary Differential Equations
	du[1]= kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr
	du[2]= gamma/nx*rmm-lam*em
	du[3]= kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq
	du[4]= kb*r*mt-ku*rmt-gamma/nx*rmt-lam*rmt
	du[5]= gamma/nx*rmt-lam*et
	du[6]= kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm
	du[7]= (we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt
	du[8]= (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm
	du[9]= gamma/nx*rmq-lam*q
	du[10]= (et*vt*s0/(Kt + s0))-nucat-lam*si
	du[11]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq
	du[12]= (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr
	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep_1*rmrep_1 + (1-kappa_ini)*kurep_2*rmrep_2+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep_1 +vrep_2-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep_1*r*mrep_1 -kbrep_2*r*mrep_2-lam*r
	du[14]= ns*nucat-ttrate-lam*a
	# ODEs to implement a NOT gate :-> contains two genes
	du[15]= vrep_1 - (lam + dprep)*rep_1
	du[16]= wrep_1 - (lam + dmrep)*mrep_1 + vrep_1 - kbrep_1*r*mrep_1 + (1-kappa_ini)*kurep_1*rmrep_1
	du[17]= -lam*rmrep_1 + kbrep_1*r*mrep_1 - (1-kappa_ini)*kurep_1*rmrep_1 - vrep_1
	du[18]= vrep_2 - (lam + dprep)*rep_2
	du[19]= wrep_2 * (1 / (1 + (rep_1/Kq_rep_1)^nq_rep_1)) - (lam + dmrep)*mrep_2 + vrep_2 - kbrep_2*r*mrep_2 + (1-kappa_ini)*kurep_2*rmrep_2
	du[20]= -lam*rmrep_2 + kbrep_2*r*mrep_2 - (1-kappa_ini)*kurep_2*rmrep_2 - vrep_2

end

###############
# host coupled with a AND gate
"""
    AND_gate_ODE_model!(du,u,p,t)

Defines ODEs for extension of the growth model by Weisse et al, 2015, 
with effects of translation initiation on the expression of a AND gate synthetic circuit.

Created by Evangelos-Marios Nikolados.
"""
function AND_gate_ODE_model!(du,u,p,t)

	thetar= p[1]  # ribosome transcription threshold
	k_cm= p[2]    # chloramphenicol-binding rate
	s0= p[3]      # external nutrient
	gmax= p[4]    # maximal translation elongation rate
	cl= p[5]      # 
	thetax= p[6]  # non-ribosomal transcription threshold
	Kt= p[7]      # nutrient import threshold
	M= p[8]       # total cell mass
	we= p[9]      # max. enzyme transcription rate we = wt = wm
	Km= p[10]     # enzymatic threshold 
    #f= cl*k_cm
	vm= p[11]     # max. enzymatic rate
	nx= p[12]     # length of non-ribosomal proteins
	Kq= p[13]     # q-autoinhibition threshold
	Kp= p[14]     # 
	vt= p[15]     # max. nutrient import rate (something to do with nuimpt)
	wr= p[16]     # maximal transcription rate of r
	wq= p[17]     # maximal transcription rate of q 
	wp= p[18]     # maximal transcription rate of p
	nq= p[19]     # q-autoinhibition Hill coefficient
	nr= p[20]     # ribosome length
	ns= p[21]     # nutrient efficiency

	dm= 0.1        # endogenous mRNA degradation rate
    kb= 0.0095     # endogenous mRNA-ribosome binding rate
    ku= 1          # endogenous mRNA-ribosome unbinding rate

    dmrep = p[22]     # heterologous mRNA degradation rate
    dprep = p[23]     # active protein degradation rate
    kappa_ini = p[24] # translation initiation efficiency

	wmaxrep_1 = p[25] # maximal transcription rate for heterologous gene 1
	kbrep_1= p[26]    # heterologous mRNA-ribosome 1 binding rate
	kurep_1= p[27]	  # heterologous mRNA-ribosome 1 unbinding rate

	wmaxrep_2 = p[28] # maximal transcription rate for heterologous gene 2
	kbrep_2= p[29]    # heterologous mRNA-ribosome 2 binding rate
	kurep_2= p[30]	  # heterologous mRNA-ribosome 2 unbinding rate

	wmaxrep_3 = p[31] # maximal transcription rate for heterologous gene 3
	kbrep_3= p[32]    # heterologous mRNA-ribosome 3 binding rate
	kurep_3= p[33]	  # heterologous mRNA-ribosome 3 unbinding rate

    Kq_rep_1 = p[34]  # heterologous protein 1 inhibition threshold
    nq_rep_1 = p[35]  # heterologous protein 1 inhibition Hill coefficient

    Kq_rep_2 = p[36]  # heterologous protein 2 inhibition threshold
    nq_rep_2 = p[37]  # heterologous protein 2 inhibition Hill coefficient 

	rmr= u[1]   # ribosomal mRNA :: ribosome complex
	em= u[2]    # metabolic enzyme
	rmq= u[3]   # housekeeping protein mRNA :: ribosome complex
	rmt= u[4]   # transporter protein mRNA :: ribosome complex
	et= u[5]    # transporter protein
	rmm= u[6]   # metabolic enzyme mRNA :: ribosome complex
	mt= u[7]    # transporter protein mRNA
	mm= u[8]    # metabolic enzyme mRNA
	q= u[9]     # house keeping proteins
	si= u[10]   # internalized nutrient
	mq= u[11]   # housekeeping protein mRNA
	mr= u[12]   # ribosomal mRNA
	r= u[13]    # free ribosomes
	a= u[14]    # energy

	# AND gate
	rep_1= u[15]    # reporter heterologous protein 1
	mrep_1 = u[16]  # reporter heterologous protein 1 mRNA 
	rmrep_1 = u[17] # reporter heterologous protein 1 mRNA :: ribosome complex

	rep_2= u[18]    # reporter heterologous protein 2
	mrep_2 = u[19]  # reporter heterologous protein 2 mRNA 
	rmrep_2 = u[20] # reporter heterologous protein 2 mRNA :: ribosome complex
		
	rep_3= u[21]    # reporter heterologous protein 3
	mrep_3 = u[22]  # reporter heterologous protein 3 mRNA 
	rmrep_3 = u[23] # reporter heterologous protein 3 mRNA :: ribosome complex

	Kgamma= gmax/Kp            # translation elongation threshold
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation

	ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep_1 + kappa_ini*rmrep_2 + kappa_ini*rmrep_3)*gamma # total translation rate
				
	vrep_1 = gamma/nx*(kappa_ini*rmrep_1); # heterologous protein 1 translation rate
	wrep_1 = wmaxrep_1*a/(thetax+a);       # heterologous transcription rate 1

	vrep_2 = gamma/nx*(kappa_ini*rmrep_2); # heterologous protein 2 translation rate
	wrep_2 = wmaxrep_2*a/(thetax+a);       # heterologous transcription rate 2

	vrep_3 = gamma/nx*(kappa_ini*rmrep_3); # heterologous protein 3 translation rate
	wrep_3 = wmaxrep_3*a/(thetax+a);       # heterologous transcription rate 3

	lam= ttrate/M             # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of internalized nutrient
 
	# Ordinary Differential Equations
	du[1]= kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr
	du[2]= gamma/nx*rmm-lam*em
	du[3]= kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq
	du[4]= kb*r*mt-ku*rmt-gamma/nx*rmt-lam*rmt
	du[5]= gamma/nx*rmt-lam*et
	du[6]= kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm
	du[7]= (we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt
	du[8]= (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm
	du[9]= gamma/nx*rmq-lam*q
	du[10]= (et*vt*s0/(Kt + s0))-nucat-lam*si
	du[11]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq
	du[12]= (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr
	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep_1*rmrep_1 + (1-kappa_ini)*kurep_2*rmrep_2 + (1-kappa_ini)*kurep_3*rmrep_3 +gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep_1 +vrep_2 +vrep_3 -kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep_1*r*mrep_1 -kbrep_2*r*mrep_2 -kbrep_3*r*mrep_3 -lam*r
	du[14]= ns*nucat-ttrate-lam*a
	# ODEs to implement a AND gate :-> contains three genes
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

###############
# host coupled with a NAND gate
"""
    NAND_gate_ODE_model!(du,u,p,t)

Defines ODEs for extension of the growth model by Weisse et al, 2015, 
with effects of translation initiation on the expression of a NAND gate synthetic circuit.

Created by Evangelos-Marios Nikolados.
"""
function NAND_gate_ODE_model!(du,u,p,t)

	thetar= p[1]  # ribosome transcription threshold
	k_cm= p[2]    # chloramphenicol-binding rate
	s0= p[3]      # external nutrient
	gmax= p[4]    # maximal translation elongation rate
	cl= p[5]      # 
	thetax= p[6]  # non-ribosomal transcription threshold
	Kt= p[7]      # nutrient import threshold
	M= p[8]       # total cell mass
	we= p[9]      # max. enzyme transcription rate we = wt = wm
	Km= p[10]     # enzymatic threshold 
    #f= cl*k_cm
	vm= p[11]     # max. enzymatic rate
	nx= p[12]     # length of non-ribosomal proteins
	Kq= p[13]     # q-autoinhibition threshold
	Kp= p[14]     # 
	vt= p[15]     # max. nutrient import rate (something to do with nuimpt)
	wr= p[16]     # maximal transcription rate of r
	wq= p[17]     # maximal transcription rate of q 
	wp= p[18]     # maximal transcription rate of p
	nq= p[19]     # q-autoinhibition Hill coefficient
	nr= p[20]     # ribosome length
	ns= p[21]     # nutrient efficiency

	dm= 0.1        # endogenous mRNA degradation rate
    kb= 0.0095     # endogenous mRNA-ribosome binding rate
    ku= 1          # endogenous mRNA-ribosome unbinding rate
	
    dmrep = p[22]     # heterologous mRNA degradation rate
    dprep = p[23]     # active protein degradation rate
    kappa_ini = p[24] # translation initiation efficiency

	wmaxrep_1 = p[25] # maximal transcription rate for heterologous gene 1
	kbrep_1= p[26]    # heterologous mRNA-ribosome 1 binding rate
	kurep_1= p[27]	  # heterologous mRNA-ribosome 1 unbinding rate

	wmaxrep_2 = p[28] # maximal transcription rate for heterologous gene 2
	kbrep_2= p[29]    # heterologous mRNA-ribosome 2 binding rate
	kurep_2= p[30]	  # heterologous mRNA-ribosome 2 unbinding rate

	wmaxrep_3 = p[31] # maximal transcription rate for heterologous gene 3
	kbrep_3= p[32]    # heterologous mRNA-ribosome 3 binding rate
	kurep_3= p[33]	  # heterologous mRNA-ribosome 3 unbinding rate

	wmaxrep_4 = p[34] # maximal transcription rate for heterologous gene 4
	kbrep_4= p[35]    # heterologous mRNA-ribosome 3 binding rate
	kurep_4= p[36]    # heterologous mRNA-ribosome 3 unbinding rate

	Kq_rep_1 = p[37]  # heterologous protein 1 inhibition threshold
    nq_rep_1 = p[38]  # heterologous protein 1 inhibition Hill coefficient

    Kq_rep_2 = p[39]  # heterologous protein 2 inhibition threshold
    nq_rep_2 = p[40]  # heterologous protein 2 inhibition Hill coefficient 

    Kq_rep_3 = p[41]  # heterologous protein 3 inhibition threshold
    nq_rep_3 = p[42]  # heterologous protein 3 inhibition Hill coefficient

	rmr= u[1]   # ribosomal mRNA :: ribosome complex
	em= u[2]    # metabolic enzyme
	rmq= u[3]   # housekeeping protein mRNA :: ribosome complex
	rmt= u[4]   # transporter protein mRNA :: ribosome complex
	et= u[5]    # transporter protein
	rmm= u[6]   # metabolic enzyme mRNA :: ribosome complex
	mt= u[7]    # transporter protein mRNA
	mm= u[8]    # metabolic enzyme mRNA
	q= u[9]     # house keeping proteins
	si= u[10]   # internalized nutrient
	mq= u[11]   # housekeeping protein mRNA
	mr= u[12]   # ribosomal mRNA
	r= u[13]    # free ribosomes
	a= u[14]    # energy

	# NAND gate
	rep_1= u[15]    # reporter heterologous protein 1
	mrep_1 = u[16]  # reporter heterologous protein 1 mRNA 
	rmrep_1 = u[17] # reporter heterologous protein 1 mRNA :: ribosome complex

	rep_2 = u[18]   # reporter heterologous protein 2
	mrep_2 = u[19]  # reporter heterologous protein 2 mRNA 
	rmrep_2 = u[20] # reporter heterologous protein 2 mRNA :: ribosome complex
		
	rep_3 = u[21]   # reporter heterologous protein 3
	mrep_3 = u[22]  # reporter heterologous protein 3 mRNA 
	rmrep_3 = u[23] # reporter heterologous protein 3 mRNA :: ribosome complex
		
	rep_4 = u[24]   # reporter heterologous protein 4
	mrep_4 = u[25]  # reporter heterologous protein 4 mRNA
	rmrep_4 = u[26] # reporter heterologous protein 4 mRNA :: ribosome complex	

	Kgamma= gmax/Kp            # translation elongation threshold
	gamma= gmax*a/(Kgamma + a) # effective rate of translational elongation

	ttrate= (rmq + rmr + rmt + rmm + kappa_ini*rmrep_1 + kappa_ini*rmrep_2 + kappa_ini*rmrep_3 + kappa_ini*rmrep_4)*gamma # total translation rate

	vrep_1 = gamma/nx*(kappa_ini*rmrep_1); # heterologous protein 1 translation rate
	wrep_1 = wmaxrep_1*a/(thetax+a);       # heterologous transcription rate 1

	vrep_2 = gamma/nx*(kappa_ini*rmrep_2); # heterologous protein 2 translation rate
	wrep_2 = wmaxrep_2*a/(thetax+a);       # heterologous transcription rate 2

	vrep_3 = gamma/nx*(kappa_ini*rmrep_3); # heterologous protein 3 translation rate
	wrep_3 = wmaxrep_3*a/(thetax+a);       # heterologous transcription rate 3

	vrep_4 = gamma/nx*(kappa_ini*rmrep_4); # heterologous protein 4 translation rate
	wrep_4 = wmaxrep_4*a/(thetax+a);       # heterologous transcription rate 4

	lam= ttrate/M             # growth rate
	nucat= em*vm*si/(Km + si) # rate of metabolism of internalized nutrient
 
	# Ordinary Differential Equations
	du[1]= kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr
	du[2]= gamma/nx*rmm-lam*em
	du[3]= kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq
	du[4]= kb*r*mt-ku*rmt-gamma/nx*rmt-lam*rmt
	du[5]= gamma/nx*rmt-lam*et
	du[6]= kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm
	du[7]= (we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt
	du[8]= (we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm
	du[9]= gamma/nx*rmq-lam*q
	du[10]= (et*vt*s0/(Kt + s0))-nucat-lam*si
	du[11]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq
	du[12]= (wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr
	du[13]= ku*rmr+ku*rmt+ku*rmm+ku*rmq+(1-kappa_ini)*kurep_1*rmrep_1 + (1-kappa_ini)*kurep_2*rmrep_2 + (1-kappa_ini)*kurep_3*rmrep_3 + (1-kappa_ini)*kurep_4*rmrep_4+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq+vrep_1 +vrep_2 +vrep_3 +vrep_4-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-kbrep_1*r*mrep_1 -kbrep_2*r*mrep_2 -kbrep_3*r*mrep_3 -kbrep_4*r*mrep_4-lam*r
	du[14]= ns*nucat-ttrate-lam*a
	# ODEs to implement a NAND gate :-> contains four genes
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