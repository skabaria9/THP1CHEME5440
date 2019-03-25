#Problem 2
    #-----------------------------------------------------------------
    #Set-up of Parameters
    #-----------------------------------------------------------------

    #Input the known values from problem statement
    Lx_char = 1000; #nt
    Lt_char = 330; #AA, 333 before
    Lx_1 = 1000; #1200; #nt
    Lx_2 = 1400; #2400; #nt
    Lx_3 = 1000; #600; #nt
    Ll_1 = Lx_1/3; #AA
    Ll_2 = Lx_2/3; #AA
    Ll_3 = Lx_3/3; #AA
    doubletime = 0.66*60; #minutes
    dil_rate = .693/doubletime; #min^-1
    mass_cell_water = 0.7 #percent
    copies_per_cell = 200; #200 copies
    total_cell_mass = 2.8e-13; #g/cell
    dryweight_cell = total_cell_mass*(1-mass_cell_water)

    #Degradation rate of mRNA (kdeg_m)
    halflife_m =  2.1; #min
    #convert to degradation rate, assuming first order kinetics
    kdeg_m = .693/halflife_m #min-1
    #Degradation rate of protein (kdeg_p), bionumbers
    halflife_p = 24*60; #min
    kdeg_p = .693/halflife_p #min-1

    #Concentration of gene (G) #nmol/gDW
    mol_gene = copies_per_cell / (6.02e23) #mol
    nmol_gene = mol_gene*10.0^9;
    G_gDW = nmol_gene/dryweight_cell #nmol/gDW
    #G_mM = mol_gene*1000/(vol_cell*mass_cell_water) #mM

    #kej for each mi
    e_X = 60.0*60; # nt/min, rate #42 before
    ke_mchar = e_X/Lx_char #min-1
    ke_m1 = ke_mchar * Lx_char/Lx_1 #min^-1
    ke_m2 = ke_mchar * Lx_char/Lx_2 #min^-1
    ke_m3 = ke_mchar * Lx_char/Lx_3 #min^-1
    e_L = 16.5*60; #AA/min #14.5 before
    ke_pchar = e_L/Lt_char #min^-1
    ke_p1 = ke_pchar * Lt_char/Ll_1
    ke_p2 = ke_pchar * Lt_char/Ll_2
    ke_p3 = ke_pchar * Lt_char/Ll_3

    #R_X and R_L
    RNAP_per_cell = 1150; #number per cell
    RNAP_per_cell_nmol = RNAP_per_cell/(6.02e23)*10.0^9; #nmol
    R_X_gDW = RNAP_per_cell_nmol/dryweight_cell #nmol/gDW
    rib_per_cell = 45000; #ribosomes/cell
    rib_per_cell_nmol = rib_per_cell/(6.02e23)*10.0^9 #nmol/cell
    R_L_gDW = rib_per_cell_nmol/dryweight_cell #nmol/gDW

    #Saturation Constants
    K_X = 0.24 #nmol/gDW
    K_T = 454.64; #nmol/gDW

    #set up tau_m (for mRNA) and tau_p (for protein)
    tau_m1 = 2.7 #ke_m1/kI_m
    tau_m2 = 2.7 #ke_m2/kI_m
    tau_m3 = 2.7 #ke_m3/kI_m
    tau_p1 = 0.8 #ke_p1/kI_p
    tau_p2 = 0.8 #ke_p2/kI_p
    tau_p3 = 0.8 #ke_p3/kI_p

    # #-----------------------------------------------------------------
    # #Set-up Inducer Concentration
    # #-----------------------------------------------------------------
    #
    # #Set up time dependence on I
    # if t<310 #before 60 min
    #    I = 0.0 #mM
    # else #after 60 min
    #     I = 0.0 #mM
    # end
    #
    # #Make I constant (for troublshooting)
    # #I=10

    #-----------------------------------------------------------------
    #Set-up promoter control rate parameters
    #-----------------------------------------------------------------
    #Promoter weights
    #PWeights -Inducer
    W_I1 = 100;
    #Weights - Protein P1
    W_11 = 0.000001;
    W_12 = 10.0;
    W_13 = 5.0;
    #Weights - Protein P2
    W_22 = 0.000001;
    W_23 = 50.0;
    #Weights - Protein 3 P3
    W_33 = 0.000001;

    #Promoter Binding Parameters
    #Effector=Inducer, Target=P1
    k_I1 = 0.30; #mM
    n_I1 = 1.5;
    #Effector=P1, Target=P2
    k_12 = 1000.0; #nmol/gDW
    n_12 = 1.5;
    #Effector=P1, Target=P3
    k_13 = 1000.0; #nmol/gDW
    n_13 = 1.5;
    #Effector=P2, Target=P3
    k_23 = 10000.0; #nmol/gDW
    n_23 = 10.0;
