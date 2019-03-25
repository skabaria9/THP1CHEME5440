#Use this file to output the x_matrices when there is a parameter change

#Included packages
using Pkg
Pkg.add("LinearAlgebra")
using LinearAlgebra
#Included sub-parameters and functions
include("ps2_parameters.jl")

function parameter_change(factors)
    #factors = [factor_K_L,factor_K_X,factor_eX,factor_eL,factor_RX,factor_RL,factor_tauX,factor_tauL,
        #factor_wI1,factor_w12,factor_w13,factor_w23,factor_wjj]

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

    factor_K_L = factors[1]
    factor_K_X = factors[2]
    factor_eX = factors[3]
    factor_eL = factors[4]
    factor_RX = factors[5]
    factor_RL = factors[6]
    factor_tauX = factors[7]
    factor_tauL = factors[8]
    factor_wI1 = factors[9]
    factor_w12 = factors[10]
    factor_w13 = factors[11]
    factor_w23 = factors[12]
    factor_wjj = factors[13]
    factor_k_I1 = factors[14]
    factor_n_I1 = factors[15]
    factor_k_12 = factors[16]
    factor_n_12 = factors[17]
    factor_k_13 = factors[18]
    factor_n_13 = factors[19]
    factor_G_gDW =factors[20]
    factor_Lx_1 = factors[21]
    factor_Lx_2 = factors[22]
    factor_Lx_3 = factors[23]
    factor_kdeg_m = factors[24]
    factor_kdeg_p = factors[25]
    factor_dil_rate = factors[26]
    #maybe add the k-binding coefficients, and the n

    parameters_original = [
    K_T;
    K_X;
    e_X;
    e_L;
    R_X_gDW;
    R_L_gDW;
    tau_m1;
    tau_p1;
    W_I1;
    W_12;
    W_13;
    W_23;
    W_11;
    k_I1;
    n_I1;
    k_12;
    n_12;
    k_13;
    n_13;
    G_gDW;
    Lx_1;
    Lx_2;
    Lx_3;
    kdeg_m;
    kdeg_p;
    dil_rate;
    ]

    #Calculate the h in all cases
    h = zeros(length(factors))
    h[1] = abs(K_T - factor_K_L*K_T)
    h[2] = abs(K_X - factor_K_X*K_X)
    h[3] = abs(e_X - factor_eX*e_X)
    h[4] = abs(e_L - factor_eL*e_L)
    h[5] = abs(R_X_gDW - factor_RX*R_X_gDW)
    h[6] = abs(R_L_gDW - factor_RL*R_L_gDW)
    h[7] = abs(tau_m1 - factor_tauX*tau_m1)
    h[8] = abs(tau_p1 - factor_tauL*tau_p1)
    h[9] = abs(W_I1 - factor_wI1*W_I1)
    h[10] = abs(W_12 - factor_w12*W_12)
    h[11] = abs(W_13 - factor_w13*W_13)
    h[12] = abs(W_23 - factor_w23*W_23)
    h[13] = abs(W_11 - factor_wjj*W_11)
    h[14] = abs(k_I1 - factor_k_I1*k_I1)
    h[15] = abs(n_I1 - factor_n_I1*n_I1)
    h[16] = abs(k_12 - factor_k_12*k_12)
    h[17] = abs(n_12 - factor_n_12*n_12)
    h[18] = abs(k_13 - factor_k_13*k_13)
    h[19] = abs(n_13 - factor_n_13*n_13)
    h[20] = abs(G_gDW - factor_G_gDW*G_gDW)
    h[21] = abs(Lx_1 - factor_Lx_1*Lx_1)
    h[22] = abs(Lx_2 - factor_Lx_2*Lx_2)
    h[23] = abs(Lx_3 - factor_Lx_3*Lx_3)
    h[24] = abs(kdeg_m - factor_kdeg_m*kdeg_m)
    h[25] = abs(kdeg_p - factor_kdeg_p*kdeg_p)
    h[26] = abs(dil_rate - factor_dil_rate*dil_rate)


    #---------------------------------------------------------------------------------
    #Change all the parameters by the input factors
    #--------------------------------------------------------------------------------

    #Binding COefficients
    K_T = K_T*factor_K_L
    K_X = K_X*factor_K_X

    #Elongation Rates
    ke_m1 = ke_m1*factor_eX;
    ke_m2 = ke_m2*factor_eX;
    ke_m3 = ke_m3*factor_eX;
    ke_p1 = ke_p1*factor_eL;
    ke_p2 = ke_p2*factor_eL;
    ke_p3 = ke_p3*factor_eL;

    #Numbers of Ribosomes/RNAP
    R_X_gDW = R_X_gDW*factor_RX
    R_L_gDW = R_L_gDW*factor_RL

    #Change the tau's
    tau_original_X = tau_m1
    tau_original_L = tau_p1
    tau_m1 = tau_m2 = tau_m2 =  tau_original_X*factor_tauX
    tau_p1 = tau_p2 = tau_p3 = tau_original_L*factor_tauL

    #Change the weights
    W_I1 = factor_wI1*W_I1
    W_12 = factor_w12*W_12
    W_13 = factor_w13*W_13
    W_23 = factor_w23*W_23
    W_11 = factor_wjj*W_11

    #Changing the binding k and the n
    k_I1 = factor_k_I1*k_I1;
    n_I1 = factor_n_I1*n_I1;
    k_12 = factor_k_12*k_12;
    n_12 = factor_n_12*n_12;
    k_13 = factor_k_13*k_13;
    n_13 = factor_n_13*n_13;
    G_gDW = factor_G_gDW*G_gDW;

    #Change lengths
    Lx_1 = factor_Lx_1*Lx_1
    Lx_2 = factor_Lx_2*Lx_2
    Lx_3 = factor_Lx_3*Lx_3

    #Change degradation rates
    kdeg_m = factor_kdeg_m*kdeg_m
    kdeg_p = factor_kdeg_p*kdeg_p
    dil_rate = factor_dil_rate*dil_rate

    #-----------------------------------------------------------------
    #Set-up Times for Discrete Computation
    #-----------------------------------------------------------------
    #Time ranges
    tstart = 0; #min
    tstep = 1; #min
    tend = 560; #min
    tSim = collect(tstart:tstep:tend)

    #Empty array for outputs
    x_total = zeros(length(tSim),6);

    #Set-up A
    A0 = zeros(6,6)
    A0[1,1]=A0[2,2]=A0[3,3] = -kdeg_m-dil_rate
    A0[4,4]=A0[5,5]=A0[6,6] = -kdeg_p-dil_rate

    #Set-up S
    S0 = zeros(6,6)
    S0[1,1]=S0[2,2]=S0[3,3]=S0[4,4]=S0[5,5]=S0[6,6]=1

    Identity_matrix = S0;

    #-----------------------------------------------------------------
    #Set-up Initial Values
    #-----------------------------------------------------------------

    #Set initial values for x and I
    x0 = [0.0; #m1
          0.0; #m2
          0.0; #m3
          0.0; #p1
          0.0; #p2
          0.0; #p3
          ]
    m1= x_total[1,1]
    m2 = x_total[1,2]
    m3 = x_total[1,3]
    p1 = x_total[1,4]
    p2 = x_total[1,5]
    p3 = x_total[1,6]
    I = 0.0;

    #Calculate initial f
        #Calculate the fraction f as a function of input concentrations
        f(p,k,n) = (p^n)/(k^n + p^n);
        f_I1 = f(I,k_I1,n_I1);
        f_12 = f(p1,k_12,n_12);
        f_13 = f(p1,k_13,n_13);
        f_23 = f(p2,k_23,n_23);
        #print(f_13)

    #Calculate initial rates
        #Calculate the u control functions for m1, m2, and m3
        u_m1 = (W_11 + W_I1*f_I1)/(1+W_11 + W_I1*f_I1);
        u_m2 = (W_22 + W_12*f_12)/(1+W_22 + W_12*f_12);
        u_m3 = (W_33 + W_13*f_13)/(1 + W_33 + W_23*f_23 + W_13*f_13);
        #for proteins assume translation happens at the kinetic limit
        u_p1 = 1;
        u_p2 = 1;
        u_p3 = 1;
        #Set-up rX and rL
        r_X1 = ke_m1*R_X_gDW*(G_gDW/(tau_m1*K_X + (tau_m1+1)*G_gDW));
        r_X2 = ke_m2*R_X_gDW*(G_gDW/(tau_m2*K_X + (tau_m2+1)*G_gDW));
        r_X3 = ke_m3*R_X_gDW*(G_gDW/(tau_m3*K_X + (tau_m3+1)*G_gDW));
        r_L1 = ke_p1*R_L_gDW*(m1/(tau_p1*K_T + (tau_p1+1)*m1));
        r_L2 = ke_p2*R_L_gDW*(m2/(tau_p2*K_T + (tau_p2+1)*m2));
        r_L3 = ke_p3*R_L_gDW*(m3/(tau_p3*K_T + (tau_p3+1)*m3));
        #Set-up TX and TL rates
        TX_1 = r_X1 * u_m1;
        TX_2 = r_X2 * u_m2;
        TX_3 = r_X3 * u_m3;
        TL_1 = r_L1 * u_p1;
        TL_2 = r_L2 * u_p2;
        TL_3 = r_L3 * u_p3;
        r0 = [TX_1; #m1
              TX_2; #m2
              TX_3; #m3
              TL_1; #p1
              TL_2; #p2
              TL_3; #p3
              ]

        #Set up array to store rates at eery given timestep
        r_total = zeros(length(tSim),6)
        r_total[1,1] = r0[1];
        r_total[1,2] = r0[2];
        r_total[1,3] = r0[3];
        r_total[1,4] = r0[4];
        r_total[1,5] = r0[5];
        r_total[1,6] = r0[6];


    #----------------------------------------------------------------------------------------------
    #For - Loop for Discrete Calculation
    #----------------------------------------------------------------------------------------------
    #Define your A_hat and S_hat
    A_hat = exp((A0*tstep))
    S_hat = inv(A0)*(A_hat-Identity_matrix)*S0

    #Designate when to start the Inducer
    I_start = 260; #min
    #How much inducer
    I_level = 10.0; #mM

    f_total = zeros(length(tSim),4);
    u_total = zeros(length(tSim),6);
    dxdt_total = zeros(length(tSim),6);

    for k in 1:(length(tSim)-1)
          #set time
          t = tSim[k];

          #-----------------------------------------------------------------
          #Set-up Inducer Concentration
          #-----------------------------------------------------------------
          if t < I_start #less than I start
                I = 0.0 #mM
          else #after I_start
                I = I_level #mM
          end

          #---------------------------------------------------------------------
          #Use current values (at k) to calculate the x for the next time step
          #---------------------------------------------------------------------

          #get current x and r
          x_k = zeros(6,1); r_k = zeros(6,1);
          x_k[1]=x_total[k,1];
          x_k[2]=x_total[k,2];
          x_k[3]=x_total[k,3];
          x_k[4]=x_total[k,4];
          x_k[5]=x_total[k,5];
          x_k[6]=x_total[k,6];

          r_k[1]=r_total[k,1];
          r_k[2]=r_total[k,2];
          r_k[3]=r_total[k,3];
          r_k[4]=r_total[k,4];
          r_k[5]=r_total[k,5];
          r_k[6]=r_total[k,6];

          #calculate the new x_kplus1
          x_kplus1 = A_hat*x_k + S_hat*r_k

          #update x_total
          x_total[(k+1),1]=x_kplus1[1]
          x_total[(k+1),2]=x_kplus1[2]
          x_total[(k+1),3]=x_kplus1[3]
          x_total[(k+1),4]=x_kplus1[4]
          x_total[(k+1),5]=x_kplus1[5]
          x_total[(k+1),6]=x_kplus1[6]

          #---------------------------------------------------------------------
          #Update the next steps
          #---------------------------------------------------------------------
          #Extract variables
          m1=x_kplus1[1]
          m2=x_kplus1[2]
          m3=x_kplus1[3]
          p1=x_kplus1[4]
          p2=x_kplus1[5]
          p3=x_kplus1[6]

          f_I1 = f(I,k_I1,n_I1);
          f_12 = f(p1,k_12,n_12);
          f_13 = f(p1,k_13,n_13);
          f_23 = f(p2,k_23,n_23);

          #Save f for troublshooting
          f_total[(k+1),1] = f_I1;
          f_total[(k+1),2] = f_12;
          f_total[(k+1),3] = f_13;
          f_total[(k+1),4] = f_23;

          #print(f_12)

          #Calculate the u control functions for m1, m2, and m3
          u_m1 = (W_11 + W_I1*f_I1)/(1+W_11 + W_I1*f_I1);
          u_m2 = (W_22 + W_12*f_12)/(1+W_22 + W_12*f_12);
          u_m3 = (W_33 + W_13*f_13)/(1 + W_33 + W_23*f_23 + W_13*f_13);
          #for proteins assume translation happens at the kinetic limit
          u_p1 = 1;
          u_p2 = 1;
          u_p3 = 1;

          u_total[(k+1),1] = u_m1
          u_total[(k+1),2] = u_m2
          u_total[(k+1),3] = u_m3
          u_total[(k+1),4] = u_p1
          u_total[(k+1),5] = u_p2
          u_total[(k+1),6] = u_p3

          #Update r_X and r_L
          r_X1 = ke_m1*R_X_gDW*(G_gDW/(tau_m1*K_X + (tau_m1+1)*G_gDW));
          r_X2 = ke_m2*R_X_gDW*(G_gDW/(tau_m2*K_X + (tau_m2+1)*G_gDW));
          r_X3 = ke_m3*R_X_gDW*(G_gDW/(tau_m3*K_X + (tau_m3+1)*G_gDW));
          r_L1 = ke_p1*R_L_gDW*(m1/(tau_p1*K_T + (tau_p1+1)*m1));
          r_L2 = ke_p2*R_L_gDW*(m2/(tau_p2*K_T + (tau_p2+1)*m2));
          r_L3 = ke_p3*R_L_gDW*(m3/(tau_p3*K_T + (tau_p3+1)*m3));
          #Update TX and TL
          TX_1 = r_X1 * u_m1;
          TX_2 = r_X2 * u_m2;
          TX_3 = r_X3 * u_m3;
          TL_1 = r_L1 * u_p1;
          TL_2 = r_L2 * u_p2;
          TL_3 = r_L3 * u_p3;
          r_kplus1 = [TX_1;TX_2;TX_3;TL_1;TL_2;TL_3]

          #put into r total
          r_total[(k+1),1]=TX_1
          r_total[(k+1),2]=TX_2
          r_total[(k+1),3]=TX_3
          r_total[(k+1),4]=TL_1
          r_total[(k+1),5]=TL_2
          r_total[(k+1),6]=TL_3

          dxdt_kplus1 = A0*x_kplus1 + S0*r_kplus1;
          dxdt_total[(k+1),1] = dxdt_kplus1[1]
          dxdt_total[(k+1),2] = dxdt_kplus1[2]
          dxdt_total[(k+1),3] = dxdt_kplus1[3]
          dxdt_total[(k+1),4] = dxdt_kplus1[4]
          dxdt_total[(k+1),5] = dxdt_kplus1[5]
          dxdt_total[(k+1),6] = dxdt_kplus1[6]

      end

    return (x_total,h,parameters_original);


end
