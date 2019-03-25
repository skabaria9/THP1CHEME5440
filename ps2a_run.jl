#Problem 2 - THP1 - Sneha Kabaria

#Included packages
using Pkg
Pkg.add("LinearAlgebra")
Pkg.add("PyPlot")
using LinearAlgebra
using PyPlot
#Included sub-parameters and functions
include("ps2_parameters.jl")

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

#----------------------------------------------------------------------
#Plotting
#----------------------------------------------------------------------

  #print(x_total)

  m1 = x_total[:,1] #nmol/gDW
  m2 = x_total[:,2] #nmol/gDW
  m3 = x_total[:,3] #nmol/gDW
  p1 = x_total[:,4] #nmol/gDW
  p2 = x_total[:,5] #nmol/gDW
  p3 = x_total[:,6] #nmol/gDW

  #Make I vector for graphing
  I_mM_graph = zeros(length(tSim)) #set default all zeros in I
  for i in 1:length(tSim)
        if tSim[i] < I_start
              #print(tSim[i])
              I_mM_graph[i] = 0.0
        else
              I_mM_graph[i] = I_level
        end
  end

  #in mM Plotting
  using PyPlot
  figure(figsize=(4*4,3*4))
  figure(1)
  plot(tSim,m1,color="black")
  plot(tSim,m2,color="blue")
  plot(tSim,m3,color="red")
  plot(tSim,I_mM_graph,color="green")
  xlabel("time (min)")
  ylabel("Concentration mRNA (nmol/gDW)")
  legend(["m1","m2","m3"])

  figure(figsize=(4*4,3*4))
  figure(2)
  plot(tSim,p1,color="black")
  plot(tSim,p2,color="blue")
  plot(tSim,p3,color="red")
  plot(tSim,I_mM_graph*100,color="green")
  xlabel("time (min)")
  ylabel("Concentration protein (nmol/gDW)")
  legend(["p1","p2","p3"])
