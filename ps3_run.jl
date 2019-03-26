#THP1 Problem 3
#Sneha Kabaria

#Include flux function
include("flux.jl")
#Plotting Package
using Pkg
Pkg.add("PyPlot")
using PyPlot

#-----------------------------------------------------------------------------
#Given Information
#----------------------------------------------------------------------------

#Given Information
L_p= 308; #AA, length protein
L_g = 924; #nt, length gene
gene_conc = 5e-3; #uM
rxn_volume = 15; #uL
R_X = 0.15; #uM
R_L = 1.6; #uM
e_X = 60; #nt/sec
e_L = 16.5; #nt/sec
K_X = 0.3; #uM
K_L = 57.0; #uM
tau_X = 2.7;
tau_L = 0.8;
kdeg_m = 8.35/3600; #s^-1
kdeg_p = 9.9e-3/3600; #s^-1
L_g_char = 1000; #nt
L_p_char = 330; #aa

#-------------------------------------------------------------------------------
#Stoichiometric Matrix
#-------------------------------------------------------------------------------

n = L_g_char; #nt
a = L_p_char; #aa

#Set-up Stoichiometric Matrix
S = zeros(17,15);
#v1-v6, b1-b9
S[1,1]=-1; S[1,2] = 1; #G
S[2,1]=-1; S[2,2] = 1; #RNAP
S[3,1] = 1; S[3,2] = -1; #G*
S[4,2] = -n; S[4,8] = 1; #NTP
S[5,2] = 1; S[5,3] = -1; S[5,4]=-1; S[5,5]=1; #mRNA
S[6,2]=2*n; S[6,5]=2*a; S[6,6]=2; S[6,15]=-1; #Pi
S[7,3]=n; S[7,10] = -1; #NMP
S[8,4]=-1; S[8,5] = 1; #rib
S[9,4]=1; S[9,5]=-1; #rib*
S[10,5] = -a; S[10,6]=1; #AAtRNA
S[11,5]=-2*a; S[11,13]=1; #GTP
S[12,5] = 2*a; S[12,14] = -1; #GDP
S[13,6] = -1; S[13,7] = 1; #AA
S[14,5] = a; S[14,6] = -1;#tRNA
S[15,6] =-1; S[15,11] = 1; #ATP
S[16,6] = 1; S[16,12] = -1;#AMP
S[17,5] = 1; S[17,9] = -1; #protein

#Print "S" in the command line to see the stoichiometric matrix

#-------------------------------------------------------------------------------
#Set up control function parameters and #Calculate r_X and r_L
#-------------------------------------------------------------------------------

#actual elongation rates for mRNA and protein
v_X_char = e_X/L_g_char;
v_X = v_X_char*L_g_char/L_g;
v_L_char = e_L/L_p_char;
v_L = v_L_char*L_p_char/L_p;

#Set-up Rates
#Transcription rate r_X
r_X = v_X*R_X * ( gene_conc / ( tau_X*K_X + (tau_X+1)*gene_conc ) );

#Find the control function
w1 = 0.26; #mM
w2 = 300.0; #mM
K = 0.30; #mM
n_binding = 1.5;
#Input an I
fraction(I) = (I^n_binding) / (K^n_binding+I^n_binding);
#Define u as a function of I
u(I) = (w1 + fraction(I)*w2)/(1+w1 + fraction(I)*w2);
TX(I) = r_X*u(I);
#Find steady-state m1 (mRNA) concentration
mRNA(I) = TX(I)/kdeg_m;
mRNA_max = r_X/kdeg_m
#Translation rate r_L
r_L(I) = v_L*R_L * ( mRNA(I) / (tau_L*K_L + (tau_L+1)*mRNA(I) ) );
r_L_max = v_L*R_L * ( mRNA_max / (tau_L*K_L + (tau_L+1)*mRNA_max ) );

#---------------------------------------------------------------------------
#Set I for testing
#---------------------------------------------------------------------------

#I = 10;

#-----------------------------------------------------------------------------
#Set Inputs for flux function and iterate through different I
#-----------------------------------------------------------------------------
I_lower = 0.0001; #mM
step1 = .0001; #mM
I_mid1 = 0.001; #mM
I_Sim1 = collect(I_lower:step1:I_mid1)
step2 = .001; #mM
I_mid2 = .01; #mM
I_Sim2 = collect(I_mid1:step2:I_mid2)
step3 = .01; #mM
I_mid3 = 0.1; #mM
I_Sim3 = collect(I_mid2:step3:I_mid3)
step4 = 0.1 #mM
I_upper = 10.0; #mM
I_Sim4 = collect(I_mid3:step4:I_upper)
I_Sim = [I_Sim1;I_Sim2;I_Sim3;I_Sim4]

#Gather v5 rates (make empty array)
v5_total = zeros(length(I_Sim));
#Gather protein at SS (make empty array); protein = v5/kdeg_p
protein_total = zeros(length(I_Sim));


for i in 1:(length(I_Sim))
    I = I_Sim[i]
    #Define Bounds
    b_lowerbound = -100000.0; #uM/hr
    b_upperbound = 100000.0; #uM/hr
    v_lowerbound = 0.0; #uM/hr
    v_upperbound = Inf;
    v5_lowerbound = 0.0;
    v5_upperbound = r_L(I);
    v2_lowerbound = TX(I);
    v2_upperbound = TX(I);
    #Make the bounds vector
    default_bounds_array = zeros(15,2);
    #Fill the bounds array
    for k in 1:(length(default_bounds_array[:,1]))
        if k == 2
            default_bounds_array[k,1] = v2_lowerbound
            default_bounds_array[k,2] = v2_upperbound
        elseif k == 5
            default_bounds_array[k,1] = v5_lowerbound
            default_bounds_array[k,2] = v5_upperbound
        elseif k > 7
            default_bounds_array[k,1] = b_lowerbound
            default_bounds_array[k,2] = b_upperbound
        else
            default_bounds_array[k,1] = v_lowerbound
            default_bounds_array[k,2] = v_upperbound
        end
    end

    objective_coefficient_array= zeros(15);
    objective_coefficient_array[5] = -1.0;

    species_bounds_array = zeros(17,2);

    flux_answer = calculate_optimal_flux_distribution(S,default_bounds_array,species_bounds_array,objective_coefficient_array)
    #Pull out the answers
    objective_value= flux_answer[1] #uM per second
    calculated_flux_array= flux_answer[2] #uM per second
    dual_value_array = flux_answer[3]
    uptake_array2 = flux_answer[4]
    print(I)
    if I == 10.0
        calculated_flux_array10 = calculated_flux_array
    else
        #nothing
    end

    #Calculate SS protein
    v5 = calculated_flux_array[5];
    protein_current = v5/kdeg_p;
    #Record SS protein
    v5_total[i] = v5;
    protein_total[i] = protein_current #uM

end

#----------------------------------------------------------------------
#Plotting
#----------------------------------------------------------------------
figure(figsize=(4*4,3*4))
figure(1)
semilogx(I_Sim,protein_total,color="black")
xlabel("I (mM))")
ylabel("protein (uM)")

figure(figsize=(4*4,3*4))
figure(2)
semilogx(I_Sim,v5_total,color="black")
xlabel("I (mM))")
ylabel("v5 (uM/s)")

#------------------------------------------------------------------------------
#part C -- get the values for when I = 0.01 so that you can decide which is the most important metabolite
#------------------------------------------------------------------------------

#If I = 0.01 what is the flux opitmization

I =0.01
#Define Bounds
b_lowerbound = -100000.0; #uM/hr
b_upperbound = 100000.0; #uM/hr
v_lowerbound = 0.0; #uM/hr
v_upperbound = Inf;
v5_lowerbound = 0.0;
v5_upperbound = r_L(I);
v2_lowerbound = TX(I);
v2_upperbound = TX(I);
#Make the bounds vector
default_bounds_array = zeros(15,2);
#Fill the bounds array
for k in 1:(length(default_bounds_array[:,1]))
    if k == 2
        default_bounds_array[k,1] = v2_lowerbound
        default_bounds_array[k,2] = v2_upperbound
    elseif k == 5
        default_bounds_array[k,1] = v5_lowerbound
        default_bounds_array[k,2] = v5_upperbound
    elseif k > 7
        default_bounds_array[k,1] = b_lowerbound
        default_bounds_array[k,2] = b_upperbound
    else
        default_bounds_array[k,1] = v_lowerbound
        default_bounds_array[k,2] = v_upperbound
    end
end

objective_coefficient_array= zeros(15);
objective_coefficient_array[5] = -1.0;

species_bounds_array = zeros(17,2);

flux_answer = calculate_optimal_flux_distribution(S,default_bounds_array,species_bounds_array,objective_coefficient_array)
#Pull out the answers
objective_value_pt1= flux_answer[1] #uM per second
calculated_flux_array_pt1= flux_answer[2] #uM per second
dual_value_array_pt1 = flux_answer[3]
uptake_array2_pt1 = flux_answer[4]

#------------------------------------------------------------------------------
#part C -- get the values for when I = 0.01 so that you can decide which is the most important metabolite
#------------------------------------------------------------------------------

#If I = 0.01 what is the flux opitmization

I =5
#Define Bounds
b_lowerbound = -100000.0; #uM/hr
b_upperbound = 100000.0; #uM/hr
v_lowerbound = 0.0; #uM/hr
v_upperbound = Inf;
v5_lowerbound = 0.0;
v5_upperbound = r_L(I);
v2_lowerbound = TX(I);
v2_upperbound = TX(I);
#Make the bounds vector
default_bounds_array = zeros(15,2);
#Fill the bounds array
for k in 1:(length(default_bounds_array[:,1]))
    if k == 2
        default_bounds_array[k,1] = v2_lowerbound
        default_bounds_array[k,2] = v2_upperbound
    elseif k == 5
        default_bounds_array[k,1] = v5_lowerbound
        default_bounds_array[k,2] = v5_upperbound
    elseif k > 7
        default_bounds_array[k,1] = b_lowerbound
        default_bounds_array[k,2] = b_upperbound
    else
        default_bounds_array[k,1] = v_lowerbound
        default_bounds_array[k,2] = v_upperbound
    end
end

objective_coefficient_array= zeros(15);
objective_coefficient_array[5] = -1.0;

species_bounds_array = zeros(17,2);

flux_answer = calculate_optimal_flux_distribution(S,default_bounds_array,species_bounds_array,objective_coefficient_array)
#Pull out the answers
objective_value_5= flux_answer[1] #uM per second
calculated_flux_array_5= flux_answer[2] #uM per second
dual_value_array_5 = flux_answer[3]
uptake_array2_5 = flux_answer[4]
