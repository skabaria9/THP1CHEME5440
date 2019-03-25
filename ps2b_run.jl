#Run 2b

#Import Statistics Tools -> Statistics.mean(A)
import Statistics

#Import Function for Iteration
include("2b_parameter_change_function.jl")
include("trapezoid_rule.jl")


#Define Factors of Change
#factors = [factor_K_L,factor_K_X,factor_eX,factor_eL,factor_RX,factor_RL,factor_tauX,factor_tauL,factor_wI1,factor_w12,factor_w13,factor_w23,factor_wjj]
percent_change = 0.05;
factor_increase = 1+percent_change;
factor_decrease = 1-percent_change;

#Percent Change that want to test
factors = ones(26);

#original_x
x_total, h_not, parameters_original = parameter_change(factors)

# #----------------------------------------------------------------------
# #Plotting
# #----------------------------------------------------------------------
#
#   #print(x_total)
#
#   m1 = x_total[:,1] #nmol/gDW
#   m2 = x_total[:,2] #nmol/gDW
#   m3 = x_total[:,3] #nmol/gDW
#   p1 = x_total[:,4] #nmol/gDW
#   p2 = x_total[:,5] #nmol/gDW
#   p3 = x_total[:,6] #nmol/gDW
#
#   #Make I vector for graphing
#   I_mM_graph = zeros(length(tSim)) #set default all zeros in I
#   for i in 1:length(tSim)
#         if tSim[i] < I_start
#               #print(tSim[i])
#               I_mM_graph[i] = 0.0
#         else
#               I_mM_graph[i] = I_level
#         end
#   end
#
#   #in mM Plotting
#   using PyPlot
#   figure(figsize=(4*4,3*4))
#   figure(1)
#   plot(tSim,m1,color="black")
#   plot(tSim,m2,color="blue")
#   plot(tSim,m3,color="red")
#   plot(tSim,I_mM_graph,color="green")
#   xlabel("time (min)")
#   ylabel("Concentration mRNA (nmol/gDW)")
#
#   figure(figsize=(4*4,3*4))
#   figure(2)
#   plot(tSim,p1,color="black")
#   plot(tSim,p2,color="blue")
#   plot(tSim,p3,color="red")
#   plot(tSim,I_mM_graph*100,color="green")
#   xlabel("time (min)")
#   ylabel("Concentration protein (nmol/gDW)")

sens_phase1_total =1 #dummy variable

#-------------------------------------------------------------------------------
#Sensitivity Analysis for All the Variables
#-------------------------------------------------------------------------------

#Sensitivity Array Description
#Rows = the x variables (m1,m2,m3,p1,p2,p3)
#Columns = factors changed
# parameters_original = [
# K_T;
# K_X;
# e_X;
# e_L;
# R_X_gDW;
# R_L_gDW;
# tau_m1;
# tau_p1;
# W_I1;
# W_12;
# W_13;
# W_23;
# W_11;
# k_I1;
# n_I1;
# k_12;
# n_12;
# k_13;
# n_13;
# G_gDW;
# ]

#Set-up Sensitivity Arrays
sensitivity_phase1 = zeros(6,26);
sensitivity_phase2early = zeros(6,26);
sensitivity_phase2late = zeros(6,26);

#Definte the different Phase Times that were time-averaged:
# The results were time-averaged:
#20 min windows were chosen for each of the following within the specified ranges
#--Phase I: 200-260 minutes ;
# --Phase II early: 261-350 minutes ;
# --Phase II late: 351-560 minutes ;
#Phase 1
phase1_start = 230; #min
phase1_end = 250; #min
sens_phase1_total = zeros(length(phase1_start:phase1_end),6) #Empty Sensitivty Matrix for Phase 1
#Phase2-early
phase2early_start = 290; #min
phase2early_end = 310; #min
sens_phase2early_total = zeros(length(phase2early_start:phase2early_end),6) #Empty Sensitivty Matrix for Phase 2early
#Phase2-late
phase2late_start = 470; #min
phase2late_end = 490; #min
sens_phase2late_total = zeros(length(phase2late_start:phase2late_end),6) #Empty Sensitivty Matrix for Phase 2early

species_bounds_array = 1; #holding variable

h_record = zeros(1,26); #somewhere to record all of the h's to report at the end

for i in 1:(length(factors))

    #Change the factor that is effected
    factors_plus_h = ones(26)
    factors_plus_h[i] = factor_increase;
    factors_minus_h = ones(26)
    factors_minus_h[i] = factor_decrease;

    #Get the f(p+h) and f(p-h) for all the parameters x
    x_plus_h, plus_h, no_use = parameter_change(factors_plus_h)
    x_minus_h, minus_h, no_use= parameter_change(factors_minus_h)


    #Subtract the two matrices to get the numerator for central difference
    x_difference = x_plus_h - x_minus_h;
    #print(x_difference[1:20,:])

    #print(plus_h[i])

    #Divide x_difference by 2*h to get the partial derivative
    partial_derivative = x_difference/(2*plus_h[i])

    #Record the h to report later
    h_record[i] = plus_h[i];

    #Determine p/x for given parameter
    p_current = parameters_original[i]
    #Calculate 1/x for every term in original matrix
    invert_matrix = x_total
    for term in 1:(length(invert_matrix))
        invert_matrix[term] = 1/invert_matrix[term]
    end
    #Multiply by the current parameter
    p_over_x = p_current*invert_matrix
    #Get dimensions
    row,col = size(p_over_x);

    #Phase 1
    phase1_partials = partial_derivative[phase1_start:phase1_end,:]
    phase1_p_over_x = p_over_x[phase1_start:phase1_end,:]
    for timepoint in 1:length(phase1_start:phase1_end) #Designate timepoint (row)
        for x_var in 1:6 #designate the xvar (column) that will do one-by-one
            #the sen_phase#_total matrix contains all of the s_ij (not absolute valued)
            sens_phase1_total[timepoint,x_var] = phase1_partials[timepoint,x_var] #*p_over_x[timepoint,x_var]
        end
    end
    #Fill in the phase1 sensitivity matrix
    #Trapezoid rule applied the integral over the range with the timestep 1, it also absolute values all of the inputs as part of the N-matrix formula
    for z in 1:6
        sensitivity_phase1[z,i] = 1/(phase1_end-phase1_start)*trapezoid_rule(sens_phase1_total[:,z],1)
    end

    #Phase 2 - early
    phase2early_partials = partial_derivative[phase2early_start:phase2early_end,:]
    phase2early_p_over_x = p_over_x[phase2early_start:phase2early_end,:]
    for timepoint in 1:length(phase2early_start:phase2early_end) #Designate timepoint (row)
        for x_var in 1:6 #designate the xvar (column) that will do one-by-one
            #the sen_phase#_total matrix contains all of the s_ij (not absolute valued)
            sens_phase2early_total[timepoint,x_var] = phase2early_partials[timepoint,x_var]*phase2early_p_over_x[timepoint,x_var]
        end
    end
    #Fill in the phase1 sensitivity matrix
    #Trapezoid rule applied the integral over the range with the timestep 1, it also absolute values all of the inputs as part of the N-matrix formula
    for z in 1:6
        sensitivity_phase2early[z,i] = 1/(phase2early_end-phase2early_start)*trapezoid_rule(sens_phase2early_total[:,z],1)
    end


    #Phase 2 - late
    phase2late_partials = partial_derivative[phase2late_start:phase2late_end,:]
    phase2late_p_over_x = p_over_x[phase2late_start:phase2late_end,:]
    for timepoint in 1:length(phase2late_start:phase2late_end) #Designate timepoint (row)
        for x_var in 1:6 #designate the xvar (column) that will do one-by-one
            #the sen_phase#_total matrix contains all of the s_ij (not absolute valued)
            sens_phase2late_total[timepoint,x_var] = phase2late_partials[timepoint,x_var]*phase2late_p_over_x[timepoint,x_var]
        end
    end
    #Fill in the phase1 sensitivity matrix
    #Trapezoid rule applied the integral over the range with the timestep 1, it also absolute values all of the inputs as part of the N-matrix formula
    for z in 1:6
        sensitivity_phase2late[z,i] = 1/(phase2late_end-phase2late_start)*trapezoid_rule(sens_phase2late_total[:,z],1)
    end

end

#-------------------------------------------------------------------------------
#Write Answers to Files
#-------------------------------------------------------------------------------
import Pkg
using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
using CSV, DataFrames
CSV.write("sensitivity_phase1_file.csv",  DataFrame(sensitivity_phase1), writeheader=false)
CSV.write("sensitivity_phase2early_file.csv",  DataFrame(sensitivity_phase2early), writeheader=false)
CSV.write("sensitivity_phase2late_file2.csv",  DataFrame(sensitivity_phase2late), writeheader=false)
CSV.write("sensitivity_h_list.csv",  DataFrame(h_record), writeheader=false)
