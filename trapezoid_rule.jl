#Define the trapezoid rule

function trapezoid_rule(values_vector,time_step)

    number_of_values = length(values_vector)

    integral_total = 0; #Set this values for the summation

    for i in 1:(length(values_vector)-1)
        #Absolute values as part of the N-matrix formulas 
        b1 = abs(values_vector[i])
        b2 = abs(values_vector[i+1])
        height = time_step

        area = (b1+b2)/2*height

        integral_total = integral_total+area
    end

    return (integral_total)
end
