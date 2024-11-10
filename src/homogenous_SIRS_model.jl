module homogenous_SIRS_model

using DifferentialEquations
using Plots
using Measurements
using homogenous_SIR_model

export town_SIRS!, town_SIRS_Intervention!, error_squares
#SIR model that now incorporates re-infection and severe illness state
function town_SIRS!(dpop, pop, param, t)
    N = sum(pop)
    # alpha represents the re-infection rate
    # gamma_s represents the recovery rate of severe infection
    # p_s represents the Probability of developing severe infection
    c, Beta_c, gamma, alpha, p_s, gamma_s = param
    S, I, Is, R = pop
    lambda = c * Beta_c * I / N
    R_0 = c * 1/gamma * Beta_c # Reproduction number
    
    dpop[1] = -lambda * S + alpha * R# dS = -lambda*S
    dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
    #Severe infection
    dpop[3] = gamma * p_s * I - gamma_s * Is 
    dpop[4] = (1 - p_s) * gamma * I + gamma_s * Is - alpha * R # dR = gamma * R

end

#SIRS model that integrates an intervention
#= 
  It is anticipated that this intervention will reduce the probability of transmission from an infected person
   to a susceptible person by about 30% (this is the efficacy of the intervention). 
   The Department of Health expects approximately 
   80% of the population will comply with their directive to implement the public health intervention.
=#
function town_SIRS_Intervention!(dpop, pop, param, t)
    N = sum(pop)
    c, Beta_c, gamma, alpha, p_s, gamma_s, epsilon, phi = param
    #epsilon represents efficacy of intervention
    #phi represents proportion of intervention take-up
    S, I, Is, R = pop
    #Lambda: Force of Infection
    lambda = c * (1 - epsilon * phi) * Beta_c * I / N
    R_0 = c * 1/gamma * Beta_c # Reproduction number
    
    dpop[1] = -lambda * S + alpha * R# dS = -lambda*S
    dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
    #Severe infection
    dpop[3] = gamma * p_s * I - gamma_s * Is 
    dpop[4] = (1 - p_s) * gamma * I + gamma_s * Is - alpha * R # dR = gamma * R

end

#Calculate the sum of the error at each datapoint between the real values and the model values
#Squares error between values, useful for evaluating the efficacy of the beta parameter
function error_squares(model, data)
    err_sum = 0
    for i in 1:(length(data))
        err_sum += (model[i] - data[i])^2
    end

    return err_sum
end

end