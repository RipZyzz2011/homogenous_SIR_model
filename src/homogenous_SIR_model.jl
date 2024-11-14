module homogenous_SIR_model

using Plots
using DifferentialEquations
#using Test

export define_town_model, solve_system, plot_model_solution, errors, town_SIRS!, town_SIRS_Intervention!

# Defines an SIR model for the town based on the user input of model type given as a symbol.
# The input vectors define the parameters for operation (rate of infection, contact, etc), the initial 
# popultations of each category, and the time span at which the model is run for in days
function define_town_model(model::Symbol, parameters::Vector, initial_pop::Vector, t_span::Tuple)
    #=Function: town_SIR
    # Creates a Susceptible, Infected, Recovered (SIR) mathematical model of a town
    on the assumption and basis of a homogenous population, and that recovered people
    aren't susceptible to infection again

    dpop: vector containing the rate of change of each population category
    pop: vector containing the current number of each population category
    param: vector containing the parameters of the differential DifferentialEquations
    t: Beginning and starting time that the model operates from in days

    Each version of the function operates differently based on the size of the param vector, 
    =#

    # Operates with a specified lambda
    function town_SIR_basic!(dpop, pop, param::Vector, t)
        N = sum(pop)
        lambda, gamma = param
    
        S, I, R = pop
    
        dpop[1] = -lambda * S # dS = -lambda*S
        dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
        dpop[3] = gamma * I # dR = gamma * R
    
    end


    # Operates with the force of infection equation for lambda instead of a provided rate
    function town_SIR_foi!(dpop, pop, param::Vector, t)
        N = sum(pop)
        c, Beta_c, gamma = param 
        #c: Number of daily contacts
        #Beta_c: Chance of contracting disease from infected contact
        #gamma: Probability of recovering from infection each day
    
    
        S, I, R = pop
        lambda = c * Beta_c * I / N #I(t) / N: Chance of infectious contact given homogenous population 
        #R_0 = c * 1/gamma * Beta_c # Reproduction number
    
        dpop[1] = -lambda * S # dS = -lambda*S
        dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
        dpop[3] = gamma * I # dR = gamma * R
    
    end

    #Operates using R0 and the herd immunity threshold criteria
    function town_SIR_herd!(dpop, pop, param::Vector, t)
        N = sum(pop)
        c, Beta_c, gamma = param 
        #c: Number of daily contacts
        #Beta_c: Chance of contracting disease from infected contact
        #gamma: Probability of recovering from infection each day
        S, I, R = pop
        
        
        R_0 = c * 1/gamma * Beta_c # Reproduction number
        p_c = 1 - 1/R_0

        #If the herd immunity threshold is reached, set Beta_c to 0)
        if((I/N) >= p_c)
            Beta_c = 0
        end
        
        lambda = c * Beta_c * I / N #I(t) / N: Chance of infectious contact given homogenous population 
        
        dpop[1] = -lambda * S # dS = -lambda*S
        dpop[2] = gamma * I * (c* 1/gamma * Beta_c * S/N - 1)
        #dpop[2] = gamma * I * (c* 1/gamma * Beta_c * S/N - 1)
        dpop[3] = gamma * I # dR = gamma * R
    
    
    end
    
    # Multiple dispatch selection to choose which model of SIR to pick
    if(model == :basic)    
        return ODEProblem(town_SIR_basic!, initial_pop, t_span, parameters)
    end
    if(model == :foi)    
        return ODEProblem(town_SIR_foi!, initial_pop, t_span, parameters)
    end
    if(model == :herd)    
        return ODEProblem(town_SIR_herd!, initial_pop, t_span, parameters)
    else 
        error("Unknown Model: $model")
    end
    
end


#= 
Function: solve_system
Produces the solution to the model created by define_town_model

=#
function solve_system(prob::ODEProblem)
    return solve(prob, saveat = 1)
end
#=
Function: plot_model_solution
Plots the solved model over the time span it was solved from
=#

function plot_model_solution(sol::ODESolution)
    plot(sol, xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "SIR Model")
end

################################################################################################
################################################################################################
#= The following set of functions operate without the multiple dispatch and plotting functionality
of the previous set to allow them to be implemented independently

Refer to the report for the explanation and breakdown of how the equations for the ODE system
are defined and the relevant parameter meanings
=#

#=
    Function: town_SIRS!
    ODE system that implements the established SIR model, but now incorporates a fourth, severe illness state
    as well as the capability for reinfection to occur
=# 
function town_SIRS!(dpop, pop, param, t)
    N = sum(pop)
    # alpha represents the re-infection rate
    # gamma_s represents the recovery rate of severe infection
    # p_s represents the Probability of developing severe infection
    c, Beta_c, gamma, alpha, p_s, gamma_s = param
    S, I, Is, R = pop # Current population stats
    lambda = c * Beta_c * I / N # Force of infection
    R_0 = c * 1/gamma * Beta_c # Reproduction number
    
    dpop[1] = -lambda * S + alpha * R# dS = -lambda*S
    dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
    #Severe infection
    dpop[3] = gamma * p_s * I - gamma_s * Is 
    dpop[4] = (1 - p_s) * gamma * I + gamma_s * Is - alpha * R # dR = gamma * R

end

#=
    Function: town_SIRS_Intervention!
    The same operating principle as town_SIRS however now incorporates an intervention
    defined by epsilon and phi that reduces the rate of a susceptible person becoming 
    infected
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

#Calculate the sum of the linear-least square error at each datapoint between the real values and the model values
#Squares error between values, useful for evaluating the efficacy of the beta parameter
function error_squares(model, data)
    err_sum = 0
    for i in 1:(length(data))
        err_sum += (model[i] - data[i])^2
    end

    return err_sum
end

end

