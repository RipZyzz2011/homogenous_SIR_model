module homogenous_SIR_model

using Plots
using DifferentialEquations
#using Test

export define_town_model, solve_system, plot_model_solution, errors

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

function errors(model, data)
    err_sum = 0
    for i in 1:(length(data))
        err_sum += (model[i] - data[i])^2
    end

    return err_sum
end

end

