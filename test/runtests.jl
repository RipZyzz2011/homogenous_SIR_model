using homogenous_SIR_model
using Test
using Pkg
Pkg.add("Revise")
Pkg.add("DifferentialEquations")
using Revise
using DifferentialEquations


@testset "townSIR_basic_set" begin
    #Test that no one gets infected when a population has I = 0
    S = 100
    I = 0
    R = 0
    N_initial = S + I + R
    lambda = 0.1
    gamma = 0.1

    t_span = (0, 50)

    model = define_town_model(:basic, [lambda, gamma], [S, I, R], t_span)
    solve_system(model)
    
    @test I == 0

    #Test that the population remains the same size that it begun at
    S = 90
    I = 10
    R = 0
    N_initial = S + I + R
    lambda = 0.1
    gamma = 0.1

    t_span = (0, 50)

    model = define_town_model(:basic, [lambda, gamma], [S, I, R], t_span)
    solve_system(model)
    
    @test (S + I + R) == N_initial
end

@testset "townSIR_herd_set" begin
    # Test that the town stops getting infected once herd immunity threshold is met
    S = 90
    I = 10
    R = 0
    N_initial = S + I + R
    c = 10
    Beta_c = 0.10
    gamma = 0.05
    R_0 = Beta_c * c / gamma
    p_c = 1 - 1/R_0
    param = [c, Beta_c, gamma] #c, Beta_c, gamma
 

    t_span = (0, 50)

    model = define_town_model(:herd, param, [S, I, R], t_span)
    sol = solve_system(model)
    #Test that the rate of infection goes to 0
    @test sol[end][2] < sol[1][2]  # Infections should decrease when herd immunity threshold is met
    
end

# Test the implementation of the severe illness and resusceptance
@testset "townSIRS_set" begin
    # Make use of the data parameters from the project task
    #Town Population, assume that there was a single person infected on day 1
    N = 6000
    I = 1
    S = N - I
    I_s = 0
    R = 0

    c = 8 #Number of daily contacts on average
    gamma = 1/7 # Daily rate of recovery if it takes 7 days to recover typically
    gamma_s = 1/14 # Daily rate of recovery from severe illness
    p_s = 0.20 # Average probability of severe infection
    alpha = 1/30 # Daily rate of resusceptance if the average time for it is a month

    Beta = 0.03697
    param = [c, Beta, gamma, alpha, p_s, gamma_s]
    t_span = (0, 150)
    pop0 = [S, I, I_s, R]

    # Edge case testing
    # See that the number of infections goes to 0 if alpha = 0
    alpha = 0
    model =  ODEProblem(town_SIRS!, pop0, t_span, param)
    sol = solve(model, saveat = 1)
    I_model = [u[2] for u in sol.u]
    Is_model = [u[3] for u in sol.u]

    @test I_model[150] == 0 && Is_model[150] == 0
     

    # Likewise, test that if alpha is large enough, the number of infected never goes
    # to 0
    alpha = 1/30
    model =  ODEProblem(town_SIRS!, pop0, t_span, param)
    sol = solve(model, saveat = 1)
    I_model = [u[2] for u in sol.u]
    Is_model = [u[3] for u in sol.u]

    @test I_model[150] != 0 && Is_model[150] != 0
     
end

# test the intervention ODE functions as intended
@testset "townSIRSIntervention_set" begin

    # See that the steady-state infection values are less when the intervention
    # is implemented
    N = 6000
    I = 1
    S = N - I
    I_s = 0
    R = 0


    # Constructing parameters from the data
    c = 8 #Number of daily contacts on average
    gamma = 1 / 7 # Daily rate of recovery if it takes 7 days to recover typically
    gamma_s = 1 / 14
    p_s = measurement(0.20) # Average probability of severe infection
    alpha = 1 / 30 # Daily rate of resusceptance if the average time for it is a month
    epsilon = 0.30 # Efficacy of intervention
    phi = 0.65 # Proportion of population that will adhere to the intervention

    #Using the value of beta that best matches the current data as of 21/10/2024
    Beta = measurement(0.035)
    param_no_int = [c, Beta, gamma, alpha, p_s, gamma_s]
    param_int = [c, Beta, gamma, alpha, p_s, gamma_s, epsilon, phi]
    # Implement the intervention at 30 days, simulate for another 30
    t_span_half_1 = (0, 30)
    t_span_half_2 = (30, 365)
    pop0_no_int = [S, I, I_s, R]
    pop0_int = [S, I, I_s, R]

    # Simulate two sets of models, with one set implementing the intervention after day 30
    # The model needs to be broken into two time sets in order to run properly

    model_no_int1 = ODEProblem(town_SIRS!, pop0_no_int, t_span_half_1, param_no_int)
    sol_no_int1 = solve(model_no_int1, saveat=1)
    model_no_int2 = ODEProblem(town_SIRS!, sol_no_int1.u[31], t_span_half_2, param_no_int)
    sol_no_int2 = solve(model_no_int2, saveat=1)

    model_int1 = ODEProblem(town_SIRS!, pop0_int, t_span_half_1, param_no_int)
    sol_int1 = solve(model_int1, saveat=1)
    model_int2 = ODEProblem(town_SIRS_Intervention!, sol_int1.u[31], t_span_half_2, param_int)
    sol_int2 = solve(model_int2, saveat=1)

    # The data of interest is the number of infected, obtain from solution as so
    I_model_no_int1 = [u[2].val for u in sol_no_int1.u]
    I_model_no_int2 = [u[2].val for u in sol_no_int2.u]
    append!(I_model_no_int1, I_model_no_int2)

    I_model_int1 = [u[2].val for u in sol_int1.u]
    I_model_int2 = [u[2].val for u in sol_int2.u]
    append!(I_model_int1, I_model_int2)
    
    @test I_model_int1[300] < I_model_no_int1[300]
end
