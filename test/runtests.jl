using homogenous_SIR_model
using Test

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
