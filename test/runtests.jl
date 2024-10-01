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
