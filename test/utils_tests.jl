@testset "Testing the estimation of firing probabilities" begin
    config="""
    neurons:
      - id: n
        θ_syn: 5
        branches:
          - id: seg
            θ_syn: 1
    """

    (net,objects) = load_network(YAML_source=config)

    #Dict(id => (probabilities, weights))
    volleys = Dict(:n => (fill(0.5, 5), ones(5)), :seg => (fill(0.5, 5), ones(5)))

    @test P_superthreshold(3, [0.25,0.5,0.75],[1,2,3]) ≈ 0.25*0.5*0.75 + 0.25*0.5*(1.0-0.75) + 0.25*(1.0-0.5)*0.75 + (1.0-0.25)*0.5*0.75 + (1.0-0.25)*(1.0-0.5)*0.75
    @test P_fire(objects[:n], volleys) ≈ 0.5^5 * (1-0.5^5)
end