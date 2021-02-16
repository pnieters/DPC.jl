using ADSP, Test

@testset "All unit tests" begin
    include("loading_tests.jl")
    include("running_tests.jl")
end