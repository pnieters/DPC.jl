using ADSP, Test

@testset "All unit tests" begin
    @testset "YAML-loading tests" begin
        include("loading_tests.jl")
    end

    @testset "Functionality tests" begin
        include("running_tests.jl")
    end

    @testset "Other utility functions tests" begin
        include("utils_tests.jl")
    end
end