using Aqua
using Documenter
using HashCode2014
using JADE191
using JuliaFormatter
using Test

DocMeta.setdocmeta!(JADE191, :DocTestSetup, :(using JADE191); recursive=true)

@testset "JADE191.jl" begin
    @testset verbose = true "Code quality (Aqua.jl)" begin
        Aqua.test_all(JADE191; ambiguities=false)
    end

    @testset verbose = true "Code formatting (JuliaFormatter.jl)" begin
        @test format(JADE191; verbose=true, overwrite=false)
    end

    @testset verbose = true "Doctests (Documenter.jl)" begin
        doctest(JADE191)
    end

    @testset verbose = true "My own tests" begin
        @test 1 + 1 == 2
    end

    @testset verbose = true "Check adjacency matrix" begin
        input_path = joinpath(@__DIR__, "data", "example_input.txt")
        output_path = joinpath(@__DIR__, "data", "example_output.txt")
        city = read_city(input_path)
        expected = Dict(2 => (3, 45, 200), 3 => (2, 45, 200), 1 => (2, 30, 250))
        @test buildAdjacencyMatrix(city) == expected
    end
end
