using Aqua
using Documenter
using JADE191
using JuliaFormatter
using Test

DocMeta.setdocmeta!(
    JADE191,
    :DocTestSetup,
    :(using JADE191);
    recursive=true
)

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
end
