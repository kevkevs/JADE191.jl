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
        city = read_city(input_path)
        expected = Dict(
            2 => Any[(3, 45, 200)], 3 => Any[(2, 45, 200)], 1 => Any[(2, 30, 250)]
        )
        @test JADE191.buildAdjacencyMatrix(city) == expected
    end

    @testset verbose = true "Check boundary false 1" begin
        bounding_longitude_data, bounding_latitude_data = JADE191.read_from_file(
            "../src/dividing-line-coordinates.txt"
        )

        city = read_city()
        result = JADE191.is_out_of_bounds(
            city.starting_junction, 3, bounding_longitude_data, bounding_latitude_data, city
        )
        @test result == false
    end

    @testset verbose = true "Check boundary false 2" begin
        bounding_longitude_data, bounding_latitude_data = JADE191.read_from_file(
            "../src/dividing-line-coordinates.txt"
        )

        city = read_city()
        result = JADE191.is_out_of_bounds(
            city.starting_junction,
            668,
            bounding_longitude_data,
            bounding_latitude_data,
            city,
        )
        @test result == false
    end

    @testset verbose = true "Check boundary true" begin
        bounding_longitude_data, bounding_latitude_data = JADE191.read_from_file(
            "../src/dividing-line-coordinates.txt"
        )

        city = read_city()
        result = JADE191.is_out_of_bounds(
            city.starting_junction, 2, bounding_longitude_data, bounding_latitude_data, city
        )
        @test result == true
    end

    @testset verbose = true "Check getNextVertex is None" begin
        bounding_longitude_data, bounding_latitude_data = JADE191.read_from_file(
            "../src/dividing-line-coordinates.txt"
        )

        city = read_city()
        car = 1
        vertex = 2
        duration = 0
        adjacencyMatrix = JADE191.buildAdjacencyMatrix(city)
        visited = Set{Tuple{Int64,Int64}}()
        new_visited = Set{Tuple{Int64,Int64}}()
        should_stay_above = false
        bounding_longitude_data, bounding_latitude_data = JADE191.read_from_file(
            "../src/dividing-line-coordinates.txt"
        )

        result = JADE191.getNextVertex(
            vertex,
            duration,
            visited,
            new_visited,
            adjacencyMatrix,
            bounding_longitude_data,
            bounding_latitude_data,
            city,
            should_stay_above,
        )
        @test result == (nothing, nothing)
    end

    @testset verbose = true "Check getNextVertex is None 2" begin
        bounding_longitude_data, bounding_latitude_data = JADE191.read_from_file(
            "../src/dividing-line-coordinates.txt"
        )

        city = read_city()
        car = 1
        vertex = 1
        duration = 5000
        adjacencyMatrix = JADE191.buildAdjacencyMatrix(city)
        visited = Set{Tuple{Int64,Int64}}()
        new_visited = Set{Tuple{Int64,Int64}}()
        should_stay_above = true
        bounding_longitude_data, bounding_latitude_data = JADE191.read_from_file(
            "../src/dividing-line-coordinates.txt"
        )

        result = JADE191.getNextVertex(
            vertex,
            duration,
            visited,
            new_visited,
            adjacencyMatrix,
            bounding_longitude_data,
            bounding_latitude_data,
            city,
            should_stay_above,
        )
        @test result == (nothing, nothing)
    end

    @testset verbose = true "Check getNextVertex passes" begin
        bounding_longitude_data, bounding_latitude_data = JADE191.read_from_file(
            "../src/dividing-line-coordinates.txt"
        )

        city = read_city()
        car = 1
        vertex = 4
        duration = 5000
        adjacencyMatrix = JADE191.buildAdjacencyMatrix(city)
        visited = Set{Tuple{Int64,Int64}}()
        new_visited = Set{Tuple{Int64,Int64}}()
        should_stay_above = true
        bounding_longitude_data, bounding_latitude_data = JADE191.read_from_file(
            "../src/dividing-line-coordinates.txt"
        )

        result = JADE191.getNextVertex(
            vertex,
            duration,
            visited,
            new_visited,
            adjacencyMatrix,
            bounding_longitude_data,
            bounding_latitude_data,
            city,
            should_stay_above,
        )
        @test result == (2122, 9)
    end
end
