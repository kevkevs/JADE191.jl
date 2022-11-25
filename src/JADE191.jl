module JADE191

using HashCode2014

"""
    buildAdjacencyMatrix(city)

Compute an adjaceny matrix to represent the connections between junctions.

Returns a Dict D where D[A] is a vector of all
junctions connected to A by a street.
"""
function buildAdjacencyMatrix(city)
    adjacencyMatrix = Dict()

    for street in city.streets
        push!(
            get!(adjacencyMatrix, street.endpointA, []),
            (street.endpointB, street.duration, street.distance),
        )

        if !street.bidirectional
            continue
        end

        push!(
            get!(adjacencyMatrix, street.endpointB, []),
            (street.endpointA, street.duration, street.distance),
        )
    end

    return adjacencyMatrix
end

"""
    getNextVertex(vertex, duration_remaining, visited, adjacencyMatrix)

Compute the next vertex in a given walk.

Returns a Tuple (N, T) where both elements are nothing if the function failed
to find a street which is crossable in the amount of time left. If the
function is able to find such a street, then N is the next vector in the path
and T is the amount of time it takes to travel from the current vertex to N.
"""
function getNextVertex(vertex, duration_remaining, visited, adjacencyMatrix)
    current_max_distance = -Inf
    current_next_vertex = nothing
    current_duration = nothing
    neighbor_data = get(adjacencyMatrix, vertex, nothing)
    for movement_data in neighbor_data
        other_vertex, duration, distance = movement_data
        if (
            duration > duration_remaining ||
            (vertex, other_vertex) in visited ||
            (other_vertex, vertex) in visited
        )
            continue
        end

        if distance > current_max_distance
            current_max_distance = distance
            current_next_vertex = other_vertex
            current_duration = duration
        end
    end

    if !isnothing(current_next_vertex)
        return (current_next_vertex, current_duration)
    end

    NUM_RANDOM_TRIALS = 10
    for _ in 1:NUM_RANDOM_TRIALS
        other_vertex, duration, _ = rand(neighbor_data)
        if duration <= duration_remaining
            return (other_vertex, duration)
        end
    end

    return (nothing, nothing)
end

"""
    getSinglePath!(adjacencyMatrix, start_vertex, total_duration, visited)

Compute a path in the city which can be traversed in total_duration time.

Returns a Vector representing the path.
"""
function getSinglePath!(adjacencyMatrix, start_vertex, total_duration, visited)
    output = [start_vertex]
    vertex = start_vertex
    duration_remaining = total_duration
    while duration_remaining > 0
        res = getNextVertex(vertex, duration_remaining, visited, adjacencyMatrix)
        next_vertex, travel_time = res
        if isnothing(next_vertex)
            break
        end

        push!(visited, (vertex, next_vertex))
        push!(visited, (next_vertex, vertex))

        push!(output, next_vertex)

        vertex = next_vertex
        duration_remaining -= travel_time
    end

    return output
end

"""
    findExtremelyNaiveSolution(start_vertex, total_duration, adjacencyMatrix)

Compute paths in the city for all 8 cars.

Returns a Vector containing all the paths. We compute the path using a greedy
algorithm, where at each point of the path we aim to travel as much unique
new distance as possible.
"""
function findExtremelyNaiveSolution(
    start_vertex::Int64, total_duration::Int64, adjacencyMatrix
)
    output = []
    visited = Set{Tuple{Int64,Int64}}()
    for car in 1:8
        push!(
            output, getSinglePath!(adjacencyMatrix, start_vertex, total_duration, visited)
        )
    end

    return output
end

"""
    CityWalk(city, adjacencyMatrix)

Compute paths in the city for all 8 cars and return a solution.

Returns a Solution object representing the computed paths.
"""
function CityWalk(city, adjacencyMatrix)
    solution = findExtremelyNaiveSolution(
        city.starting_junction, city.total_duration, adjacencyMatrix
    )

    solution = Solution(solution)

    return solution
end

function read_from_file(filename)
    xdata = Vector{Float64}()
    ydata = Vector{Float64}()
    open(filename) do f
        for line in eachline(f)
            xcor = parse(Float64, split(line, ",")[1])
            ycor = parse(Float64, split(line, ",")[2])

            push!(xdata, xcor)
            push!(ydata, ycor)

        end
    end

    return xdata, ydata
end

"""
    main()

Loads the city data, computes solutions, and stores them.
"""
function main()
    city = read_city()
    adjacencyMatrix = buildAdjacencyMatrix(city)

    solution = CityWalk(city, adjacencyMatrix)

    println("Solution is feasible: ", is_feasible(solution, city))
    println("Distance covered by solution: ", total_distance(solution, city))
    println(
        "Distance covered by default random walk: ", total_distance(random_walk(city), city)
    )

    semi_random_walk_dir = "found-solutions/semi-random-walk/"
    solution_path = string(semi_random_walk_dir, "most-recent-semi-random.txt")
    plot_path = string(semi_random_walk_dir, "plots/most-recent-semi-random-plot.html")
    write_solution(solution, solution_path)
    plot_streets(city, solution; path=plot_path)

    return solution
end

end
