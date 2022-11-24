module JADE191

using HashCode2014

function buildAdjacencyMatrix(city)
    adjacencyMatrix = Dict()

    for street in city.streets
        push!(
            get!(adjacencyMatrix, street.endpointA, []),
            (street.endpointB, street.duration, street.distance)
        )

        if !street.bidirectional
            continue
        end

        push!(
            get!(adjacencyMatrix, street.endpointB, []),
            (street.endpointA, street.duration, street.distance)
        )
    end

    return adjacencyMatrix
end

function getNextVertex(vertex, duration_remaining, visited, adjacencyMatrix)
    current_max_distance = -Inf
    current_next_vertex = nothing
    current_duration = nothing
    neighbor_data = get(adjacencyMatrix, vertex, nothing)
    for movement_data in neighbor_data
        other_vertex, duration, distance = movement_data
        if (duration > duration_remaining ||
            (vertex, other_vertex) in visited ||
            (other_vertex, vertex) in visited)
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

function getSinglePath!(adjacencyMatrix, start_vertex, total_duration, visited)
    output = [start_vertex]
    vertex = start_vertex
    duration_remaining = total_duration
    while duration_remaining > 0
        res = getNextVertex(vertex, duration_remaining, visited,
                            adjacencyMatrix)
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

function findExtremelyNaiveSolution(start_vertex::Int,
                                    total_duration::Int,
                                    adjacencyMatrix)
    output = []
    visited = []
    for car in 1:8
        push!(
            output,
            getSinglePath!(
                adjacencyMatrix, start_vertex, total_duration, visited
            )
        )
    end

    return output
end

function CityWalk(city, adjacencyMatrix)
    solution = findExtremelyNaiveSolution(
        city.starting_junction,
        city.total_duration,
        adjacencyMatrix
    )

    solution = Solution(solution)

    return solution
end

function main()
    city = read_city()
    adjacencyMatrix = buildAdjacencyMatrix(city)

    solution = CityWalk(city, adjacencyMatrix)

    println("Solution is feasible: ", is_feasible(solution, city))
    println("Distance covered by solution: ", total_distance(solution, city))
    println(
        "Distance covered by default random walk: ",
        total_distance(random_walk(city), city)
    )

    semi_random_walk_dir = "found-solutions/semi-random-walk/"
    solution_path = semi_random_walk_dir + "most-recent-semi-random.txt"
    plot_path = (
        semi_random_walk_dir +
        "plots/most-recent-semi-random-plot.html"
    )
    write_solution(solution, solution_path)
    plot_streets(city, solution; path=plot_path)

    return solution
end

end
