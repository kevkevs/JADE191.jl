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

function is_out_of_bounds(vertex, other_vertex, bounding_longitude_data, bounding_latitude_data, city)
    longitude = city.junctions[other_vertex].longitude
    latitude = city.junctions[other_vertex].latitude

    longitude_lower_index = nothing
    for i in 1:(length(bounding_longitude_data) - 1)
        if bounding_longitude_data[i] <= longitude < bounding_longitude_data[i + 1]
            longitude_lower_index = i
            break
        end
    end

    slope = (
        bounding_latitude_data[longitude_lower_index + 1] -
        bounding_latitude_data[longitude_lower_index]
    )
    slope /= (
        bounding_longitude_data[longitude_lower_index + 1] -
        bounding_longitude_data[longitude_lower_index]
    )

    latitude_interpolation = slope * (
        longitude - bounding_longitude_data[longitude_lower_index]
    ) + bounding_latitude_data[longitude_lower_index]

    if latitude < latitude_interpolation
        return true
    end

    c_longitude = city.junctions[vertex].longitude
    if c_longitude < longitude
        return false
    end

    return false
end

"""
    getNextVertex(vertex, duration_remaining, visited, adjacencyMatrix)

Compute the next vertex in a given walk.

Returns a Tuple (N, T) where both elements are nothing if the function failed
to find a street which is crossable in the amount of time left. If the
function is able to find such a street, then N is the next vector in the path
and T is the amount of time it takes to travel from the current vertex to N.
"""
function getNextVertex(vertex, duration_remaining, visited, adjacencyMatrix, bounding_longitude_data, bounding_latitude_data, city, car)
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

        if is_out_of_bounds(vertex, other_vertex, bounding_longitude_data, bounding_latitude_data, city)
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
        if is_out_of_bounds(vertex, other_vertex, bounding_longitude_data, bounding_latitude_data, city)
            continue
        end

        if duration <= duration_remaining
            return (other_vertex, duration)
        end
    end

    return (nothing, nothing)
end

function compute_total_distance_traveled(adjacencyMatrix, path)
    total_distance = 0
    for i in 1:(length(path) - 1)
        start_vertex = path[i]
        next_vertex = path[i + 1]

        for data in adjacencyMatrix[start_vertex]
            other_vertex, duration, distance = data
            if other_vertex == next_vertex
                total_distance += distance
                break
            end
        end
    end

    return total_distance
end

"""
    getSinglePath!(adjacencyMatrix, start_vertex, total_duration, visited)

Compute a path in the city which can be traversed in total_duration time.

Returns a Vector representing the path.
"""
function getSinglePath!(adjacencyMatrix, start_vertex, total_duration, visited, bounding_longitude_data, bounding_latitude_data, city, car)
    output = [start_vertex]
    vertex = start_vertex
    duration_remaining = total_duration
    while duration_remaining > 0
        res = getNextVertex(vertex, duration_remaining, visited, adjacencyMatrix, bounding_longitude_data, bounding_latitude_data, city, car)
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

    # println("Total distance traveled: ", compute_total_distance_traveled(adjacencyMatrix, output))
    return output
end

function calculate_new_distance(output_so_far, new_path, car, city)
    default_value = new_path
    if car != 1
        default_value = output_so_far[1]
    end

    new_output = []
    for i in 1:(car - 1)
        push!(new_output, output_so_far[i])
    end

    push!(new_output, new_path)

    for _ in (car + 1):8
        push!(new_output, default_value)
    end

    return total_distance(Solution(new_output), city)
end

"""
    findExtremelyNaiveSolution(start_vertex, total_duration, adjacencyMatrix)

Compute paths in the city for all 8 cars.

Returns a Vector containing all the paths. We compute the path using a greedy
algorithm, where at each point of the path we aim to travel as much unique
new distance as possible.
"""
function findExtremelyNaiveSolution(
    start_vertex::Int64, total_duration::Int64, adjacencyMatrix,
    bounding_longitude_data, bounding_latitude_data, city
)
    output = []
    visited = Set{Tuple{Int64,Int64}}()
    NUM_CARS = 8
    prev_distance = 0
    expected_diff_per_car = 200_000
    still_apply_expectation_max_value = 800_000
    MAX_NUM_RETRIES = 30
    for car in 1:NUM_CARS
        potential_new_path = getSinglePath!(adjacencyMatrix, start_vertex, total_duration, visited, bounding_longitude_data, bounding_latitude_data, city, car)
        new_total_distance = calculate_new_distance(output, potential_new_path, car, city)
        max_new_total_distance = new_total_distance
        max_new_path = potential_new_path
        num_retries = 0
        while prev_distance < still_apply_expectation_max_value && (max_new_total_distance - prev_distance < expected_diff_per_car) && num_retries < MAX_NUM_RETRIES
            potential_new_path = getSinglePath!(adjacencyMatrix, start_vertex, total_duration, visited, bounding_longitude_data, bounding_latitude_data, city, car)
            new_total_distance = calculate_new_distance(output, potential_new_path, car, city)
            if new_total_distance > max_new_total_distance
                max_new_total_distance = new_total_distance
                max_new_path = potential_new_path
            end
            num_retries += 1
        end

        prev_distance = max_new_total_distance
        println("Distance from cars 1-", car, ": ", prev_distance)
        push!(output, max_new_path)
    end

    for car in 1:NUM_CARS
        distance_traveled_by_car = calculate_new_distance([], output[car], 1, city)
        println("Distance traveled by car ", car, ": ", distance_traveled_by_car)
    end

    for _ in (NUM_CARS + 1):8
        push!(output, output[1])
    end

    return output
end

"""
    CityWalk(city, adjacencyMatrix)

Compute paths in the city for all 8 cars and return a solution.

Returns a Solution object representing the computed paths.
"""
function CityWalk(city, adjacencyMatrix, bounding_longitude_data, bounding_latitude_data)
    solution = findExtremelyNaiveSolution(
        city.starting_junction, city.total_duration, adjacencyMatrix,
        bounding_longitude_data, bounding_latitude_data, city
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

    bounding_longitude_data, bounding_latitude_data = read_from_file("src/dividing-line-coordinates.txt")
    solution = CityWalk(city, adjacencyMatrix, bounding_longitude_data, bounding_latitude_data)

    println("Solution is feasible: ", is_feasible(solution, city))
    println("Distance covered by solution: ", total_distance(solution, city))
    # println(
    #     "Distance covered by default random walk: ", total_distance(random_walk(city), city)
    # )

    semi_random_walk_dir = "found-solutions/semi-random-walk/"
    solution_path = string(semi_random_walk_dir, "most-recent-semi-random.txt")
    plot_path = string(semi_random_walk_dir, "plots/most-recent-semi-random-plot.html")
    write_solution(solution, solution_path)
    plot_streets(city, solution; path=plot_path)

    return nothing
end

end
