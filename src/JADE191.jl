module JADE191

using HashCode2014
using PlotlyJS


# implement Dijkstra's algorithm to find the shortest path in time from the starting junction to the goal junction
function Dijkstra_longest_path(adjacencyMatrix, city::City, start::Int, goal_func, time::Int)
    # initialize the distances to all junctions to be infinity
    num_junctions = length(city.junctions)
    distances = fill(Inf, num_junctions)

    # initialize the distances to the starting junction to be 0
    distances[start] = 0

    # initialize the unvisited set to be the set of all junctions
    unvisited = Set{Int}(1:num_junctions)

    # initialize the visited set to be empty
    visited = Set{Int}()

    # initialize the previous junction for each junction to be 0
    previous = zeros(Int, num_junctions)

    # while there are unvisited junctions
    while !isempty(unvisited)
        # find the unvisited junction with the smallest distance
        current_junction = Inf
        current = 1
        for junction in unvisited
            if distances[junction] < current_junction
                current_junction = distances[junction]
                current = junction
            end
        end


        # if the current junction is the goal junction, return the distance
        if goal_func(current)
            # initialize the path to be empty
            path = Vector{Int}()

            # while the current junction is not the starting junction
            while current != start
                # add the current junction to the path
                pushfirst!(path, current)

                # set the current junction to be the previous junction in the path
                current = previous[current]
            end

            pushfirst!(path, start)
            # return the path
            return path
        end

        # move the current junction from the unvisited set to the visited set
        delete!(unvisited, current)
        push!(visited, current)

        # for each neighbor of the current junction
        for data in adjacencyMatrix[current]
            # if the neighbor is not in the visited set

            neighbor, street_duration, _ = data
            if !(neighbor in visited)

                # if the street can be traversed in the given amount of time
                if street_duration <= time
                    # update the distance to the neighbor using the street travel time as the edge weight
                    if distances[neighbor] > distances[current] + street_duration
                        distances[neighbor] = distances[current] + street_duration
                        previous[neighbor] = current
                    end
                end
            end
        end
    end

    # if the goal junction was not reached, return 0
    return []
end


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

function is_above_line(vertex, bounding_longitude_data, bounding_latitude_data, city)
    longitude = city.junctions[vertex].longitude
    latitude = city.junctions[vertex].latitude

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

    return false
end

function is_out_of_bounds(vertex, other_vertex, bounding_longitude_data, bounding_latitude_data, city)
    if is_above_line(other_vertex, bounding_longitude_data, bounding_latitude_data, city)
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
function getNextVertex(vertex, duration_remaining, visited, new_visited, adjacencyMatrix, bounding_longitude_data, bounding_latitude_data, city, car)
    current_max_distance = -Inf
    current_next_vertex = nothing
    current_duration = nothing
    neighbor_data = get(adjacencyMatrix, vertex, nothing)
    for movement_data in neighbor_data
        other_vertex, duration, distance = movement_data
        if (
            duration > duration_remaining ||
            (vertex, other_vertex) in visited ||
            (other_vertex, vertex) in visited ||
            (vertex, other_vertex) in new_visited ||
            (other_vertex, vertex) in new_visited
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

    new_visited = Set{Tuple{Int64,Int64}}()
    while duration_remaining > 0
        res = getNextVertex(vertex, duration_remaining, visited, new_visited, adjacencyMatrix, bounding_longitude_data, bounding_latitude_data, city, car)
        next_vertex, travel_time = res
        if isnothing(next_vertex)
            break
        end

        push!(new_visited, (vertex, next_vertex))
        push!(new_visited, (next_vertex, vertex))

        push!(output, next_vertex)

        vertex = next_vertex
        duration_remaining -= travel_time
    end

    return output, new_visited
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

function plot_path(p, car, path, city)
    for vertex in path
        extendtraces!(
            p,
            Dict(
                :lat=>Vector[[city.junctions[vertex].latitude]],
                :lon=>Vector[[city.junctions[vertex].longitude]]
            ),
            [car],
            -1
        )

        sleep(0.01)
    end
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

    NUM_CARS = 4

    paris_layout = get_paris_layout()
    traces = [get_car_trace(i) for i in 1:NUM_CARS]
    p = plot(
        traces,
        paris_layout
    )
    display(p)

    expected_diff_per_car = 200_000
    still_apply_expectation_max_value = 800_000
    current_max_paths = nothing
    current_max_dist = -Inf
    NUM_TRIALS = 5
    MAX_ATTEMPTS = 30
    for trial_num in 1:NUM_TRIALS
        trial_paths = []
        prev_distance = 0
        visited = Set{Tuple{Int64,Int64}}()
        for car in 1:NUM_CARS
            potential_new_path, add_to_visited = getSinglePath!(adjacencyMatrix, start_vertex, total_duration, visited, bounding_longitude_data, bounding_latitude_data, city, car)
            new_total_distance = calculate_new_distance(trial_paths, potential_new_path, car, city)

            max_new_total_distance = new_total_distance
            max_new_path = potential_new_path
            max_add_to_visited = add_to_visited

            num_attempts = 1
            while prev_distance < still_apply_expectation_max_value && (max_new_total_distance - prev_distance < expected_diff_per_car) && num_attempts < MAX_ATTEMPTS
                potential_new_path, add_to_visited = getSinglePath!(adjacencyMatrix, start_vertex, total_duration, visited, bounding_longitude_data, bounding_latitude_data, city, car)
                new_total_distance = calculate_new_distance(trial_paths, potential_new_path, car, city)
                if new_total_distance > max_new_total_distance
                    max_new_total_distance = new_total_distance
                    max_new_path = potential_new_path
                    max_add_to_visited = add_to_visited
                end
                num_attempts += 1
            end

            push!(trial_paths, max_new_path)
            prev_distance = max_new_total_distance
            for thing in max_add_to_visited
                push!(visited, thing)
            end
        end

        trial_paths_dist = calculate_new_distance(trial_paths, trial_paths[1], 5, city)
        if trial_paths_dist > current_max_dist
            current_max_paths = trial_paths
            current_max_dist = trial_paths_dist
        end

        println("trial_num: ", trial_num)
        println("trial distance: ", trial_paths_dist)
    end

    output = current_max_paths

    println(
        "Total distance traveled: ",
        calculate_new_distance(output , output[1], 5, city)
    )

    for car in 1:NUM_CARS
        path = output[car]
        distance_traveled_by_car = calculate_new_distance([], path, 1, city)
        println("Distance traveled by car ", car, ": ", distance_traveled_by_car)
        plot_path(p, car, path, city)
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

function get_paris_layout()
    PARIS_LAT = 48.859716
    PARIS_LON = 2.349014
    layout = Layout(
                    mapbox=attr(
                        style="open-street-map",
                        center=attr(lat=PARIS_LAT, lon=PARIS_LON),
                        zoom=10,
                    ),
                    autosize=true,
                    margin=attr(l=0, r=0, t=0, b=0),
                    showlegend=true,
                    legend=attr(
                        x=0, y=1,
                        font=attr(
                            family="Courier",
                            size=12,
                            color="black"
                        )
                    )
    )

    return layout
end

function get_car_trace(car)
    trace = scattermapbox(
        ;mode="lines",
        marker=attr(size=2),
        name="car " * string(car),
        lat=[],
        lon=[]
    )

    return trace
end

function main()
    city = read_city()

    adjacencyMatrix = buildAdjacencyMatrix(city)

    bounding_longitude_data, bounding_latitude_data = read_from_file("src/dividing-line-coordinates.txt")

    function make_goal_func(bounding_longitude_data, bounding_latitude_data, city)
        function goal_func(vertex)
            return is_above_line(vertex, bounding_longitude_data, bounding_latitude_data, city)
        end

        return goal_func
    end

    goal_func = make_goal_func(bounding_longitude_data, bounding_latitude_data, city)


    paris_layout = get_paris_layout()
    traces = [get_car_trace(i) for i in 1:8]
    p = plot(
        traces,
        paris_layout
    )
    display(p)


    adjacencyMatrix = buildAdjacencyMatrix(city)
    path = Dijkstra_longest_path(adjacencyMatrix, city, city.starting_junction, goal_func, city.total_duration)
    println(path)

    plot_path(p, 1, path, city)

    return nothing

    solution = CityWalk(city, adjacencyMatrix, bounding_longitude_data, bounding_latitude_data)

    println("Solution is feasible: ", is_feasible(solution, city))
    println("Distance covered by solution: ", total_distance(solution, city))
    # println(
    #     "Distance covered by default random walk: ", total_distance(random_walk(city), city)
    # )

    semi_random_walk_dir = "found-solutions/semi-random-walk/"
    solution_path = string(semi_random_walk_dir, "most-recent-semi-random.txt")
    path_to_save_plot = string(semi_random_walk_dir, "plots/most-recent-semi-random-plot.html")
    write_solution(solution, solution_path)
    plot_streets(city, solution; path=path_to_save_plot)

    return nothing
end

end
