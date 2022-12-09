using HashCode2014
using DataStructures

function CityWalk()
    city = read_city()

    # adjacency matrix: a[i,j] = street if street from junction i to junction j
    adjacencyMatrix = Matrix{Union{Street,Nothing}}(
        undef, size(city.junctions, 1), size(city.junctions, 1)
    )

    for street in city.streets
        adjacencyMatrix[street.endpointA, street.endpointB] = street
        if street.bidirectional
            adjacencyMatrix[street.endpointB, street.endpointA] = street
        end
    end

    solutions_junctions = Dict{Int,Any}()
    solutions_streets = Dict()

    #=
    current implementation doesn't allow bi directional travel. A different approach might be to define some
    cardinal heuristic points
    =#
    function findCycles(
        vertex::Int,
        visited::Array{<:Real},
        stack::Array{<:Real},
        duration::Int,
        stack_streets::Array{<:Street},
    )
        append!(visited, [vertex])
        append!(stack, [vertex])
        #print(stack)
        for i in 1:size(adjacencyMatrix, 1)
            neighbor = adjacencyMatrix[vertex, i] # neighbor is of type Street
            if !isnothing(neighbor) && neighbor.duration <= duration #cut search short if duration too long
                if neighbor.endpointB == city.starting_junction #only accepted cycle if start == end
                    append!(visited, [neighbor.endpointB])
                    append!(stack, [neighbor.endpointB])
                    append!(stack_streets, [neighbor])
                    n = size(stack, 1)
                    if !haskey(solutions_junctions, n)
                        solutions_junctions[n] = []
                        solutions_streets[n] = []
                        #println("n: ", n)
                    end
                    append!(solutions_junctions[n], [stack])
                    append!(solutions_streets[n], [stack_streets])

                    pop!(stack)
                    pop!(stack_streets)
                elseif !(neighbor.endpointB in visited)
                    append!(stack_streets, [neighbor])
                    findCycles(
                        neighbor.endpointB,
                        visited,
                        stack,
                        duration - neighbor.duration,
                        stack_streets,
                    )
                    pop!(stack)
                    pop!(stack_streets)
                end
            end
        end
        return nothing
    end

    findCycles(
        city.starting_junction,
        Array{Int}(undef, 0),
        Array{Int}(undef, 0),
        city.total_duration,
        Array{Street}(undef, 0),
    )

    #println("solutions' lengths: ", length(solutions_junctions))
    #println(keys(solutions_junctions))

    # as a starting point it returns an itinerary where all cars have the same path
    function make_itinerary()
        return itinerary = [
            solutions_junctions[maximum(collect(keys(solutions_junctions)))][1] for
            i in 1:(city.nb_cars)
        ]
    end

    solution = Solution(make_itinerary())
    write_solution(solution, "solution.txt")
    return is_feasible(solution, city)
end

function find_cycle(path::Vector{T}, x::T) where T
    # Initialize a map to store the index of each element in the path
    indices = Dict{T, Int}()
    
    # Iterate over the elements in the path and store their indices
    for (i, y) in enumerate(path)
        # If the element is the value we are looking for, it means we have found a cycle
        if y == x
            # Return the portion of the path from the first occurrence of the element to the end
            return path[indices[x]:end]
        end
        # Store the index of the element in the map
        indices[y] = i
    end
    
    # If we reach this point, it means the value is not in the path or there is no cycle containing the value
    return []
end


# define a function to find the optimal routes for the cars
function find_optimal_routes(junctions::Vector{Junction}, streets::Vector{Street}, n_cars::Int64, time_limit::Int64, start::Int64)
    # create a complete graph of the city, where each node represents a junction
    # and each edge represents a street connecting two junctions
    # the weights of the edges are the time it takes to traverse the street
    graph = Dict{Junction, Dict{Junction, Street}}()
    for street in streets
        endpointA = junctions[street.endpointA]
        endpointB = junctions[street.endpointB]

        if !haskey(graph, endpointA)
            graph[endpointA] = Dict(endpointB => street)
        else
            graph[endpointA][endpointB] = street
        end

        if street.bidirectional == true
            if !haskey(graph, endpointB)
                graph[endpointB] = Dict(endpointA => street)
            else
                graph[endpointB][endpointA] = street
            end
        end
    end

    # initialize the priority queue with the starting junction as the first node
    # and a value of 0

    #found_cycles = Dict{Junction, PriorityQueue{Vector{Junction}, Float64}}
    visited = Vector{Vector{Junction}}()
    finished_routes = PriorityQueue{Vector{Junction}, Float64}(Base.Order.Reverse)
    queue = PriorityQueue{Vector{Junction}, Float64}()
    enqueue!(queue, [junctions[start]], 0)

    # initialize the set of optimal routes with an empty list
    optimal_routes = Vector{Vector{Junction}}()

    # continue until the priority queue is empty
    i = 0
    while !isempty(queue)
        print(" ", i)
        i += 1
        # pop the route with the lowest total time from the queue
        (route, time) = dequeue_pair!(queue)

        # check if the route is complete (i.e., it includes all junctions)
        if length(route) == length(junctions)
            # add the route to the list of optimal routes
            push!(optimal_routes, route)

            # continue to the next iteration if the list of optimal routes
            # is already long enough
            if length(optimal_routes) >= n_cars
                continue
            end
        end

        # get the last junction in the route
        last_junction = route[end]

        found_next_junction = false

        # loop through all the possible next steps from the current junction
        for (next_junction, next_street) in graph[last_junction]
            # calculate the time it would take to reach the next junction
            next_time = time + next_street.duration
            sub_visited = false

            if length(route) >= 3
                sub_arr = [route[length(route)-2:length(route)]..., next_junction]
                println(sub_arr)
                if sub_arr ∉ visited
                    sub_visited = true
                    push!(visited, sub_arr)
                end
            end

            # check if the route is within the time limit
            if next_time <= time_limit && !sub_visited
                found_next_junction = true

                # create a new route by appending the next junction to the current route
                next_route = [route..., next_junction]

                # add the new route to the priority queue
                enqueue!(queue, next_route, next_time)
            end
        end

        if !found_next_junction
            enqueue!(finished_routes, route, length(route))
        end
    end

    # return the list of optimal routes
    return optimal_routes
end

function main()
    city = read_city()
    junctions = city.junctions
    streets = city.streets
    n_cars = city.nb_cars
    time_limit = city.total_duration
    start = city.starting_junction
    solution = find_optimal_routes(junctions, streets, n_cars, time_limit, start)

    println(solution)

    #=
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
    =#

    return solution
end

main()