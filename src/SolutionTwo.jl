using HashCode2014
using DataStructures
using Distributed
using ConcaveHull
using Plots


global_visited = Set{Tuple{Int64,Int64}}()

# define a function to calculate the euclidean distance between two junctions
function distance(junction1::Junction, junction2::Junction)
    dx = junction1.latitude - junction2.latitude
    dy = junction1.longitude - junction2.longitude
    return sqrt(dx * dx + dy * dy)
end

function pick_junctions_on_hull(junctions::Vector{Junction}, n_cars::Int)
    # Compute the convex hull of the junctions
    junctions_dict = Dict([j.latitude, j.longitude] => j for j in junctions)
    hull = concave_hull(collect(keys(junctions_dict)))

    #scatter(x,y,ms=1,label="",axis=false,grid=false,markerstrokewidth=0.0)
    

    hull = hull.vertices

    # Compute the total length of the hull
    total_length = sum(distance(junctions_dict[hull[i]], junctions_dict[hull[i+1]]) for i in 1:length(hull) - 1)
    #println(distance(junctions_dict[hull[1]], junctions_dict[hull[2]]))
    #println(total_length)

    # Compute the spacing between the junctions on the hull
    spacing = 8 * total_length / (n_cars + 1)

    # Create a list to store the selected junctions
    selected = Vector{Junction}()
    x = []
    y = []

    

    # Create a variable to store the current position on the hull
    current_pos = 0

    prev_junction = hull[1]

    # Iterate over the junctions on the hull
    for junction in hull
        # Update the current position on the hull
        current_pos += distance(junctions_dict[junction], junctions_dict[prev_junction])

        # If the current position is greater than the desired spacing
        if current_pos >= spacing
            # Add the junction to the list of selected junctions
            push!(selected, junctions_dict[junction])
            append!(x, junction[1])
            append!(y, junction[2])

            # Update the current position on the hull
            current_pos = 0
            prev_junction = junction
        end
    end

    println("hull points: ", length(x))
    #p = scatter(x,y,ms=1,label="",axis=false,grid=false,markerstrokewidth=0.0)
    #plot!(p, hull)
    #display(p)

    # Return the list of selected junctions
    return selected
    
    # Find the incenter and inradius of the concave hull
    #incenter = incenter(hull)
    #inradius = inradius(hull)

    # Generate n_cars points on the circumference of the inscribed circle
    #circumference_points = circumference(incenter, inradius, n_cars)

    #return [ for i in circumference_points]
end


# implement Dijkstra's algorithm to find the longest path from the starting junction to the goal junction
function Dijkstra_longest_path(graph, city::City, start::Int, goal::Int, time::Int)
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
        current = 0
        for junction in unvisited
            if distances[junction] < current_junction
                current_junction = distances[junction]
                current = junction
            end
        end


        # if the current junction is the goal junction, return the distance
        if current == goal
            # initialize the path to be empty
            path = Vector{Int}()

            # initialize the current junction to be the goal junction
            current = goal

            # while the current junction is not the starting junction
            while current != start
                # add the current junction to the path
                pushfirst!(path, current)

                # set the current junction to be the previous junction in the path
                current = previous[current]
            end

            # add the starting junction to the path
            pushfirst!(path, start)
            total_time = 0
            for i in 1:length(path)-1
                total_time += graph[city.junctions[path[i]]][city.junctions[path[i+1]]].duration
                push!(global_visited, (path[i], path[i+1]))
                push!(global_visited, (path[i+1], path[i]))
            end
            #println("duration: ", total_time)

            # return the path and the distance
            return path, total_time
        end

        # move the current junction from the unvisited set to the visited set
        delete!(unvisited, current)
        push!(visited, current)

        # for each neighbor of the current junction
        for (neighbor_junc, street) in graph[city.junctions[current]]
            # if the neighbor is not in the visited set
            neighbor = street.endpointB
            if !(neighbor in visited)

                # if the street can be traversed in the given amount of time
                if street.duration <= time
                    # update the distance to the neighbor using the street length as the edge weight
                    if distances[neighbor] > distances[current] + street.duration
                        distances[neighbor] = distances[current] + street.duration
                        previous[neighbor] = current
                    end
                end
            end
        end
    end

    # if the goal junction was not reached, return 0
    return [], time
end

"""
    getNextVertex(vertex, duration_remaining, visited, adjacencyMatrix)

Compute the next vertex in a given walk.

Returns a Tuple (N, T) where both elements are nothing if the function failed
to find a street which is crossable in the amount of time left. If the
function is able to find such a street, then N is the next vector in the path
and T is the amount of time it takes to travel from the current vertex to N.
"""
function getNextVertex(city, vertex, duration_remaining, visited, graph)
    current_max_distance = -Inf
    current_next_vertex = nothing
    current_duration = nothing
    #neighbor_data = get(adjacencyMatrix, vertex, nothing)
    vertex_key = city.junctions[vertex]
    for (neighbor, street) in graph[vertex_key]
        other_vertex = street.endpointB
        if other_vertex == vertex
            other_vertex = street.endpointA
        end
        duration = street.duration
        distance = street.distance
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
        neighbor = rand(collect(keys(graph[vertex_key])))
        street = graph[vertex_key][neighbor]
        other_vertex = street.endpointB
        if other_vertex == vertex
            other_vertex = street.endpointA
        end
        duration = street.duration
        
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
function getSinglePath!(graph, city, start_vertex, total_duration, visited)
    output = [start_vertex]
    vertex = start_vertex
    duration_remaining = total_duration
    while duration_remaining > 0
        res = getNextVertex(city, vertex, duration_remaining, visited, graph)
        next_vertex, travel_time = res
        if next_vertex == vertex
            #println(vertex, " ")
        end
        if isnothing(next_vertex)
            break
        end

        push!(visited, (vertex, next_vertex))
        #push!(visited, (next_vertex, vertex))

        push!(output, next_vertex)

        vertex = next_vertex
        duration_remaining -= travel_time
    end

    return output
end



function parallel_Dijkstra_longest_path(graph, city::City, start::Int, time::Int, n_cars::Int, hull::Vector{Junction})
    tasks = [Dijkstra_longest_path(graph, city, start, findfirst(isequal(hull[i]), city.junctions), time) for i in 1:n_cars]

    solution = Solution([i[1] for i in tasks])

    println("half Solution is feasible: ", is_feasible(solution, city))
    println("half Distance covered by solution: ", total_distance(solution, city))

    semi_random_walk_dir = "found-solutions/semi-random-walk/"
    solution_path = string(semi_random_walk_dir, "half-most-recent-semi-random.txt")
    plot_path = string(semi_random_walk_dir, "plots/half-most-recent-semi-random-plot.html")
    write_solution(solution, solution_path)
    plot_streets(city, solution; path=plot_path)

    output = []
    visited = global_visited
    for (path, t) in tasks
        single_path = getSinglePath!(graph, city, path[end], city.total_duration - t, visited)[2:end]
        new_path = [path; single_path]
                
        total_time = 0
        for i in 1:length(new_path)-1
            a = new_path[i]
            b = new_path[i+1]
            s = graph[city.junctions[a]][city.junctions[b]]
            total_time += s.duration
            if !HashCode2014.is_street(a, b, s)
                println("bad street: ", a, " -> ", b)
            end
        end
        if total_time > city.total_duration
            println("over time")
        end
        

        push!(output, new_path)
    end

    # Wait for all tasks to complete and return the results
    return output
end

function make_graph(junctions::Vector{Junction}, streets::Vector{Street})
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
    return graph
end


city = read_city()
junctions = city.junctions
streets = city.streets
n_cars = city.nb_cars
time_limit = city.total_duration
start = city.starting_junction

hull = pick_junctions_on_hull(junctions, n_cars)
graph = make_graph(junctions, streets)

solution = parallel_Dijkstra_longest_path(graph, city, start, time_limit, n_cars, hull)
#println(solution)

#results = parallel_max_street_length(graph, junctions[start], time_limit, n_cars, hull, junctions)
#total_length = sum(results)
#println(results)


solution = Solution(solution)

println("Solution is feasible: ", is_feasible(solution, city))
println("Distance covered by solution: ", total_distance(solution, city))

semi_random_walk_dir = "found-solutions/semi-random-walk/"
solution_path = string(semi_random_walk_dir, "most-recent-semi-random.txt")
plot_path = string(semi_random_walk_dir, "plots/most-recent-semi-random-plot.html")
write_solution(solution, solution_path)
plot_streets(city, solution; path=plot_path)

