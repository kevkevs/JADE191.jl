using HashCode2014
using DataStructures
using Distributed
using ConcaveHull

# define a function to calculate the euclidean distance between two junctions
function distance(junction1::Junction, junction2::Junction)
    dx = junction1.latitude - junction2.latitude
    dy = junction1.longitude - junction2.longitude
    return sqrt(dx * dx + dy * dy)
end

function pick_junctions_on_hull(junctions::Vector{Junction}, n_cars::Int)
    # Compute the convex hull of the junctions
    junctions_dict = Dict([j.latitude, j.longitude] => j for j in junctions)
    hull = concave_hull(collect(keys(junctions_dict))).vertices

    # Compute the total length of the hull
    total_length = sum(distance(junctions_dict[hull[i]], junctions_dict[hull[i+1]]) for i in 1:length(hull) - 1)

    # Compute the spacing between the junctions on the hull
    spacing = total_length / (n_cars + 1)

    # Create a list to store the selected junctions
    selected = Vector{Junction}()

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

            # Update the current position on the hull
            current_pos -= spacing
        end
    end

    # Return the list of selected junctions
    return selected
end

# https://stackoverflow.com/questions/48104390/julias-most-efficient-way-to-choose-longest-array-in-array-of-arrays
function findlongest(A)
    idx = 0
    len = 0
    @inbounds for i in 1:length(A)
        l = length(A[i])
        l > len && (idx = i; len=l)
    end
    return A[idx]
end

function max_street_length(graph, start_junction::Junction, duration::Int, goal::Junction, junctions::Vector{Junction}; visited = Set())
    # Create a queue to store the junctions to visit lowest, priority = closest to heuristic
    queue = PriorityQueue{Vector{Junction}, Float64}()
    enqueue!(queue, [start_junction], distance(start_junction, goal))

    # Create a list to store the paths taken by each car
    paths = []

    # Create a variable to store the total length of streets traversed
    total_length = 0

    # Create a variable to store the remaining time
    remaining_time = duration

    # While there are junctions to visit and there is still time remaining
    while !isempty(queue) && remaining_time > 0
        # Get the next junction to visit from the queue
        curr_path = dequeue!(queue)
        println(typeof(curr_path))
        junction = curr_path[end]

        # If the junction has not been visited yet
        if !(junction in visited)
            # Mark the junction as visited
            push!(visited, junction)

            if junction == goal
                push!(paths, curr_path)
            end

            # Iterate over all streets connected to the junction
            for (junction, street) in graph[junction]
                # Calculate the time it takes to traverse the street
                t = street.duration

                # If there is enough time remaining to traverse the street
                if remaining_time >= t
                    # Update the total length of streets traversed
                    total_length += street.distance

                    # Update the remaining time
                    remaining_time -= t

                    found_path_len, found_path = max_street_length(graph,junction, remaining_time, goal, junctions)

                    # Add the junction at the other end of the street to the queue, using the heuristic to prioritize junctions closer to the goal
                    next_route = [found_path..., junctions[street.endpointB]]
                    enqueue!(queue, next_route, distance(junctions[street.endpointB], goal))
                end
            end
        end
    end


    longest_path = findlongest(paths)

    # Return the total length of streets traversed
    return total_length, longest_path
end

function parallel_max_street_length(graph, start_junction::Junction, duration::Int, n_cars::Int, hull::Vector{Junction}, junctions::Vector{Junction})
    # Create a list of tasks to compute the maximum street length for each car
    tasks = [max_street_length(graph, start_junction, duration, hull[i], junctions)[2] for i in 1:n_cars]

    # Wait for all tasks to complete and return the results
    return tasks
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

results = parallel_max_street_length(graph, junctions[start], time_limit, n_cars, hull, junctions)
#total_length = sum(results)
#println(results)
