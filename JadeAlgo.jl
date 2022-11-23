using HashCode2014

city = read_city()

println(city.starting_junction)
println(city.junctions[city.starting_junction])

cycles_found = nothing

#=

# Python program to detect cycle
# in a graph
 
from collections import defaultdict
 
class Graph():
    def __init__(self,vertices):
        self.graph = defaultdict(list)
        self.V = vertices
 
    def addEdge(self,u,v):
        self.graph[u].append(v)
 
    def isCyclicUtil(self, v, visited, recStack):
 
        # Mark current node as visited and
        # adds to recursion stack
        visited[v] = True
        recStack[v] = True
 
        # Recur for all neighbours
        # if any neighbour is visited and in
        # recStack then graph is cyclic
        for neighbour in self.graph[v]:
            if visited[neighbour] == False:
                if self.isCyclicUtil(neighbour, visited, recStack) == True:
                    return True
            elif recStack[neighbour] == True:
                return True
 
        # The node needs to be popped from
        # recursion stack before function ends
        recStack[v] = False
        return False
 
    # Returns true if graph is cyclic else false
    def isCyclic(self):
        visited = [False] * (self.V + 1)
        recStack = [False] * (self.V + 1)
        for node in range(self.V):
            if visited[node] == False:
                if self.isCyclicUtil(node,visited,recStack) == True:
                    return True
        return False
=#

# adjacency matrix: a[i,j] = street if street from junction i to junction j
adjacencyMatrix = Matrix{Union{Street, Nothing}}(undef, size(city.junctions,1), size(city.junctions,1))

for street in city.streets
    adjacencyMatrix[street.endpointA,street.endpointB] = street
    if street.bidirectional
        adjacencyMatrix[street.endpointB,street.endpointA] = street
    end
end


solutions_junctions = Dict{Int, Any}()
solutions_streets = Dict()

#= 
current implementation doesn't allow bi directional travel. A different approach might be to define some
cardinal heuristic points
=#
function findCycles(vertex::Int, visited::Array{<:Real}, stack::Array{<:Real}, duration::Int, stack_streets::Array{<:Street})
    append!(visited, [vertex])
    append!(stack, [vertex])
    #print(stack)
    for i in 1:size(adjacencyMatrix,1)
        neighbor = adjacencyMatrix[vertex,i] # neighbor is of type Street
        if !isnothing(neighbor) && neighbor.duration <= duration #cut search short if duration too long
            if neighbor.endpointB == city.starting_junction #only accepted cycle if start == end
                append!(visited, [neighbor.endpointB])
                append!(stack, [neighbor.endpointB])
                append!(stack_streets, [neighbor])
                n = size(stack,1)
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
                findCycles(neighbor.endpointB, visited, stack, duration - neighbor.duration, stack_streets)
                pop!(stack)
                pop!(stack_streets)
            end
        end
    end
    return
end

findCycles(city.starting_junction, Array{Int}(undef, 0), Array{Int}(undef, 0), city.total_duration, Array{Street}(undef, 0))

println("solutions' lengths: ", length(solutions_junctions))
println(keys(solutions_junctions))

# as a starting point it returns an itinerary where all cars have the same path
function make_itinerary()
    return itinerary = [solutions_junctions[maximum(collect(keys(solutions_junctions)))][1] for i in 1:city.nb_cars]
end

solution = Solution(make_itinerary())
println("is solution? ", is_feasible(solution, city))