using LightGraphs
include("./src/load_network.jl")
include("./src/misc.jl")
include("./src/frank_wolfe.jl")

ARGSNUM = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),ARGS)

ta_data = load_ta_network(ARGS[1],ARGS[2])
link_flow, link_travel_time, objective = ta_frank_wolfe(ta_data, method=:cfw, step=:newton, log=:off, tol=ARGSNUM[3], max_iter_no=ARGSNUM[4])

#flow = sum(link_flow)
#traveltime = sum(link_travel_time)
#car_hours = dot(link_travel_time, link_flow)
#car_distance = dot(ta_data.link_length,link_flow)
flow = 0.0
for i = 1:length(ta_data.capacity)
    if ta_data.link_type[i] != 9
        if ta_data.capacity[i] != 0
            flow += link_flow[i]
        end
    end
end

traveltime = 0.0
for i = 1:length(ta_data.capacity)
    if ta_data.link_type[i] != 9
        if ta_data.capacity[i] != 0
            traveltime += link_travel_time[i]
        end
    end
end

car_hours = 0.0
for i = 1:length(ta_data.capacity)
    if ta_data.link_type[i] != 9
        if ta_data.capacity[i] != 0
            car_hours += link_travel_time[i] * link_flow[i]
        end
    end
end

car_distance = 0.0
for i = 1:length(ta_data.capacity)
    if ta_data.link_type[i] != 9
        if ta_data.capacity[i] != 0
            car_distance += ta_data.link_length[i] * link_flow[i]
        end
    end
end

println(flow," ",traveltime," ",car_hours," ",car_distance)
