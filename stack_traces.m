function [traces, rebuild_map] = stack_traces(traces)


trace_grid = traces;
rebuild_map = [];
traces = [];
for i = 1:size(trace_grid,1)
    for j = 1:size(trace_grid,2)
        traces = [traces; trace_grid{i,j}];
        rebuild_map = [rebuild_map; repmat([i j],size(trace_grid{i,j},1),1)];
    end
end
