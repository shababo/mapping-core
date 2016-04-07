function trace_grid = unstack_traces(traces,rebuild_map)

trace_grid = cell(max(rebuild_map(:,1)),max(rebuild_map(:,2)));
for trace_ind = 1:size(traces,1)
    grid_inds = rebuild_map(trace_ind,:);
    trace_grid{grid_inds(1),grid_inds(2)} = [trace_grid{grid_inds(1),grid_inds(2)}; traces(trace_ind,:)];
end