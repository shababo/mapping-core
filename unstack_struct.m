function results_grid = unstack_struct(results,rebuild_map)

results_grid = cell(max(rebuild_map(:,1)),max(rebuild_map(:,2)));

for trace_ind = 1:size(results,2)
    grid_inds = rebuild_map(trace_ind,:);
    results_grid{grid_inds(1),grid_inds(2)} = [results_grid{grid_inds(1),grid_inds(2)} results(trace_ind)];
end