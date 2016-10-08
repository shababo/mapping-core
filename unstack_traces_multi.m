function trace_grid = unstack_traces_multi(traces, trial_locations, rebuild_map)

trace_grid = cell(max(rebuild_map(:,1)),max(rebuild_map(:,2)));
n_sources = size(trial_locations,2);
for trace_ind = 1:size(traces,1)
    this_trial = trial_locations(trace_ind,:);
    for ns = 1:n_sources
        grid_inds = rebuild_map(this_trial(ns),:);
        trace_grid{grid_inds(1),grid_inds(2)} = [trace_grid{grid_inds(1),grid_inds(2)}; traces(trace_ind,:)];
    end
end