function results_grid = unstack_results(results,rebuild_map,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    do_posterior = varargin{1};
else
    do_posterior = 1;
end

results_grid = cell(max(rebuild_map(:,1)),max(rebuild_map(:,2)));

for trace_ind = 1:size(results,2)
    grid_inds = rebuild_map(trace_ind,:);
    if do_posterior
        results_grid{grid_inds(1),grid_inds(2)} = [results_grid{grid_inds(1),grid_inds(2)} results(trace_ind).trials];
    else
        results_grid{grid_inds(1),grid_inds(2)} = [results_grid{grid_inds(1),grid_inds(2)} results(trace_ind)];
    end
end