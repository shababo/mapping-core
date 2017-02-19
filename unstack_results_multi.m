function results_grid = unstack_results_multi(results,map_index,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    do_posterior = varargin{1};
else
    do_posterior = 1;
end

results_grid = cell(max(max(map_index(:,1,:))),max(max(map_index(:,2,:))));

num_trials = length(results);

for i = 1:num_trials
    
    if do_posterior
        this_result = results(i).trials;
    else
        this_result = results(i);
    end
    
    
    for k = 1:size(map_index,1)
        results_grid{map_index(k,1,i),map_index(k,2,i)} = ...
            [results_grid{map_index(k,1,i),map_index(k,2,i)}; this_result];
    end
end
