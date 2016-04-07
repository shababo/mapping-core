function [results, rebuild_map] = stack_results(results)


results_grid = results;
rebuild_map = [];
results = [];
for i = 1:size(results_grid,1)
    for j = 1:size(results_grid,2)
        for k = 1:length(results_grid{i,j})
        	results = [results results_grid{i,j}(k)];
        end
        rebuild_map = [rebuild_map; repmat([i j],length(results_grid{i,j}),1)];
    end
end
