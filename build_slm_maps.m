function maps = build_slm_maps(traces,stim_starts,map_index,trial_length)

num_traces = size(traces,1);
% num_stims = sum(diff(stim) == 1);
% stim_starts = find(diff(stim) == 1);
num_stims = length(stim_starts);

maps = cell(num_traces,1);

for i = 1:num_stims
    
    trial_start = stim_starts(i) - 100;
    trial_end = trial_start + trial_length;
    
    for j = 1:num_traces
        
        if i == 1
            maps{j} = cell(max(map_index(:,1)),max(map_index(:,2)));
        end
        
        maps{j}{map_index(i,1),map_index(i,2)} = ...
            [maps{j}{map_index(i,1),map_index(i,2)}; traces(j,trial_start:trial_end)];
    end
end
