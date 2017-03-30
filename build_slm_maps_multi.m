function maps = build_slm_maps_multi(traces,stim_starts,map_index,trial_length)

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
            maps{j} = cell(max(max(map_index(:,1,:))),max(max(map_index(:,2,:))));
        end
        
        for k = 1:size(map_index,1)
            try
            maps{j}{map_index(k,1,i),map_index(k,2,i)} = ...
                [maps{j}{map_index(k,1,i),map_index(k,2,i)}; traces(j,trial_start:trial_end)];
            catch
                trial_start
                trial_end
                return
            end
        end
    end
end
