function stacks = build_stim_stack_multi(traces,stim_starts,trial_length)

num_traces = size(traces,1);
% num_stims = sum(diff(stim) == 1);
% stim_starts = find(diff(stim) == 1);
num_stims = length(stim_starts);

stacks = cell(num_traces,1);

for i = 1:num_stims
    
    trial_start = stim_starts(i);
    trial_end = trial_start + trial_length;
    
    for j = 1:num_traces
        
        if i == 1
            stacks{j} = zeros(num_stims,trial_length+1);
        end

        stacks{j}(i,:) = traces(j,trial_start:trial_end);
    end
end
