function [trace_grid1,trace_grid2] = get_maps(data,trial)

this_seq = data.trial_metadata(trial).sequence;
stim_key = data.trial_metadata(trial).stim_key;
[trace1_stack,trace2_stack] = ...
    get_stim_stack(data,trial,...
    length(this_seq));

trace_grid1 = cell(17,17);
trace_grid2 = cell(17,17);
for i = 1:289
    inds = stim_key(i,[1 2])/20 + 9;
    trace_grid1{inds(1),inds(2)} = trace1_stack([this_seq.precomputed_target_index] == i,:);
    trace_grid2{inds(1),inds(2)} = trace2_stack([this_seq.precomputed_target_index] == i,:);
end