function [traces_ch1, traces_ch2, full_seq, traces_ch3, full_stim_key] = ...
    get_traces(data,trials,varargin)

this_seq = cell(length(trials),1);
stim_starts = cell(length(trials),1);
% full_seq = [];
stims_per_trial = zeros(length(trials),1);
full_stim_key = [];
for i = 1:length(trials)
    
    cur_trial = trials(i);
    this_seq{i} = data.trial_metadata(cur_trial).sequence;
    stim_starts{i} =  20*[data.trial_metadata(cur_trial).sequence.start]; % MAGIC NUMBER  %20*[data.sequence.start]
    stims_per_trial(i) = length(stim_starts{i});
    if isfield(this_seq{i},'group')
        this_seq{i} = rmfield(this_seq{i},'group');
    end
    for j = 1:length(this_seq{i})
        
        if i == 1 && j == 1
            full_seq(1) = this_seq{i}(j);
        else
            full_seq(end+1) = this_seq{i}(j);
        end
        full_seq(end).precomputed_target_index = ...
            full_seq(end).precomputed_target_index + size(full_stim_key,1);
    end
    full_stim_key = [full_stim_key; data.trial_metadata(cur_trial).stim_key];
end

% max_trial = 1200;
[traces_ch1,traces_ch2, traces_ch3] = ...
    get_stim_stack(data,trials,...
        stims_per_trial,stim_starts,varargin{:});
