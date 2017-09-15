function [traces_ch1, traces_ch2, full_seq] = ...
    get_traces(data,trials)



this_seq = cell(length(trials),1);
this_stim_key = cell(length(trials),1);
power_curve_num = cell(length(trials),1);
stim_starts = cell(length(trials),1);

full_stim_key = [];
% clear full_seq
% full_seq = struct();
for i = 1:length(trials)
    cur_trial = trials(i);
    this_seq{i} = data.trial_metadata(cur_trial).sequence;
    stims_per_trial(i) = length(this_seq{i});
    this_stim_key{i} = data.trial_metadata(cur_trial).stim_key;
    power_curve_num{i} = unique([this_seq{i}.target_power]);
    stim_starts{i} = [data.trial_metadata(cur_trial).sequence.start];
    for j = 1:length(this_seq{i})
        if i == 1 && j == 1
            full_seq(1) = this_seq{i}(j);
        else
            full_seq(end+1) = this_seq{i}(j);
        end
        full_seq(end).precomputed_target_index = ...
            full_seq(end).precomputed_target_index + size(full_stim_key,1);
    end
    full_stim_key = [full_stim_key; this_stim_key{i}];
end
power_curve_num = unique([power_curve_num{:}]);
maps = cell(length(power_curve_num),1);
mpp_maps = cell(length(power_curve_num),1);
map_index = cell(length(power_curve_num),1);
corr_maps = cell(length(power_curve_num),1);
stddev_maps = cell(length(power_curve_num),1);
this_seq = [this_seq{:}];
max_trial = length(this_seq);
% max_trial = 1200;
[traces_ch1,traces_ch2] = ...
    get_stim_stack(data,trials,...
        stims_per_trial,stim_starts);
