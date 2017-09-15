function [maps,mpp_maps,map_index,corr_maps,stddev_maps] = ...
    summarize_map(data,trials,show_raw_data,do_stdmap,do_corrmap)



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
stim_inds = [full_seq.precomputed_target_index];
% on_cell_trials = isnan(full_stim_key(stim_inds,1,2));
on_cell_trials = ones(size([full_seq.target_power]))';
for i = 1:length(power_curve_num)
    
    
%     this_seq = this_seq(1:max_trial);
    traces_pow{1} = traces_ch1(on_cell_trials' & [full_seq.target_power] == power_curve_num(i),:);
    traces_pow{2} = traces_ch2(on_cell_trials' & [full_seq.target_power] == power_curve_num(i),:);
    this_seq_power = full_seq(on_cell_trials' & [full_seq.target_power] == power_curve_num(i));
%     mpp_pow = mpp(on_cell_trials' & [full_seq.target_power] == power_curve_num(i));
    [maps{i},mpp_maps{i},map_index{i},corr_maps{i},stddev_maps{i}] = ...
        see_grid_multi(traces_pow,[],this_seq_power,full_stim_key,1,show_raw_data,do_stdmap,do_corrmap);
    if show_raw_data
        title(['Power = ' num2str(power_curve_num(i)) ' mW'])
    end
%     get_mean_events_image(mpp_maps{i}, 2000, 1, 1);
%     title(['Event Counts, Power = ' num2str(power_curve_num(i)) ' mW'])
end

return 

this_seq_plot = acq_gui_data.data.trial_metadata(cur_trial).sequence;
this_stim_key = acq_gui_data.data.trial_metadata(cur_trial).stim_key;
power_curve_num = unique([this_seq_plot.target_power]);
show_raw_data = 0;
for j = 1:length(power_curve_num)
    [traces_ch1,traces_ch2] = ...
    get_stim_stack(acq_gui_data.data,cur_trial,...
        length(this_seq_plot),[this_seq_plot.start]);

    traces_pow{1} = traces_ch1([this_seq_plot.target_power] == power_curve_num(j),:);
    traces_pow{2} = traces_ch2([this_seq_plot.target_power] == power_curve_num(j),:);
    this_seq_power = this_seq_plot([this_seq_plot.target_power] == power_curve_num(j));
    maps,mpp_maps,map_index,corr_maps,stddev_maps
    [maps{j},mpp_maps{j},map_index{j},corr_maps{j},stddev_maps{j}] = ...
        see_grid_multi(traces_pow,this_seq_power,this_stim_key,5,show_raw_data);
    if show_raw_data
        title(['Power = ' num2str(power_curve_num(j)) ' mW'])
    end
end