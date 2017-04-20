function [spike_data, voltage_data, current_data, intrinsics] = preprocess_opto_response(data,start_trial,spike_thresh,num_spike_locs,do_cc,do_vc,stim_start,cell_pos)

% power response curves
spike_data = struct();
voltage_data = struct();
current_data = struct();
intrinsics = struct();

for j = 1:num_spike_locs
    
    trial = j+start_trial-1;
    this_seq = data.trial_metadata(trial).sequence;
    powers = unique([this_seq.target_power]);
    [trace_stack] = ...
        get_stim_stack(data,trial,...
        length(this_seq));
    trace_grid = cell(length(powers),1);
    spike_times_grid = cell(length(powers),1);
    for i = 1:length(powers)
        these_traces = trace_stack([this_seq.target_power] == powers(i),:);
        trace_grid{i} = these_traces;
        spike_times_grid{i} = detect_peaks(-bsxfun(@minus,these_traces,median(these_traces(:,60:100),2)),spike_thresh,200,0,1,-Inf,0,0,0);
    end

    spike_data(j).data = trace_grid;
    spike_data(j).spike_times = spike_times_grid;
    spike_data(j).powers = powers;
    stim_ind = this_seq(1).precomputed_target_index;
    stim_pos = data.trial_metadata(trial).stim_key(stim_ind,:) + round(data.trial_metadata(trial).relative_position) - cell_pos;
    spike_data(j).location = stim_pos;

    
    spike_data(j).num_spike_means = cell2float(cellfun(@(x) nansum(sign(cell2float(x))),spike_times_grid,'UniformOutput',0))./cellfun(@(x) length(x),spike_times_grid);
    spike_data(j).spike_times_means = cell2float(cellfun(@(x) nanmean(cell2float(x) - stim_start),spike_times_grid,'UniformOutput',0));
    spike_data(j).spike_times_std = cell2float(cellfun(@(x) nanstd(cell2float(x) - stim_start),spike_times_grid,'UniformOutput',0));
    
    if do_cc
        trial = j+start_trial-1+num_spike_locs+2;
        this_seq = data.trial_metadata(trial).sequence;
        powers = unique([this_seq.target_power]);
        [trace_stack] = ...
            get_stim_stack(data,trial,...
            length(this_seq));
        trace_grid = cell(length(powers),1);
        spike_times_grid = cell(length(powers),1);
        for i = 1:length(powers)
            trace_grid{i} = trace_stack([this_seq.target_power] == powers(i),:);
            spike_times_grid{i} = detect_peaks(trace_grid{i},-30,200,0,1,-Inf,0,0,0);
        end

        voltage_data(j).data = trace_grid;
        voltage_data(j).spike_times = spike_times_grid;
        voltage_data(j).powers = powers;

        stim_ind = this_seq(1).precomputed_target_index;
        stim_pos = data.trial_metadata(trial).stim_key(stim_ind,:) + round(data.trial_metadata(trial).relative_position) - cell_pos;
        voltage_data(j).location = stim_pos;

    end
    
    if j == 1 && do_vc
        if do_cc
            trial = start_trial+num_spike_locs*2+2;
        else
            trial = start_trial+num_spike_locs+2;
        end
        this_seq = data.trial_metadata(trial).sequence;
        powers = unique([this_seq.target_power]);
        [trace_stack] = ...
            get_stim_stack(data,trial,...
            length(this_seq));
        trace_grid = cell(length(powers),1);
        for i = 1:length(powers)
            trace_grid{i} = trace_stack([this_seq.target_power] == powers(i),:);
        end  
        stim_ind = this_seq(1).precomputed_target_index;
        stim_pos = data.trial_metadata(trial).stim_key(stim_ind,:) + round(data.trial_metadata(trial).relative_position) - cell_pos;
        current_data.power_response = trace_grid;
        current_data.powers = powers;
        current_data.power_response_pos = stim_pos;
        
        peak_currents = cellfun(@(y) cellfun(@(x) -min(x(100:200) - mean(x(80:100))),num2cell(y,2)),trace_grid,'UniformOutput',0);
        peak_currents_means = cellfun(@(x) mean(x),peak_currents);
        peak_currents_stds = cellfun(@(x) std(x),peak_currents);
        current_data.power_response_means = peak_currents_means;
        current_data.power_response_stds = peak_currents_stds;
        
    end
end

if ~do_vc
    return
end

intrinsics_trial = start_trial + num_spike_locs + 1;
intrinsics.data = data.sweeps{intrinsics_trial}(:,1)';

if do_cc
    vc_trial_1 = start_trial + num_spike_locs*2 + 3;
else
    vc_trial_1 = start_trial + num_spike_locs + 3;
end

z_depths = [data.trial_metadata.relative_position];
z_depths = unique(round(z_depths((2+start_trial):3:end),-1));
num_depths = length(z_depths);
largest_grid = 9;
all_trials = [];
all_inds = [];
all_targets = [];
shape_svd = zeros(largest_grid,largest_grid,num_depths);
shape_max = zeros(largest_grid,largest_grid,num_depths);

for j = vc_trial_1+(0:num_depths-1)
    
     
    trial = j
    this_seq = data.trial_metadata(trial).sequence;
    stim_key = data.trial_metadata(trial).stim_key([this_seq.precomputed_target_index],:);
    targets = bsxfun(@plus,stim_key,round(data.trial_metadata(trial).relative_position,-1) - cell_pos);
    inds = [targets(:,[1 2])/10 + 5 repmat(find(z_depths == targets(1,3)),size(targets,1),1)]; % DON'T HARD CODE THIS
    [trace_stack] = ...
        get_stim_stack(data,trial,...
        length(this_seq));   
    all_trials = [all_trials; trace_stack];
    all_inds = [all_inds; inds];
    all_targets = [all_targets; targets];
    
end

current_data.shape_data = all_trials;
current_data.shape_locations = all_targets;
current_data.shape_inds = all_inds;

peak_currents = cellfun(@(x) -min(x(100:200) - mean(x(80:100))),num2cell(all_trials,2));
current_data.peak_currents = peak_currents;

[u,s,v] = svd(-all_trials);
if mean(v(:,1)) < 0
    svd_weights = -u(:,1);
else
    svd_weights = u(:,1);
end
for i = 1:length(peak_currents)
    shape_max(all_inds(i,1),all_inds(i,2),all_inds(i,3)) = peak_currents(i);
    shape_svd(all_inds(i,1),all_inds(i,2),all_inds(i,3)) = svd_weights(i);
end
current_data.shape_max = shape_max;
current_data.shape_svd = shape_svd;

all_trials = all_trials(peak_currents < 1200,:);
[u,s,v] = svd(-all_trials);
if mean(v(:,1)) < 0
    current_shape = -v(:,1);
else
    current_shape = v(:,1);
end

current_shape = current_shape - current_shape(1);
current_data.current_shape = current_shape/max(current_shape);
    
    
    
    