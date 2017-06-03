function [spike_data, voltage_data, current_data, intrinsics] = ...
    preprocess_opto_response(data,start_trial,cc_spike_thresh,...
    ca_num_spike_locs,ca_spike_thresh,do_cc,cc_num_spike_locs,do_vc,...
    stim_start,cell_pos,first_spike_only,trial_dur,do_hpf)

% power response curves
spike_data = struct();
voltage_data = struct();
current_data = struct();
intrinsics = struct();

for j = 1:ca_num_spike_locs
    
    trial = j+start_trial-1;
    this_seq = data.trial_metadata(trial).sequence;
    powers = unique([this_seq.target_power]);
    expected_stim_starts = {[this_seq.start]};
    [trace_stack] = ...
        get_stim_stack(data,trial,...
        length(this_seq),expected_stim_starts);
    trace_grid = cell(length(powers),1);
    spike_times_grid = cell(length(powers),1);
    spike_data(j).num_spike_means = zeros(1,length(powers));
    spike_data(j).spike_times_means = zeros(1,length(powers));
    spike_data(j).spike_times_std = zeros(1,length(powers));
    for i = 1:length(powers)
        these_traces = trace_stack([this_seq.target_power] == powers(i),:);
        trace_grid{i} = these_traces;
        spike_times = [];
        spike_times_first = [];
        spike_times_grid{i} = detect_peaks(-bsxfun(@minus,these_traces,median(these_traces(:,1:10),2)),ca_spike_thresh,30,0,trial_dur,-Inf,0,do_hpf,first_spike_only);
        for  k = 1:length(spike_times_grid{i})
            if ~isempty(spike_times_grid{i}{k})
                spike_times = [spike_times spike_times_grid{i}{k}];
                spike_times_first = [spike_times_first spike_times_grid{i}{k}(1)];
            end
        end
        spike_data(j).num_spike_means(i) = length(spike_times)/length(spike_times_grid{i});
        spike_data(j).spike_times_means(i) = mean(spike_times_first);
        spike_data(j).spike_times_std(i) = std(spike_times_first);
    end

    spike_data(j).data = trace_grid;
    spike_data(j).spike_times = spike_times_grid;
    spike_data(j).powers = powers;
    stim_ind = this_seq(1).precomputed_target_index;
    stim_pos = data.trial_metadata(trial).stim_key(stim_ind,:) + round(data.trial_metadata(trial).relative_position) - cell_pos;
    spike_data(j).location = stim_pos;

    
%     spike_data(j).num_spike_means = cell2float(cellfun(@(x) sum(cellfun(@(y) ~isempty(y),x)),spike_times_grid,'UniformOutput',0))./cellfun(@(x) length(x),spike_times_grid);
%     spike_data(j).spike_times_means = cell2float(cellfun(@(x) nanmean(cell2float(x)),spike_times_grid,'UniformOutput',0));
%     spike_data(j).spike_times_std = cell2float(cellfun(@(x) nanstd(cell2float(x)),spike_times_grid,'UniformOutput',0));
    
    
end

if do_cc
    for j = 1:cc_num_spike_locs
        
        

        trial = j+start_trial-1+ca_num_spike_locs+2;
        this_seq = data.trial_metadata(trial).sequence;
        powers = unique([this_seq.target_power]);
        expected_stim_starts = {[this_seq.start]};
        [trace_stack] = ...
            get_stim_stack(data,trial,...
            length(this_seq),expected_stim_starts);
        trace_grid = cell(length(powers),1);
        spike_times_grid = cell(length(powers),1);
        voltage_data(j).num_spike_means = zeros(1,length(powers));
        voltage_data(j).spike_times_means = zeros(1,length(powers));
        voltage_data(j).spike_times_std = zeros(1,length(powers));
        
        for i = 1:length(powers)
            trace_grid{i} = trace_stack([this_seq.target_power] == powers(i),:);
            spike_times_grid{i} = detect_peaks(trace_grid{i},cc_spike_thresh,30,0,trial_dur,-Inf,0,0,first_spike_only);
            spike_times = [];
            spike_times_first = [];
            for  k = 1:length(spike_times_grid{i})
                if ~isempty(spike_times_grid{i}{k})
                    these_spikes = spike_times_grid{i}{k};
%                     these_spikes = these_spikes(these_spikes < 400);
                    if ~isempty(these_spikes)
                        spike_times = [spike_times these_spikes];
                        spike_times_first = [spike_times_first these_spikes(1)];
                    end
                end
            end
            voltage_data(j).num_spike_means(i) = length(spike_times)/length(spike_times_grid{i});
            voltage_data(j).spike_times_means(i) = mean(spike_times_first);
            voltage_data(j).spike_times_std(i) = std(spike_times_first);
        end

        voltage_data(j).data = trace_grid;
        voltage_data(j).spike_times = spike_times_grid;
        voltage_data(j).powers = powers;

        stim_ind = this_seq(1).precomputed_target_index;
        stim_pos = data.trial_metadata(trial).stim_key(stim_ind,:) + round(data.trial_metadata(trial).relative_position) - cell_pos;
        voltage_data(j).location = stim_pos;
        
%         voltage_data(j).num_spike_means = cell2float(cellfun(@(x) sum(cellfun(@(y) ~isempty(y),x)),spike_times_grid,'UniformOutput',0))./cellfun(@(x) length(x),spike_times_grid);
%         voltage_data(j).spike_times_means = cell2float(cellfun(@(x) nanmean(cell2float(x)),spike_times_grid,'UniformOutput',0));
%         voltage_data(j).spike_times_std = cell2float(cellfun(@(x) nanstd(cell2float(x)),spike_times_grid,'UniformOutput',0));
    end
end

if do_vc
        
    trial = start_trial+ca_num_spike_locs+cc_num_spike_locs+2;

    this_seq = data.trial_metadata(trial).sequence;
    powers = unique([this_seq.target_power]);
    [trace_stack] = ...
        get_stim_stack(data,trial,...
        length(this_seq),{[this_seq.start]});
    trace_grid = cell(length(powers),1);
    for i = 1:length(powers)
        trace_grid{i} = trace_stack([this_seq.target_power] == powers(i),:);
    end  
    stim_ind = this_seq(1).precomputed_target_index;
    stim_pos = data.trial_metadata(trial).stim_key(stim_ind,:) + round(data.trial_metadata(trial).relative_position) - cell_pos;
    current_data.power_response = trace_grid;
    current_data.powers = powers;
    current_data.power_response_pos = stim_pos;

    peak_currents = cellfun(@(y) cellfun(@(x) -min(x(2:200)) - mean(x(1:5)),num2cell(y,2)),trace_grid,'UniformOutput',0);
    peak_currents_means = cellfun(@(x) mean(x),peak_currents);
    peak_currents_stds = cellfun(@(x) std(x),peak_currents);
    current_data.power_response_means = peak_currents_means;
    current_data.power_response_stds = peak_currents_stds;

end

if do_vc || do_cc
    intrinsics_trial = start_trial + ca_num_spike_locs + 1;
    intrinsics.data = data.sweeps{intrinsics_trial}(:,1)';
end

if do_vc

    vc_trial_1 = start_trial + ca_num_spike_locs+cc_num_spike_locs + 3;

    z_depths = [data.trial_metadata(vc_trial_1:end).relative_position];
    z_depths = unique(round(z_depths(3:3:end),-1))
    num_depths = length(z_depths);
    largest_grid = 9;
    all_trials = [];
    all_inds = [];
    all_targets = [];
    shape_svd = zeros(largest_grid,largest_grid,num_depths);
    shape_max = zeros(largest_grid,largest_grid,num_depths);

    for j = vc_trial_1+(0:num_depths-1)


        trial = j;
        this_seq = data.trial_metadata(trial).sequence;
        stim_key = data.trial_metadata(trial).stim_key([this_seq.precomputed_target_index],:);
        targets = bsxfun(@plus,stim_key,round(data.trial_metadata(trial).relative_position,-1) - cell_pos);
        inds = [targets(:,[1 2])/10 + 5 repmat(find(z_depths == targets(1,3)),size(targets,1),1)]; % DON'T HARD CODE THESE VALUES
        [trace_stack] = ...
            get_stim_stack(data,trial,...
            length(this_seq),{[this_seq.start]});   
        all_trials = [all_trials; trace_stack];
        all_inds = [all_inds; inds];
        all_targets = [all_targets; targets];

    end

    current_data.shape_data = all_trials;
    current_data.shape_locations = all_targets;
    current_data.shape_inds = all_inds;

    peak_currents = cellfun(@(x) -min(x(5:50)) + min(x(1:5)),num2cell(all_trials,2));
    current_data.peak_currents = peak_currents;

    [u,s,v] = svd(-all_trials(:,1:50));
    if mean(v(:,1)) < 0
        svd_weights = -u(:,1);
    else
        svd_weights = u(:,1);
    end
    for i = 1:length(peak_currents)
        shape_max(all_inds(i,1),all_inds(i,2),all_inds(i,3)) = peak_currents(i);
        shape_svd(all_inds(i,1),all_inds(i,2),all_inds(i,3)) = svd_weights(i);
    end
    
    % TO DO CLEAN UP SMALL VALUES AND ZEROS AT NO MEASUREMENTS
    
    current_data.shape_max = shape_max;%/max(shape_max(:));
    current_data.shape_svd = shape_svd/max(shape_svd(:));

    all_trials = all_trials(peak_currents < 1200,:);
    [u,s,v] = svd(-all_trials);
    if mean(v(:,1)) < 0
        current_shape = -v(:,1);
    else
        current_shape = v(:,1);
    end

    current_shape = current_shape - current_shape(1);
    current_data.current_shape = current_shape/max(current_shape);
    current_data.current_shape_sv = s(1)/sum(diag(s));

    % interpolate current shape

    minx = min(current_data.shape_locations(:,1));
    maxx = max(current_data.shape_locations(:,1));
    miny = min(current_data.shape_locations(:,2));
    maxy = max(current_data.shape_locations(:,2));
    minz = min(current_data.shape_locations(:,3));
    maxz = max(current_data.shape_locations(:,3));
    current_data.upres = 1;
    current_data.upres_x_locs = minx:current_data.upres:maxx;
    current_data.upres_y_locs = miny:current_data.upres:maxy;
    current_data.upres_z_locs = minz:current_data.upres:maxz;

    start_x = sort(unique(current_data.shape_locations(:,1)));
    start_y = sort(unique(current_data.shape_locations(:,2)));
    start_z = sort(unique(current_data.shape_locations(:,3)));

    
    [X, Y, Z] = meshgrid(start_x,start_y,start_z);
    [X_upres, Y_upres, Z_upres] = ...
        meshgrid(current_data.upres_x_locs, current_data.upres_y_locs, current_data.upres_z_locs);
    current_data.upres_shape = ...
        interp3(X,Y,Z,current_data.shape_max,X_upres,Y_upres,Z_upres,'linear');
%     current_data.upres_shape = current_data.upres_shape/max(current_data.upres_shape(:));
    
%     test_points = -20:10:20;
%     x_inds = ceil((test_points - current_data.upres_x_locs(1))/current_data.upres);
%     y_inds = ceil((test_points - current_data.upres_y_locs(1))/current_data.upres);
%     z_inds = ceil((test_points - current_data.upres_z_locs(1))/current_data.upres);
%     current_data.test_shape = current_data.upres_shape(x_inds,y_inds,z_inds);
%     size(current_data.test_shape)

end

















    
    