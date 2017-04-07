function [spike_data, voltage_data, current_data] = view_resp_protocol(data,start_trial,spike_thresh,num_spike_locs,stim_start)

% power response curves
spike_data = struct();
voltage_data = struct();
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
        spike_times_grid{i} = detect_peaks(-bsxfun(@minus,these_traces,median(these_traces(:,60:100),2)),spike_thresh,200,0,1,10,0,0,1);
    end
    figure
    if j == 1
        subplot(221)
    else
        subplot(121)
    end
    spike_data(j).data = trace_grid;
    spike_data(j).spike_times = spike_times_grid;
    spike_data(j).powers = [10 25 50 100 150];
    stim_ind = this_seq(1).precomputed_target_index;
    stim_pos = data.trial_metadata(trial).stim_key(stim_ind,:) + round(data.trial_metadata(trial).relative_position) - [50 50 0];
    spike_data(j).location = stim_pos;
    plot_trace_stack_grid(trace_grid,Inf,1,0,[],[],[],spike_times_grid);
    title(['Cell attached: stim location ' mat2str(stim_pos)])
    if j == 1
        subplot(222)
    else
        subplot(122)
    end
    
%     peak_currents = cellfun(@(y) cellfun(@(x) -min(x(100:200) - mean(x(80:100))),num2cell(y,2)),trace_grid,'UniformOutput',0);
    num_spike_means = cell2float(cellfun(@(x) nansum(sign(cell2float(x))),spike_times_grid,'UniformOutput',0))./cellfun(@(x) length(x),spike_times_grid);
    spike_times_means = cell2float(cellfun(@(x) nanmean(cell2float(x) - stim_start),spike_times_grid,'UniformOutput',0));
%     peak_currents_stds = cellfun(@(x) std(x),peak_currents);
%     yyaxis left
%     plot(powers,num_spike_means)
%     yyaxis right
%     plot(powers,spike_times_means)
    ax = plotyy(powers,num_spike_means,powers,spike_times_means/20)
    title('spike stats')
    xlabel('power (mW)')
    ylabel(ax(1),'Prob. Spike')
    ylabel(ax(2),'Mean 1st Spike Time (ms)')
    
%     trial = j+start_trial-1+num_spike_locs+2;
%     this_seq = data.trial_metadata(trial).sequence;
%     powers = unique([this_seq.target_power]);
%     [trace_stack] = ...
%         get_stim_stack(data,trial,...
%         length(this_seq));
%     trace_grid = cell(length(powers),1);
%     for i = 1:length(powers)
%         trace_grid{i} = trace_stack([this_seq.target_power] == powers(i),:);
%     end
%     if j == 1
%         subplot(132)
%     else
%         subplot(122)
%     end
%     voltage_data(j).data = trace_grid;
%     voltage_data(j).powers = [10 25 50 100 150];
%     plot_trace_stack_grid(flipud(trace_grid),Inf,1,0);
%     stim_ind = this_seq(1).precomputed_target_index;
%     stim_pos = data.trial_metadata(trial).stim_key(stim_ind,:) + round(data.trial_metadata(trial).relative_position) - [50 50 0];
%     voltage_data(j).location = stim_pos;
%     title(mat2str(stim_pos))
%     title(['Current clamp'])
    
    if j == 1
        trial = j+start_trial-1+num_spike_locs+2;
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
        stim_pos = data.trial_metadata(trial).stim_key(stim_ind,:) + round(data.trial_metadata(trial).relative_position) - [50 50 0];
        current_data.power_response = trace_grid;
        current_data.powers = powers;
        subplot(223)
        plot_trace_stack_grid(trace_grid,Inf,1,0);
        title(mat2str(stim_pos))
        title(['Voltage clamp'])
        
        peak_currents = cellfun(@(y) cellfun(@(x) -min(x(100:200) - mean(x(80:100))),num2cell(y,2)),trace_grid,'UniformOutput',0);
        peak_currents_means = cellfun(@(x) mean(x),peak_currents);
        peak_currents_stds = cellfun(@(x) std(x),peak_currents);
        subplot(224)
        errorbar(powers,peak_currents_means,peak_currents_stds)
        hold on
        
        title('Current vs. Power Response')
        xlabel('Power (mW)')
        ylabel('Peak Current (pA)')
    end
end

vc_trial_1 = start_trial + num_spike_locs + 3;

% trial = vc_trial_1;
% this_seq = data.trial_metadata(trial).sequence;
% powers = unique([this_seq.target_power]);
% [trace_stack] = ...
%     get_stim_stack(data,trial,...
%     length(this_seq));
% current_power_grid = cell(length(powers),1);
% for i = 1:length(powers)
%     these_traces = trace_stack([this_seq.target_power] == powers(i),:);
%     current_power_grid{i} = these_traces;
% end


% z_depths = [-60 -40 -20 -10 0 10 20 40 60];
z_depths = [-90 -50 -20 0 20 50 90];
num_depths = length(z_depths);
largest_grid = 9;
all_trials = [];
all_inds = [];
all_targets = [];
shape_svd = zeros(largest_grid,largest_grid,num_depths);
shape_max = zeros(largest_grid,largest_grid,num_depths);
count = 1;
for j = vc_trial_1+(0:num_depths-1)
    
     
    trial = j;
    this_seq = data.trial_metadata(trial).sequence;
    stim_key = data.trial_metadata(trial).stim_key([this_seq.precomputed_target_index],:);
    targets = bsxfun(@plus,stim_key,round(data.trial_metadata(trial).relative_position,-1) - [50 50 0]);
    inds = [targets(:,[1 2])/10 + 5 repmat(find(z_depths == targets(1,3)),size(targets,1),1)];
    [trace_stack] = ...
        get_stim_stack(data,trial,...
        length(this_seq));   
    all_trials = [all_trials; trace_stack];
    all_inds = [all_inds; inds];
    all_targets = [all_targets; targets];
end
current_data.data = all_trials;
current_data.locations = all_targets;
current_data.inds = all_inds;
assignin('base','all_trials',all_trials)
assignin('base','all_inds',all_inds)
peak_currents = cellfun(@(x) -min(x(100:200) - mean(x(80:100))),num2cell(all_trials,2));
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
if any(shape_max(:) > 2000)
    shape_max(shape_max > 2000) = -Inf;
    shape_svd(shape_max > 2000) = -Inf;
    new_max1 = max(shape_max(:));
    new_max2 = max(shape_svd(:));
    shape_max(isinf(shape_max)) = new_max1;
    shape_svd(isinf(shape_svd)) = new_max2;
end
% figure; implay(shape_max)
assignin('base','shape_max',shape_max)
% grid_sizes = [1 1 2 3 4 3 2 1 1];
grid_sizes = [3 3 4 4 4 3 3];
grid_strs = {'3x3','5x5','7x7','9x9'};
figure
max_cur = max(shape_max(:));
for i = 1:num_depths
    subplot(num_depths,1,i)
    imagesc(shape_max(:,:,i))
    caxis([0 max_cur])
    title(sprintf(['z = ' num2str(z_depths(i))])) %grid size = ' grid_strs{grid_sizes(i)} '\n
    axis image
    axis off
end
set(gcf,'Position',[680 1 135 973])
% figure
% max_cur = max(shape_svd(:));
% for i = 1:num_depths
%     subplot(num_depths,1,i)
%     imagesc(shape_svd(:,:,i))
%     caxis([0 max_cur])
%     title(sprintf(['z = ' num2str(z_depths(i))])) %grid size = ' grid_strs{grid_sizes(i)} '\n
%     axis image
%     axis off
% end
    
figure;
% subplot(121)
plot((0:2000)/20,v(:,1))
title('SVD Current Template')
xlim([0 100])
xlabel('time (ms)')
ylabel('Current Amplitude (a.u.)')
% subplot(122)
% plot(u(:,1))
% title('SVD Spatial Weights')  
    
    
    
    