function results = analyze_spatial_resolution_spiking(filenames,run_count_id,trial_ids,threshold, do_plot)

num_spikes_cells_by_trial = cell(length(filenames),size(trial_ids,1));
num_spikes_cells_by_trial_mean = zeros(length(filenames),size(trial_ids,1));
num_spikes_cells_by_trial_std = zeros(length(filenames),size(trial_ids,1));
num_spikes_trial_mean = zeros(size(trial_ids,1),1);
num_spikes_trial_std = zeros(size(trial_ids,1),1);

spike_window = 20000*[.300 .450];

for i = 1:length(filenames)
    
    load(['C:\data\Shababo\' filenames{i}])
    
    [traces, traces_metadata] = get_sweeps_dir('C:\data\Shababo',filenames{i},0,1,0,Inf,'run_count',run_count_id{i});
    traces = traces{1};
    size(traces)
    assignin('base','traces_metadata',traces_metadata)
    
    params1.run_count = run_count_id{i};
    match_inds = match_trials(params1, traces_metadata{1});
    traces = traces(match_inds,:);
    temp = traces_metadata{1};
    traces_metadata = temp(match_inds);
    assignin('base','traces_metadata',traces_metadata)
    trial_types = zeros(size(traces,1),1);
    

    for j = 1:size(trial_ids,1)
        
        params2.relative_position = trial_ids(j,:);
        %params2.hologram_id = 'ROI3';
%         params3.relative_position = trial_ids2(j,:);
%         match_inds = unique([match_trials(params2, traces_metadata) match_trials(params3, traces_metadata)]);
        match_inds = unique([match_trials(params2, traces_metadata)]);
        size(match_inds)
        if isempty(match_inds)
            ['C:\data\Shababo\' filenames{i}]
%             assignin('base','whatthe',traces_metadata)
        end
        trial_types(match_inds) = j;
                
        spike_traces = find_cell_attached_spikes(traces(match_inds,:), spike_window, threshold, 20);
        size(traces)
        sum(spike_traces,2)
        if any(sum(spike_traces,2) == 2)
            figure
            j
            match_inds
            doubles = find(sum(spike_traces,2) == 2)
            plot(spike_traces(doubles(1),:))
        end
        assignin('base','spike_traces',spike_traces)
        num_spikes_cells_by_trial{i,j} = sum(spike_traces,2);
        num_spikes_cells_by_trial_mean(i,j) = mean(num_spikes_cells_by_trial{i,j});
        num_spikes_cells_by_trial_std(i,j) = std(num_spikes_cells_by_trial{i,j});
        
    end


end


for j = 1:length(trial_ids)
    
    num_spikes_trial_mean(j) = mean(num_spikes_cells_by_trial_mean(:,j));
    num_spikes_trial_std(j) = mean(num_spikes_cells_by_trial_std(:,j));
    
end

assignin('base','num_spikes_cells_by_trial_mean',num_spikes_cells_by_trial_mean)
assignin('base','num_spikes_cells_by_trial',num_spikes_cells_by_trial)
results = [];

if do_plot
    
    switch_ind = length(num_spikes_trial_std)/3;

    positions =trial_ids(1:switch_ind,1);

    
    colors = lines(length(filenames));
    figure; %h = plot(positions,num_spikes_cells_by_trial_mean(:,1:switch_ind)');
    hold on;
    for i = 1:length(filenames)
        for j = 1:switch_ind
            plot(positions,num_spikes_cells_by_trial_mean(i,1:switch_ind),'-.','color',colors(i,:))
            %hold on
            %scatter(positions(j)*ones(length(num_spikes_cells_by_trial{i,j}),1),num_spikes_cells_by_trial{i,j},[],repmat(colors(i,:),length(num_spikes_cells_by_trial{i,j}),1),'filled');
            %hold on;
        end

    end
    % legend(filenames)
    title('x')
    figure; %h = plot(positions,num_spikes_cells_by_trial_mean(:,switch_ind+1:end)');
    hold on;
    for i = 1:length(filenames)
        for j = switch_ind+1:switch_ind*2
            plot(positions,num_spikes_cells_by_trial_mean(i,switch_ind+1:switch_ind*2),'-.','color',colors(i,:))
            hold on
            %scatter(positions(j-switch_ind)*ones(length(num_spikes_cells_by_trial{i,j}),1),num_spikes_cells_by_trial{i,j},[],repmat(colors(i,:),length(num_spikes_cells_by_trial{i,j}),1),'filled');
            %hold on;
        end
    end
    % legend(filenames)
    title('y')

    figure; %h = plot(positions,num_spikes_cells_by_trial_mean(:,switch_ind+1:end)');
    hold on;
    for i = 1:length(filenames)
        for j = switch_ind*2+1:size(trial_ids,1)
            plot(positions,num_spikes_cells_by_trial_mean(i,switch_ind*2+1:end),'-.','color',colors(i,:))
            %hold on
            %scatter(positions(j-switch_ind*2)*ones(length(num_spikes_cells_by_trial{i,j}),1),num_spikes_cells_by_trial{i,j},[],repmat(colors(i,:),length(num_spikes_cells_by_trial{i,j}),1),'filled');
            %hold on;
        end
    end
    % legend(filenames)
    title('z')
end