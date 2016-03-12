function results = analyze_spatial_resolution(filenames,run_count_id,trial_ids)

peak_currents_cells_by_trial = cell(length(filenames),size(trial_ids,1));
peak_currents_cells_by_trial_mean = zeros(length(filenames),size(trial_ids,1));
peak_currents_cells_by_trial_std = zeros(length(filenames),size(trial_ids,1));
peak_currents_trial_mean = zeros(size(trial_ids,1),1);
peak_currents_trial_std = zeros(size(trial_ids,1),1);

baseline_window = 20000*[.299 .300]; measure_window = 20000*[.301 .330];

for i = 1:length(filenames)
    
    load([filenames{i}])
    
    [traces, traces_metadata] = get_sweeps_dir('',filenames{i},0,1,0,Inf,'run_count',run_count_id{i});
    traces = traces{1};
    
    params1.run_count = run_count_id{i};
    match_inds = match_trials(params1, traces_metadata{1});
    traces = traces(match_inds,:);
    temp = traces_metadata{1};
    traces_metadata = temp(match_inds);
    
    trial_types = zeros(size(traces,1),1);
    

    for j = 1:size(trial_ids,1)
        
        params2.relative_position = trial_ids(j,:);
%         params3.relative_position = trial_ids2(j,:);
%         match_inds = unique([match_trials(params2, traces_metadata) match_trials(params3, traces_metadata)]);
        match_inds = unique([match_trials(params2, traces_metadata)]);
        size(match_inds)
        if isempty(match_inds)
            ['data/' filenames{i}]
%             assignin('base','whatthe',traces_metadata)
        end
        trial_types(match_inds) = j;
        
        upper_limit = 500;
        
        peak_currents_cells_by_trial{i,j} = get_current_amp(traces(match_inds,:),baseline_window,measure_window);
        peak_currents_cells_by_trial{i,j}(peak_currents_cells_by_trial{i,j} > upper_limit) = [];
        peak_currents_cells_by_trial_mean(i,j) = mean(peak_currents_cells_by_trial{i,j});
        peak_currents_cells_by_trial_std(i,j) = std(peak_currents_cells_by_trial{i,j});
        
    end


end


for j = 1:length(trial_ids)
    
    peak_currents_trial_mean(j) = mean(peak_currents_cells_by_trial_mean(:,j));
    peak_currents_trial_std(j) = mean(peak_currents_cells_by_trial_std(:,j));
    
end