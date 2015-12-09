function scratch_the_func()

filenames = {'12_1_slice1_cell1.mat',...
'12_1_slice1_cell2.mat',...
'12_1_slice3_cell1.mat',...
'12_1_slice3_cell2.mat',...
'12_1_slice4_cell2.mat',...
'12_1_slice5_cell1.mat',...
'12_1_slice5_cell2.mat',...
'12_2_slice1_cell1.mat',...
'12_3_slice1_cell2.mat'};

%'12_3_slice2_cell1.mat',...
%'12_3_slice3_cell1.mat',...
%'12_3_slice5_cell1'

run_count_id = [9
2
10
3
13
15
2
6
6
13
7
13];




trial_ids = [-60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0  -60 -45 -30 -15 0 15 30 45 60
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';


peak_currents_cells_by_trial = cell(length(filenames),size(trial_ids,1));
peak_currents_cells_by_trial_mean = zeros(length(filenames),size(trial_ids,1));
peak_currents_cells_by_trial_std = zeros(length(filenames),size(trial_ids,1));
peak_currents_trial_mean = zeros(size(trial_ids,1),1);
peak_currents_trial_std = zeros(size(trial_ids,1),1);

baseline_window = 20000*[.299 .300]; measure_window = 20000*[.301 .330];

for i = 1:length(filenames)
    
    load(['data/' filenames{i}])
    
    [traces, traces_metadata] = get_sweeps_dir('data',filenames{i},0,1,0,Inf,'run_count',run_count_id(i));
    traces = traces{1};
    
    params1.run_count = run_count_id(i);
    match_inds = match_trials(params1, traces_metadata{1})
    traces = traces(match_inds,:);
    temp = traces_metadata{1};
    traces_metadata = temp(match_inds);
    
    trial_types = zeros(size(traces,1),1);
    
    for j = 1:size(trial_ids,1)
        
        params2.relative_position = trial_ids(j,:);
        match_inds = match_trials(params2, traces_metadata);
%         size(match_inds)
        trial_types(match_inds) = j;
        
        peak_currents_cells_by_trial{i,j} = get_current_amp(traces(match_inds,:),baseline_window,measure_window);
        peak_currents_cells_by_trial_mean(i,j) = mean(peak_currents_cells_by_trial{i,j});
        peak_currents_cells_by_trial_std(i,j) = std(peak_currents_cells_by_trial{i,j});
        
    end


end

for j = 1:length(trial_ids)
    
    peak_currents_trial_mean(j) = mean(peak_currents_cells_by_trial_mean(:,j));
    peak_currents_trial_std(j) = mean(peak_currents_cells_by_trial_std(:,j));
    
end