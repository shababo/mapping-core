function result = analyze_current_diffraction_map_rev(data,pos_order,start_trial)

num_trials = length(pos_order);
num_locs = length(unique(pos_order));

result.traces = cell(num_locs,1);
[traces, ~,~, stim_traces] = get_traces(data,start_trial:start_trial + num_trials -1);

for i = 1:num_locs
	
    i
    this_loc_trials = find(pos_order == i);
    result.traces{i} = traces(this_loc_trials,:);
    result.stim_size{i} = mean(stim_traces(this_loc_trials,.0005*20000:.0025*20000),2);
    result.max_curr{i} = traces(this_loc_trials,1) - min(traces(this_loc_trials,:),[],2);
    
end