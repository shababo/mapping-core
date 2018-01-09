function result = analyze_current_diffraction_map(data,experiment_setup,trials)

num_locs = length(trials);


for i = 1:num_locs
    
    [curr_traces, ~,~, stim_traces] = get_traces(data,trials(i));
    result.traces{experiment_setup.pos_order(i)} = curr_traces;
    result.stim_size{experiment_setup.pos_order(i)} = mean(stim_traces(:,.0005*20000:.0025*20000),2);
    result.max_curr{experiment_setup.pos_order(i)} = curr_traces(:,1) - min(curr_traces,[],2);
    
end