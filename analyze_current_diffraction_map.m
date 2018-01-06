function result = analyze_current_diffraction_map(data,experiment_setup,trials)

num_locs = length(trials);
stim_starts = 20*[data.sequence.start];

for i = 1:num_locs
    result.traces{experiment_setup.pos_order(i)} = get_traces(data,trials(i),