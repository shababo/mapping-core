function result = analyze_current_diffraction_map(data,pos_order,trials,measure_type)

num_locs = length(trials);

if all(cellfun(@isempty,trials))
    result.traces{1} = [];
    result.stim_size{1} = [];
    switch measure_type
        case 'current'
            result.max_curr{1} = [];
        case 'spikes'
        result.spike_times{1} = [];
    end
else
    for i = 1:num_locs

        [traces, ~,~, stim_traces] = get_traces(data,trials{i});
        result.traces{pos_order(i)} = traces;
        result.stim_size{pos_order(i)} = mean(stim_traces(:,.0005*20000:.0025*20000),2);
        switch measure_type
            case 'current'
                result.max_curr{pos_order(i)} = traces(:,1) - min(traces,[],2);
            case 'spikes'
                result.spike_times{pos_order(i)} = detect_peaks(-bsxfun(@minus,traces,traces(:,1)),20,30,0,Inf,-Inf,0,0,1);
        end

    end
end