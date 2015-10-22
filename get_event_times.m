function event_times = get_event_times(results,traces,traces_inds,start,stop)

event_times = cell(length(traces_inds),1);

for i = 1:length(traces_inds)
    
    trace_ind = traces_inds(i);

    trials = results(trace_ind).trials;
    times = results(trace_ind).times;
    mcmc = results(trace_ind).mcmc;
    trace = max(traces(trace_ind,:)) - traces(trace_ind,:);

    errP = zeros(1,length(trials.curves));
    for j = 1:length(trials.curves)
        errP(j) = sum((trials.curves{j}-trace).^2);
    end
    % figure;plot(errP)

    [me mi] = min(errP);
    
    times = trials.times{mi};
    times(times < start) = [];
    times(times > stop) = [];
    
    event_times{i} = times;
end
    
    
