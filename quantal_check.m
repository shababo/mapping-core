num_traces = length(results);

tau1 = [];
tau2 = [];
amps = [];


for i = 1:num_traces
    
    sample_ind = 100;
    
    num_events = length(results(i).trials.tau{sample_ind});
    
    for j = 1:num_events
        tau1 = [tau1 results(i).trials.tau{sample_ind}{j}(1)];
        tau2 = [tau2 results(i).trials.tau{sample_ind}{j}(2)];
        amps = [amps results(i).trials.amp{sample_ind}(j)];
    end
end