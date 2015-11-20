function traces = get_sweeps(filename,ch_ind,trace_inds,plot_stack)

load(filename,'sweeps','data')

start_ind = .3*20000;
end_ind = .530*20000;%.55*20000;

if ~exist('sweeps')
    sweeps = data.sweeps;
end

if ~isempty(trace_inds)
    trace_array = sweeps(trace_inds);
else
    trace_array = sweeps;
end
traces = zeros(length(trace_array),length(trace_array{1}(start_ind:end_ind,ch_ind)'));
for i = 1:length(trace_array), traces(i,:) = trace_array{i}(start_ind:end_ind,ch_ind)'; end

% traces = bsxfun(@minus,traces,mean(traces,1));

if plot_stack
    figure; plot_trace_stack(traces,100,zeros(length(traces),3),'-')
end