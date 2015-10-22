function selected_traces = create_trace_matrix(filename, start_ind, trace_inds)

load(filename,'sweeps')

trace_length = length(sweeps{1}(start_ind:end,1));

selected_traces = zeros(length(trace_inds),trace_length);

if isempty(trace_inds)
    trace_inds = 1:length(sweeps);
end

for i = 1:length(trace_inds)
    
    selected_traces(i,:) = sweeps{trace_inds(i)}(start_ind:end,1);
    
end