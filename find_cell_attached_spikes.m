function spike_traces = find_cell_attached_spikes(traces, spike_window, threshold, refract_period)

size(traces)
spike_traces = zeros(size(traces));
for i = 1:size(traces,1)
    trace = traces(i,:);
    size(trace)
    spike_trace = [0 trace(1:end-2) > trace(2:end-1) & trace(2:end-1) < trace(3:end) & trace(2:end-1) < median(trace(spike_window(1):spike_window(2)))-threshold 0];
    spike_trace(1:spike_window(1)) = 0;
    spike_trace(spike_window(2):end) = 0;
    spike_times = find(spike_trace == 1);
    
    isi = diff(spike_times);
    while any(diff(spike_times) < refract_period)
        for j = 1:length(isi)
            if isi(j) < refract_period
                spike_times(j+1) = [];
                isi = diff(spike_times);
                break
            end
        end
    end
    spike_traces(i,spike_times) = 1;
end