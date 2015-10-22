function spike_trace = find_cell_attached_spikes(trace, start_ind)

spike_trace = [0 trace(1:end-2) < trace(2:end-1) & trace(2:end-1) < trace(3:end) & trace(2:end-1) < mean(trace)-50 0];

spike_trace(1:start_ind) = 0;