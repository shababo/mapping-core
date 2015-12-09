function current_amps = get_current_amp(traces,baseline_window,measure_window)
% size(min(traces(:,baseline_window(1):baseline_window(2))))
current_amps = min(traces(:,baseline_window(1):baseline_window(2))') - min(traces(:,measure_window(1):measure_window(2))');

% size(current_amps)