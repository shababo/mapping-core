function [trace_grid, current_image] = get_mapping_data(filename,run,min_or_max,plot_grid,plot_image)%,max_traces,trace_limits, plot_grid,min_or_max,plot_image)

% trace_grid = build_trace_grid(filename, run, 10, max_traces, trace_limits, plot_grid);
% 
% current_image = get_current_image(trace_grid,min_or_max,plot_image);
figure
subplot(221)
trace_grid = build_trace_grid(filename, run, 3, [], plot_grid);
subplot(222)
current_image = get_current_image(trace_grid,min_or_max,plot_image);
subplot(223)
corr_image = get_corr_image(trace_grid,plot_image);
subplot(224)
imagesc(current_image .* (corr_image > 3.0*std(corr_image(:))))
