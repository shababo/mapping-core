function [trace_grid, current_image] = get_mapping_data(filename,run,min_or_max,plot_grid,plot_image)%,max_traces,trace_limits, plot_grid,min_or_max,plot_image)

% trace_grid = build_trace_grid(filename, run, 10, max_traces, trace_limits, plot_grid);
% 
% current_image = get_current_image(trace_grid,min_or_max,plot_image);

trace_grid = build_trace_grid(filename, run, 3, [], plot_grid);

current_image = get_current_image(trace_grid,min_or_max,plot_image);