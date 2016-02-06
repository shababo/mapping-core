function traces_by_location = build_trace_grid(filename, runs, step_size, max_traces, trace_limits, plot_grid)


[traces, traces_metadata] =  get_sweeps(filename,1,[],0,Inf,'run_count',runs);
% get_sweeps_dir(directory,file_string,0,1,0,Inf,'run_count',runs);

size(traces)

% traces = traces{1};
% traces_metadata = traces_metadata{1};

traces_by_location = cell(11,11);
if ~isempty(trace_limits)
    start_ind = trace_limits(1); end_ind = trace_limits(2);
else
    start_ind = 1; end_ind = size(traces,2);
end

% find limits of grid
x_min = Inf; x_max = -Inf; y_min = Inf; y_max = -Inf;

for i = 1:length(traces_metadata)
    this_location = traces_metadata(i).relative_position;
    if x_min > this_location(1)
        x_min = this_location(1);
    end
    if x_max < this_location(1)
        x_max = this_location(1);
    end
    if y_min > this_location(2)
        y_min = this_location(2);
    end
    if y_max < this_location(2)
        y_max = this_location(2);
    end
end

num_x_positions = (x_max - x_min)/step_size + 1
num_y_positions = (y_max - y_min)/step_size + 1


traces_by_location = cell(num_x_positions,num_y_positions);
for i = 1:size(traces,1)
    ind1 = traces_metadata(i).relative_to_start_position(1)/step_size + ceil(num_x_positions/2)
    ind2 = traces_metadata(i).relative_to_start_position(2)/step_size + ceil(num_y_positions/2)
    traces_by_location{end-ind1+1,ind2} = [traces_by_location{end-ind1+1,ind2}; traces(i,start_ind:end_ind)];
end

if plot_grid
    figure; plot_trace_stack_grid(traces_by_location,max_traces,3,0)
end