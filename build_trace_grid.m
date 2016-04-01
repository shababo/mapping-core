function traces_by_location = build_trace_grid(filename, runs, max_traces, trace_limits, plot_grid,varargin)


[traces, traces_metadata] =  get_sweeps(filename,1,[],0,Inf,'run_count',runs,varargin{:});
% get_sweeps_dir(directory,file_string,0,1,0,Inf,'run_count',runs);

size(traces)

% traces = traces{1};
% traces_metadata = traces_metadata{1};600


if ~isempty(trace_limits)
    start_ind = trace_limits(1); end_ind = trace_limits(2);
else
    start_ind = 1; end_ind = size(traces,2);
end

% find limits of grid
x_min = Inf; x_max = -Inf; y_min = Inf; y_max = -Inf;

x_locations = [];
y_locations = [];

for i = 1:length(traces_metadata)
      
    this_location = traces_metadata(i).relative_to_start_position;
    
    x_locations = [x_locations this_location(1)];
    y_locations = [y_locations this_location(2)];
%     if x_min > this_location(1)
%         x_min = this_location(1);
%     end
%     if x_max < this_location(1)
%         x_max = this_location(1);
%     end
%     if y_min > this_location(2)
%         y_min = this_location(2);
%     end
%     if y_max < this_location(2)
%         y_max = this_location(2);
%     end
end

x_locations = sort(unique(x_locations))
y_locations = sort(unique(y_locations))


num_x_positions = length(x_locations);%(x_max - x_min)/step_size + 1;
num_y_positions = length(y_locations);%(y_max - y_min)/step_size + 1;


traces_by_location = cell(num_x_positions,num_y_positions);
for i = 1:size(traces,1)

%     ind1 = traces_metadata(i).relative_to_start_position(1)/step_size + ceil(num_x_positions/2);
%     ind2 = traces_metadata(i).relative_to_start_position(2)/step_size + ceil(num_y_positions/2);
    ind1 = find(x_locations == traces_metadata(i).relative_to_start_position(1));
    ind2 = find(y_locations == traces_metadata(i).relative_to_start_position(2));
    traces_by_location{end-ind1+1,ind2} = [traces_by_location{end-ind1+1,ind2}; traces(i,start_ind:end_ind)];
    
end

if plot_grid
    figure; plot_trace_stack_grid(traces_by_location,max_traces,5,0)
end