function compare_trace_stack_grid_overlap(traces_arrays,in_max_traces,downsample_rate,in_axes,plot_avg, varargin)

if ~isempty(varargin)
    plot_names = varargin{1};
end

if ~isempty(varargin) && length(varargin) > 1
    num_plot_rows = varargin{2};
else
    num_plot_rows = 1;
end

num_arrays = length(traces_arrays);
all_axes = [];

if isempty(in_axes)
    ax = gca;

else
    ax = in_axes;

end

axes(ax)

colors = [1 0 0; 0 0 1];

for array_i = 1:num_arrays
    

    
    this_array = traces_arrays{array_i};
    
    [num_rows, num_cols] = size(this_array);

    count = 1;
    grid_offset_y_spacer = 150;
    if array_i == 1
        grid_offset_y = -grid_offset_y_spacer;
    end
    grid_offset_x = .01*20000;

    if isinf(in_max_traces)
        max_traces = 0;
        for i = 1:num_cols
            for j = 1:num_rows
                num_traces = size(this_array{j,i},2);
                if num_traces > max_traces
                    max_traces = num_traces;
                end
            end
        end
    else
        max_traces = in_max_traces;
    end


    for i = 1:num_cols



        for j = 1:num_rows

            if ~isempty(this_array{j,i})
                these_traces = this_array{j,i};
                if size(these_traces,1) < max_traces
                    these_traces = [these_traces; zeros(max_traces - size(these_traces,1),size(these_traces,2))];
                elseif size(these_traces,1) > max_traces
                    these_traces = these_traces(1:max_traces,:);    
                end
                these_traces_offset = get_trace_stack(these_traces,size(these_traces,2)-1,100,downsample_rate);
                if plot_avg
                    these_traces_offset = mean(these_traces_offset);
                end
                time = (1:size(these_traces_offset,2))*downsample_rate + (i-1)*(size(these_traces_offset,2)*downsample_rate + grid_offset_x);

                plot(repmat(time',1,size(these_traces_offset,1)),these_traces_offset' + grid_offset_y(j),'color',colors(array_i,:))
                hold on;
                if i == 1 && array_i == 1;
                    grid_offset_y(j+1) = grid_offset_y(j) - 600 - grid_offset_y_spacer;
%                     grid_offset_y(j+1) = grid_offset_y(j) - (max(max(these_traces_offset)) - min(min(these_traces_offset))) - grid_offset_y_spacer;
                end
    %         elseif j ~= 1 && i == 1
    %             grid_offset_y(j+1) = grid_offset_y(j) + (grid_offset_y(j) - grid_offset_y(j-1)) + grid_offset_y_spacer;
    %         elseif i == 1
    %             grid_offset_y(j+1) = grid_offset_y(j) + (max_traces + 1)*grid_offset_y_spacer;
            elseif i == 1 && array_i == 1;
                 grid_offset_y(j+1) = grid_offset_y(j) - 600 - grid_offset_y_spacer;
            end
        end

    end


    axis tight
    axis off
    
    if exist('plot_names','var') && i <= length(plot_names)
        title(plot_names{array_i})
    end
    
end
    hold off
% linkaxes(all_axes,'xy')

