function plot_trace_stack_grid(traces_array,in_max_traces,downsample_rate,plot_avg)


[num_rows, num_cols] = size(traces_array);

count = 1;
grid_offset_y_spacer = 150;
grid_offset_y = -grid_offset_y_spacer;
grid_offset_x = .01;

if isinf(in_max_traces)
    max_traces = 0;
    for i = 1:num_cols
        for j = 1:num_rows
            num_traces = size(traces_array{j,i},2);
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
        
        length(grid_offset_y)
        if ~isempty(traces_array{j,i})
            these_traces = traces_array{j,i};
            if size(these_traces,1) < max_traces
                these_traces = [these_traces; zeros(max_traces - size(these_traces,1),size(these_traces,2))];
            elseif size(these_traces,1) > max_traces
                these_traces = these_traces(1:max_traces,:);    
            end
            these_traces_offset = get_trace_stack(these_traces,size(these_traces,2)-1,25,downsample_rate);
            if plot_avg
                these_traces_offset = mean(these_traces_offset);
            end
            time = (1:size(these_traces_offset,2))/20000*downsample_rate + (i-1)*(size(these_traces_offset,2)/20000*downsample_rate + grid_offset_x);

            plot(repmat(time',1,size(these_traces_offset,1)),these_traces_offset' + grid_offset_y(j),'k')
            hold on;
            if i == 1
                grid_offset_y(j+1) = grid_offset_y(j) - (max(max(these_traces_offset)) - min(min(these_traces_offset))) - grid_offset_y_spacer;
            end
%         elseif j ~= 1 && i == 1
%             grid_offset_y(j+1) = grid_offset_y(j) + (grid_offset_y(j) - grid_offset_y(j-1)) + grid_offset_y_spacer;
%         elseif i == 1
%             grid_offset_y(j+1) = grid_offset_y(j) + (max_traces + 1)*grid_offset_y_spacer;
        elseif i == 1
             grid_offset_y(j+1) = grid_offset_y(j) - 600 - grid_offset_y_spacer;
        end
    end
    
end

hold off
axis tight
axis off