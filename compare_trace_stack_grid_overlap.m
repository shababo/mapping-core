function compare_trace_stack_grid_overlap(traces_arrays,in_max_traces,downsample_rate,in_axes,plot_avg, varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    plot_names = varargin{1};
end

if length(varargin) > 1 && ~isempty(varargin{2})
    num_plot_rows = varargin{2};
else
    num_plot_rows = 1;
end

if length(varargin) > 2 && ~isempty(varargin{3})
    alphas = varargin{3};
end

num_arrays = length(traces_arrays);
all_axes = [];

if isempty(in_axes)
    ax = gca;

else
    ax = in_axes;

end

axes(ax)

colors = [39 168 224; .6*255 0 255]/255;

for array_i = 1:num_arrays
    

    
    this_array = traces_arrays{array_i};
    
    [num_rows, num_cols] = size(this_array);
    if isinf(in_max_traces)
        max_traces = 0;
        for i = 1:num_cols
            for j = 1:num_rows
                num_traces = size(this_array{j,i},1);
                if num_traces > max_traces
                    max_traces = num_traces;
                end
            end
        end
    else
        max_traces = in_max_traces;
    end
    max_traces
    
    count = 1;
    grid_offset_y_spacer = 400*max_traces;
    if array_i == 1
        grid_offset_y = -grid_offset_y_spacer;
    end
    grid_offset_x = .01*20000;

    


    for i = 1:num_cols



        for j = 1:num_rows

            if ~isempty(this_array{j,i})
                these_traces = this_array{j,i};
%                 these_traces(:,1:20) = 0;
                if size(these_traces,1) < max_traces
                    these_traces = [these_traces; zeros(max_traces - size(these_traces,1),size(these_traces,2))];
                elseif size(these_traces,1) > max_traces
%                     if j > 9
%                     these_traces = these_traces([1 3 6],:);  
%                     else
%                       these_traces = these_traces([2 3 6],:);   
%                     end
                    these_traces = these_traces(1:max_traces,:);
                end
%                 size_tmp = size(these_traces);
%                 these_traces = these_traces(:);
%                 these_traces = zscore(these_traces);
%                 these_traces = reshape(these_traces,size_tmp(1),size_tmp(2));
%                 these_traces = zscore(these_traces,0,2);
                these_traces_offset = get_trace_stack(these_traces,size(these_traces,2)-1,200,downsample_rate);
                if plot_avg
                    these_traces_offset = mean(these_traces_offset);
                end
                time = (1:size(these_traces_offset,2))*downsample_rate + (i-1)*(size(these_traces_offset,2)*downsample_rate + grid_offset_x);
                if exist('alphas','var')
                    this_color = [colors(array_i,:) alphas{array_i}(j,i) > 0];
                else
                    this_color = colors(array_i,:);
                end
                plot(repmat(time',1,size(these_traces_offset,1)),these_traces_offset' + grid_offset_y(j),'color',this_color,'Linewidth',1)
                hold on;
                if i == 1 && array_i == 1;
                    grid_offset_y(j+1) = grid_offset_y(j) - grid_offset_y_spacer;
%                     grid_offset_y(j+1) = grid_offset_y(j) - (max(max(these_traces_offset)) - min(min(these_traces_offset))) - grid_offset_y_spacer;
                end
    %         elseif j ~= 1 && i == 1
    %             grid_offset_y(j+1) = grid_offset_y(j) + (grid_offset_y(j) - grid_offset_y(j-1)) + grid_offset_y_spacer;
    %         elseif i == 1
    %             grid_offset_y(j+1) = grid_offset_y(j) + (max_traces + 1)*grid_offset_y_spacer;
            elseif i == 1 && array_i == 1;
                 grid_offset_y(j+1) = grid_offset_y(j) - grid_offset_y_spacer;
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

