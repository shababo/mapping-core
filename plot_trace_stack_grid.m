function plot_trace_stack_grid(traces_array,in_max_traces,downsample_rate,plot_avg,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    grid_colors = varargin{1};
end

if length(varargin) > 1 && ~isempty(varargin{2})
    in_axes = varargin{2};
    axes(in_axes);
end

if length(varargin) > 2 && ~isempty(varargin{3})
    alphas = varargin{3};
end

if length(varargin) > 3 && ~isempty(varargin{4})
    events = varargin{4};
end


[num_rows, num_cols] = size(traces_array);

count = 1;
grid_offset_y_spacer = 200;
grid_offset_y = -grid_offset_y_spacer;
grid_offset_x = .015;

if isinf(in_max_traces)
    max_traces = 0;
    for i = 1:num_cols
        for j = 1:num_rows
            num_traces = size(traces_array{j,i},1);
            if num_traces > max_traces
                max_traces = num_traces;
            end
        end
    end
else
    max_traces = in_max_traces;
end

trace_stack_offset = 0;
for i = 1:num_cols
    
    
    
    for j = 1:num_rows
        
        if ~isempty(traces_array{j,i})
            these_traces = traces_array{j,i};
            if size(these_traces,1) < max_traces
                these_traces = [these_traces; zeros(max_traces - size(these_traces,1),size(these_traces,2))];
            elseif size(these_traces,1) > max_traces
                these_traces = these_traces(1:max_traces,:);    
            end
            [these_traces_offset, offsets] = get_trace_stack(these_traces,size(these_traces,2)-1,25,downsample_rate);
            if plot_avg
                these_traces_offset = mean(these_traces_offset);
            end
            time = (1:size(these_traces_offset,2))/20000*downsample_rate + (i-1)*(size(these_traces_offset,2)/20000*downsample_rate + grid_offset_x);

            if exist('grid_colors','var')
                if grid_colors.color_i(j,i) == 0
                    this_color = [0 0 0];
                else
                    this_color_i = fix((grid_colors.color_i(j,i)-grid_colors.clims(1))/(grid_colors.clims(2)-grid_colors.clims(1))*size(grid_colors.colormap,1))+1;
                    this_color = grid_colors.colormap(min(this_color_i,size(grid_colors.colormap,1)),:);
                end
            elseif exist('alphas','var')
                this_color = [0 0 0 alphas(j,i)];
            else
                this_color = [0 0 0];
            end
            plot(repmat(time',1,size(these_traces_offset,1)),these_traces_offset' + grid_offset_y(j),'Color',this_color)
            hold on
            if exist('events','var')
                these_events = events{j,i};
                event_times = [];
                event_pos = [];
                for ii = 1:length(offsets)
                    event_times = [event_times these_events{ii}];
                    event_pos = [event_pos (offsets(ii) + grid_offset_y(j))*ones(size(these_events{ii}))];
                end
                scatter(time(round(event_times)),event_pos,[],'filled')
                hold on;
            end
            
            if i == 1
                grid_offset_y(j+1) = grid_offset_y(j) - (max(max(these_traces_offset)) - min(min(these_traces_offset))) - grid_offset_y_spacer;
%                 grid_offset_y(j+1) = grid_offset_y(j) - trace_stack_offset - grid_offset_y_spacer;
            end
%         elseif j ~= 1 && i == 1
%             grid_offset_y(j+1) = grid_offset_y(j) + (grid_offset_y(j) - grid_offset_y(j-1)) + grid_offset_y_spacer;
%         elseif i == 1
%             grid_offset_y(j+1) = grid_offset_y(j) + (max_traces + 1)*grid_offset_y_spacer;
        elseif i == 1
             grid_offset_y(j+1) = grid_offset_y(j) - trace_stack_offset - grid_offset_y_spacer;
        end
    end
    
end
hold off
axis tight
axis off


