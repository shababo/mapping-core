function plot_timeseries_map(traces_array,in_max_traces,downsample_rate,plot_avg,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    color_map = varargin{1};
end

if length(varargin) > 1 && ~isempty(varargin{2})
    linewidth_map = varargin{2};
end

% if length(varargin) > 1 && ~isempty(varargin{2})
%     in_axes = varargin{2};
%     axes(in_axes);
% end

if length(varargin) > 2 && ~isempty(varargin{3})
    cell_map = varargin{3};
end

if length(varargin) > 3 && ~isempty(varargin{4})
    events = varargin{4};
end

if length(varargin) > 4 && ~isempty(varargin{5})
    loc_names = varargin{5};
else
    loc_names = [];
end


[num_rows, num_cols] = size(traces_array);

count = 1;
grid_offset_y_spacer = 75;
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

% for i = 1:num_cols
%     for j= 1:num_rows
%         if exist('cell_map','var') && cell_map(j,i)
%             
%             scatter((i-1)*(500/20000*downsample_rate + grid_offset_x),grid_offset_y(j),2000,[1 .6 .6],'filled')
%             hold on
%         end
%     end
% end


for i = 1:num_cols
    
    
    
    for j = 1:num_rows
        
        if exist('cell_map','var') && cell_map(j,i)
            
            scatter((i-1)*(500/20000*downsample_rate + grid_offset_x),grid_offset_y(j),2000,[1 0 cell_map(j,i)/max(cell_map(:))],'filled')
%             text((i-1)*(500/20000*downsample_rate + grid_offset_x),grid_offset_y(j),txt1)
            
            hold on
        end
        
            
        if ~isempty(traces_array{j,i})
            these_traces = traces_array{j,i};
            if size(these_traces,1) < max_traces
%                 these_traces = [these_traces; zeros(max_traces - size(these_traces,1),size(these_traces,2))];
            elseif size(these_traces,1) > max_traces
                these_traces = these_traces(1:max_traces,:);    
            end
%             for k = 1:size(these_traces,1)
%                 these_traces(k,:) = highpass_filter(these_traces(k,:),20000);
%             end
            [these_traces_offset, offsets] = get_trace_stack(these_traces,size(these_traces,2)-1,25,downsample_rate);
            
            time = (1:size(these_traces_offset,2))/20000*downsample_rate + (i-1)*(size(these_traces_offset,2)/20000*downsample_rate + grid_offset_x);

            if exist('color_map','var')
                these_colors = color_map{j,i};
            else
                these_colors = zeros(size(these_traces,1),3);
            end
            
            if exist('linewidth_map','var')
                these_lws = linewidth_map{j,i};
            else
                these_lws = zeros(size(these_traces,1),1);
            end
            these_lws(isinf(these_lws)) = 0.01;
%             elseif exist('alphas','var')
%                 this_color = [0 0 0 alphas(j,i)];
%             else
%                 this_color = [0 0 0];
%             end
%             for k = 1:size(these_traces_offset,1)
%                 plot(time,these_traces_offset(k,:) + grid_offset_y(j),'Color',these_colors(k,:))
%                 hold on
%             end

%             plot(repmat(time',1,size(these_traces_offset,1)),these_traces_offset' + grid_offset_y(j),'Color',this_color)
            for k = 1:size(these_traces_offset,1)
                plot(time,these_traces_offset(k,:) + grid_offset_y(j),'Color',these_colors(k,:),'Linewidth',these_lws(k))
            end
            hold on
            if plot_avg
                this_color = [0 0 0];
                these_traces_offset = mean(these_traces_offset);

            
                plot(repmat(time',1,size(these_traces_offset,1)),these_traces_offset' + grid_offset_y(j),'Color',this_color,'linewidth',2)
            end
            hold on
            x1 = time(20);
            y1 = these_traces_offset(1,20)' + grid_offset_y(j) + 50;
            if isempty(loc_names)
                txt1 = [num2str(j - 150) ', ' num2str(i - 150)];
            else
                txt1 = '';loc_names{j,i};
            end
            text(x1,y1,txt1)
            if exist('events','var')
                these_events = events{j,i};
                event_times = [];
                event_pos = [];
                if ~ isempty(these_events)
%                     length(offsets)
                    for ii = 1:length(offsets)
                        if iscell(these_events)
                            event_times = [event_times these_events{ii}];
                            event_pos = [event_pos (offsets(ii) + grid_offset_y(j))*ones(size(these_events{ii}))];
                        else
                            this_struct = these_events(ii); 
                            event_times = [event_times this_struct.times];
                            event_pos = [event_pos (offsets(ii) + grid_offset_y(j))*ones(size(this_struct.times))];
                        end
                    end
                    event_pos(event_times > length(time)-20) = [];
                    event_times(event_times > length(time)-20) = [];
                    event_pos(event_times < 20) = [];
                    event_times(event_times < 20) = [];
                    scatter(time(round(event_times)),event_pos+10,20*ones(size(event_pos)),[.0 .5 1],'filled')
                    hold on
%                     histtime = downsample(time,20);
%                     histedge = 0:20:length(time);
%                     event_counts = histcounts(event_times,histedge)*5;
%                     hist_pos = max(these_traces_offset' + grid_offset_y(j)) + 10;
%                     assignin('base','histtime',histtime)
%                     assignin('base','event_counts',event_counts)
%                     plot(histtime(1:end-1),event_counts + hist_pos,'linewidth',2)
                end
                hold on
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


