function std_dev_map = get_stdev_map(trace_array,do_plot,neighborhood_size)

std_dev_map = zeros(size(trace_array));

for i = 1:size(trace_array,1)
    for j = 1:size(trace_array,2)
        
        if ~isempty(trace_array{i,j})

%             mean_trace = mean(trace_array{i,j},1);
%             if min(mean_trace) < -1000
%                 mean_trace = mean(trace_array{i,j}([1 2],:),1);
%             end
        traces = [];
%         neighborhood_size = 0;
        for ii = -neighborhood_size:neighborhood_size
            for jj = -neighborhood_size:neighborhood_size
                if ~(i+ii < 1 || i+ii > size(trace_array,1) || j+jj < 1 || j+jj > size(trace_array,2))
                    traces = [traces; trace_array{i+ii,j+jj}];
                end
            end
        end
        these_stds = std(these_traces,[],2);
        std_dev_map(i,j) = mean(these_stds);
%         switch min_or_max
%                 case 'min'
%                     std_dev_map(i,j) = mean_trace(.005*20000) - min(mean_trace);
%                         
%                 case 'max'
%                     maxes = max(trace_array{i,j},[],2);
%                     starts = trace_array{i,j}(:,.005*20000);
%                     std_dev_map(i,j) = mean(max(traces));
%                     std_dev_map(i,j) =  mean(maxes - starts);
%             end
        else
            std_dev_map(i,j) = NaN;
        end
        
    end
end

if do_plot
   
%     surf(1:size(std_dev_map,1),size(std_dev_map,2):-1:1,std_dev_map)
    imagesc(std_dev_map)
    colormap hot
    colorbar
end