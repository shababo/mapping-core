function current_image = get_current_image(trace_array,min_or_max,do_plot)

current_image = zeros(size(trace_array));

for i = 1:size(trace_array,1)
    for j = 1:size(trace_array,2)
        
        if ~isempty(trace_array{i,j})

            mean_trace = mean(trace_array{i,j},1);
            if min(mean_trace) < -1000
                mean_trace = mean(trace_array{i,j}([1 2],:),1);
            end
            these_traces = [];
            neighborhood_size = 1;
            for ii = -neighborhood_size:neighborhood_size
                for jj = -neighborhood_size:neighborhood_size
                     if ~(i+ii < 1 || i+ii > size(trace_array,1) || j+jj < 1 || j+jj > size(trace_array,2))
                         these_traces = [these_traces; trace_array{i+ii,j+jj}];
                     end
                end
            end
            switch min_or_max
                case 'min'
                    mins = min(these_traces,[],2);
                    current_image(i,j) = mean(median(these_traces,2) - mins);
                        
                case 'max'
                    maxes = max(these_traces,[],2);
%                     starts = these_traces(:,.005*20000);
%                     current_image(i,j) = mean(max(traces));
                    current_image(i,j) =  mean(maxes - median(these_traces,2));
            end
        else
            current_image(i,j) = NaN;
        end
        
    end
end

if do_plot
   
%     surf(1:size(current_image,1),size(current_image,2):-1:1,current_image)
    imagesc(current_image)
    colormap hot
    colorbar
end