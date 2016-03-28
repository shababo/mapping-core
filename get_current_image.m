function current_image = get_current_image(trace_array,min_or_max,do_plot)

current_image = zeros(size(trace_array));

for i = 1:size(trace_array,1)
    for j = 1:size(trace_array,2)
        
        if ~isempty(trace_array{i,j})

            mean_trace = mean(trace_array{i,j},1);
            if min(mean_trace) < -1000
                mean_trace = mean(trace_array{i,j}([1 2],:),1);
            end
            switch min_or_max
                case 'min'
                    current_image(i,j) = mean_trace(.005*20000) - min(mean_trace);
                        
                case 'max'
%                     current_image(i,j) = mean(max(traces));
                    current_image(i,j) =  max(mean_trace) - mean_trace(.005*20000);
            end
        else
            current_image(i,j) = NaN;
        end
        
    end
end

if do_plot
    figure
    imagesc(current_image)
    colormap hot
    colorbar
end