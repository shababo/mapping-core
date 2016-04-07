function variance_image = get_corr_image(trace_array,min_or_max,do_plot)

variance_image = zeros(size(trace_array));

for i = 1:size(trace_array,1)
    for j = 1:size(trace_array,2)
        
        if ~isempty(trace_array{i,j})

            mean_trace = mean(trace_array{i,j},1);
            if min(mean_trace) < -1000
                mean_trace = mean(trace_array{i,j}([1 2],:),1);
            end
            switch min_or_max
                case 'min'
                    variance_image(i,j) = mean_trace(.005*20000) - min(mean_trace);
                        
                case 'max'
                    maxes = max(trace_array{i,j},[],2);
                    starts = trace_array{i,j}(:,.005*20000);
%                     variance_image(i,j) = mean(max(traces));
                    variance_image(i,j) =  mean(maxes - starts);
            end
        else
            variance_image(i,j) = NaN;
        end
        
    end
end

if do_plot
   
    imagesc(variance_image)
    colormap hot
    colorbar
end