function corr_image = get_corr_image(trace_array,do_plot)

corr_image = zeros(size(trace_array));

for i = 1:size(trace_array,1)
    for j = 1:size(trace_array,2)
        
        if ~isempty(trace_array{i,j})

%             mean_trace = mean(trace_array{i,j},1);
%             if min(mean_trace) < -1000
%                 mean_trace = mean(trace_array{i,j}([1 2],:),1);
%             end

        traces = trace_array{i,j};
        corrs = corr(traces');
        unique_corrs = [corrs(1,2) corrs(1,3) corrs(2,3)];
        corr_image(i,j) = mean(unique_corrs);
%         switch min_or_max
%                 case 'min'
%                     corr_image(i,j) = mean_trace(.005*20000) - min(mean_trace);
%                         
%                 case 'max'
%                     maxes = max(trace_array{i,j},[],2);
%                     starts = trace_array{i,j}(:,.005*20000);
%                     corr_image(i,j) = mean(max(traces));
%                     corr_image(i,j) =  mean(maxes - starts);
%             end
        else
            corr_image(i,j) = NaN;
        end
        
    end
end

if do_plot
   
    imagesc(corr_image)
    colormap hot
    colorbar
end