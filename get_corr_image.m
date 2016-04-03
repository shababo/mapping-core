function corr_image = get_corr_image(trace_array,do_plot)

corr_image = zeros(size(trace_array));

for i = 1:size(trace_array,1)
    for j = 1:size(trace_array,2)
        
        if ~isempty(trace_array{i,j})

%             mean_trace = mean(trace_array{i,j},1);
%             if min(mean_trace) < -1000
%                 mean_trace = mean(trace_array{i,j}([1 2],:),1);
%             end
        traces = [];
        neighborhood_size = 0;
        for ii = -neighborhood_size:neighborhood_size
            for jj = -neighborhood_size:neighborhood_size
                if ~(i+ii < 1 || i+ii > size(trace_array,1) || j+jj < 1 || j+jj > size(trace_array,2))
                    traces = [traces; trace_array{i+ii,j+jj}];
                end
            end
        end
        corrs = triu(corr(traces(:,.008*20000:end)'),1);
        unique_corrs = unique(corrs(:));
        unique_corrs = unique_corrs(2:end); % get rid of leading 0
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
   
%     surf(1:size(corr_image,1),size(corr_image,2):-1:1,corr_image)
    imagesc(corr_image)
    colormap hot
    colorbar
end