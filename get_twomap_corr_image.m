function corr_image = get_twomap_corr_image(trace_array1, trace_array2, do_plot)

corr_image = zeros(size(trace_array1));

for i = 1:size(trace_array1,1)
    for j = 1:size(trace_array1,2)
        
        if ~isempty(trace_array1{i,j}) || ~isempty(trace_array2{i,j}) && isequal(size(trace_array1{i,j}),size(trace_array2{i,j}))

%             mean_trace = mean(trace_array1{i,j},1);
%             if min(mean_trace) < -1000
%                 mean_trace = mean(trace_array1{i,j}([1 2],:),1);
%             end
        
            traces1 = trace_array1{i,j};
            traces2 = trace_array2{i,j};
%             num_traces = size(traces1,1);
%             
%             for k = 1:num_traces
%               
%                corr_image(i,j) = corr(traces1(k,:)',traces2(k,:)')/num_traces;
%                                
%             end

            tmp = cov(mean(traces1)',mean(traces2)');
            corr_image(i,j) = tmp(1,2);
            

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