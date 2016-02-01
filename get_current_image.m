function current_image = get_current_image(trace_array,do_plot)

current_image = zeros(size(trace_array));

for i = 1:size(trace_array,1)
    for j = 1:size(trace_array,2)
        
        if ~isempty(trace_array{i,j})
            mean_trace = mean(trace_array{i,j});
            current_image(i,j) = mean_trace(.005*20000) - min(mean_trace);
        else
            current_image(i,j) = NaN;
        end
        
    end
end

if do_plot
    figure
    imagesc(current_image)
    colorbar
end