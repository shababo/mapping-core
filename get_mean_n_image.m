function current_image = get_mean_n_image(mean_n_grid,do_plot)

current_image = zeros(size(mean_n_grid));

for i = 1:size(mean_n_grid,1)
    for j = 1:size(mean_n_grid,2)
        
        if ~isempty(mean_n_grid{i,j})
            current_image(i,j) = mean(mean_n_grid{i,j});
        else
            current_image(i,j) = NaN;
        end
        
    end
end

if do_plot
    imagesc(current_image)
    colorbar
end