function stacked_grid = stack_grids(grids)

num_grids = length(grids);

num_x_locs = size(grids{1},1);
num_y_locs = size(grids{1},2);

stacked_grid = cell(num_x_locs,num_y_locs);

for i = 1:num_x_locs
    for j = 1:num_y_locs
        
        for k = 1:num_grids
            stacked_grid{i,j} = [stacked_grid{i,j}; grids{k}{i,j}];
        end
    end
end
        
            