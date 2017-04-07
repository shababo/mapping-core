function maps = build_slm_maps_multi(traces,sequence,stim_key,spacing)

num_cells = size(traces,1);
maps = cell(num_cells,1);

stim_key_bin = floor(stim_key/spacing)*spacing;

% x_bins = unique(stim_key_bin(:,1,:));
% y_bins = unique(stim_key_bin(:,2,:));
% z_bins = unique(stim_key_bin(:,3,:));

stim_x_min = min(min(stim_key_bin(:,1,:)));
stim_x_max = max(max(stim_key_bin(:,1,:)));

stim_y_min = min(min(stim_key_bin(:,2,:)));
stim_y_max = max(max(stim_key_bin(:,2,:)));

stim_z_min = min(min(stim_key_bin(:,3,:)));
stim_z_max = max(max(stim_key_bin(:,3,:)));

x_bins = stim_x_min:spacing:stim_x_max;
y_bins = stim_y_min:spacing:stim_y_max;
z_bins = stim_z_min:spacing:stim_z_max;

grid_dims = [length(x_bins) length(y_bins) length(z_bins)];

min_bin = [stim_x_min stim_y_min stim_z_min];
map_index = (bsxfun(@minus,stim_key_bin,min_bin) + spacing)/spacing;


num_traces = length(sequence);



for i = 1:num_cells
    
    maps{i} = cell(grid_dims);
    
    for j = 1:num_traces
        for k = 1:size(map_index,1)
            
            maps{i}{map_index(j,1,k),map_index(j,2,k),map_index(j,3,k)} = ...
                [maps{j}{map_index(j,1,k),map_index(j,2,k),map_index(j,3,k)}; ...
                traces(j,:)];
            
        end
    end
end
