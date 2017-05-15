function [maps, map_index] = build_slm_maps_multi(traces,sequence,stim_key,spacing)

num_cells = length(traces);
maps = cell(num_cells,1);

stim_key_bin = round(stim_key/spacing)*spacing;

% x_bins = unique(stim_key_bin(:,1,:));
% y_bins = unique(stim_key_bin(:,2,:));
% z_bins = unique(stim_key_bin(:,3,:));

stim_x_min = min(min(stim_key_bin(:,1,:)));
stim_x_max = max(max(stim_key_bin(:,1,:)));

stim_y_min = min(min(stim_key_bin(:,2,:)));
stim_y_max = max(max(stim_key_bin(:,2,:)));

stim_z_min = 0;%min(min(stim_key_bin(:,3,:)));
stim_z_max = 0;%max(max(stim_key_bin(:,3,:)));

x_bins = stim_x_min:spacing:stim_x_max;
y_bins = stim_y_min:spacing:stim_y_max;
z_bins = stim_z_min:spacing:stim_z_max;

grid_dims = [length(x_bins) length(y_bins) length(z_bins)];

min_bin = [stim_x_min stim_y_min stim_z_min];
map_index = (bsxfun(@minus,stim_key_bin,min_bin) + spacing)/spacing;
map_index(:,3,:) = 1;

num_traces = length(sequence);
% size(traces,1)
assignin('base','map_index',map_index)


for i = 1:num_cells
    
    maps{i} = cell(grid_dims);
    these_traces = traces{i};
%     size(these_traces)
    for j = 1:num_traces
        j_stim = sequence(j).precomputed_target_index;
        for k = 1:size(map_index,3)
%             i
%             j
%             j_stim
%             k
            if isnan(map_index(j_stim,1,k))
                break
            end
            maps{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)} = ...
                [maps{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)}; ...
                these_traces(j,:)];
            
        end
    end
end
