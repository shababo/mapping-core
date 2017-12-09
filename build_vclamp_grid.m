
function [vclamp_map, psc_time_map, map_index, color_maps] = build_vclamp_grid(experiment_setup,trials,spacing,varargin)

num_trials = length(trials);

if length(varargin) > 2 && ~isempty(varargin{3})
    these_colors = varargin{3};
else
    for i = 1:num_trials        
        these_colors{i} = zeros(1,3);
    end
end

x_bins = stim_x_min:spacing:stim_x_max;
y_bins = stim_y_min:spacing:stim_y_max;
z_bins = 0;%stim_z_min:spacing:stim_z_max;

grid_dims = [length(x_bins) length(y_bins) length(z_bins)];

vclamp_map = cell(num_cells,1);
psc_time_map = cell(num_cells,1);
color_maps = cell(num_cells,1);




num_traces = length(sequence);
% size(traces,1)
% assignin('base','map_index',map_index)


for i = 1:num_cells
    
    vclamp_map{i} = cell(grid_dims);
    psc_time_map{i} = cell(grid_dims);
    color_maps{i} = cell(grid_dims);
    these_traces = traces{i};
%     this_mpp = mpp{i};
    this_color = these_colors{i};
%     size(these_traces)
    for j = 1:num_traces
        j_stim = sequence(j).precomputed_target_index;
        for k = 1:size(map_index,3)
%             i
%             j
%             j_stim
%             k

            min_bin = [stim_x_min stim_y_min stim_z_min];
map_index = (bsxfun(@minus,stim_key_bin,min_bin) + spacing)/spacing;

map_index(:,3,:) = 1;
            if isnan(map_index(j_stim,1,k))
                break
            end

            vclamp_map{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)} = ...
                [vclamp_map{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)}; ...
                these_traces(j,:)];

            color_maps{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)} = ...
                [color_maps{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)}; ...
                this_color(j,:)];

            if ~isempty(mpp{i})
                psc_time_map{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)} = ...
                    [psc_time_map{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)}; ...
                    mpp{i}(j)];
            end
        end
    end
end
