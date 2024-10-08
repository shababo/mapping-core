
function [maps, mpp_maps, map_index, color_maps] = build_slm_maps_multi(traces,mpp,sequence,stim_key,spacing,varargin)

assignin('base','mpp_build_slm',mpp)
if length(varargin) > 2 && ~isempty(varargin{3})
    these_colors = varargin{3};
else
    for i = 1:length(traces)        
        these_colors{i} = zeros(size(traces{i},1),3);
    end
end

num_cells = length(traces);
maps = cell(num_cells,1);
mpp_maps = cell(num_cells,1);
color_maps = cell(num_cells,1);

spacing = [spacing spacing -1];

stim_z_min = floor(min(min(stim_key(:,3,:))));
stim_z_max = ceil(max(max(stim_key(:,3,:))));
spacing(3) = ceil((stim_z_max - stim_z_min)/3);

stim_key_bin = bsxfun(@times,round(bsxfun(@rdivide,stim_key,spacing)),spacing);

% x_bins = unique(stim_key_bin(:,1,:));
% y_bins = unique(stim_key_bin(:,2,:));
% z_bins = unique(stim_key_bin(:,3,:));

if ~isempty(varargin) && ~isempty(varargin{1})
    x_range = varargin{1};
    stim_x_min = x_range(1);
    stim_x_max = x_range(2);
else
    stim_x_min = min(min(stim_key_bin(:,1,:)));
    stim_x_max = max(max(stim_key_bin(:,1,:)));
end


if length(varargin) > 1 && ~isempty(varargin{2})
    y_range = varargin{2};
    stim_y_min = y_range(1);
    stim_y_max = y_range(2);
else
    stim_y_min = min(min(stim_key_bin(:,2,:)));
    stim_y_max = max(max(stim_key_bin(:,2,:)));
end
% 
% stim_y_min = min(min(stim_key_bin(:,2,:)));
% stim_y_max = max(max(stim_key_bin(:,2,:)));

stim_z_min = floor(min(min(stim_key_bin(:,3,:))));
stim_z_max = ceil(max(max(stim_key_bin(:,3,:))));
% spacing(3) = ceil((stim_z_max - stim_z_min)/9);

x_bins = stim_x_min:spacing(1):stim_x_max;
y_bins = stim_y_min:spacing(2):stim_y_max;
z_bins = stim_z_min:spacing(3):stim_z_max;

grid_dims = [length(x_bins) length(y_bins) length(z_bins)];

min_bin = [stim_x_min stim_y_min stim_z_min];
map_index = bsxfun(@rdivide,bsxfun(@plus,bsxfun(@minus,stim_key_bin,min_bin),spacing),spacing);

% map_index(:,3,:) = 1;


num_traces = length(sequence);
% size(traces,1)
% assignin('base','map_index',map_index)


for i = 1:num_cells
    
    maps{i} = cell(grid_dims);
    mpp_maps{i} = cell(grid_dims);
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
            if isnan(map_index(j_stim,1,k))
                break
            end
            
            j_stim
            k
            map_index(j_stim,:)
            if isnan(map_index(j_stim,3))
                map_index(j_stim,3) = 1;
            end
            maps{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)} = ...
                [maps{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)}; ...
                these_traces(j,:)];

            color_maps{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)} = ...
                [color_maps{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)}; ...
                this_color(j,:)];

            if ~isempty(mpp{i})
                mpp_maps{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)} = ...
                    [mpp_maps{i}{map_index(j_stim,1,k),map_index(j_stim,2,k),map_index(j_stim,3,k)}; ...
                    mpp{i}(j)];
            end
        end
    end
end
