function [vclamp_map, psc_time_map, color_map] = build_vclamp_grid(experiment_setup,trials,spacing,varargin)

num_trials = length(trials);

% if length(varargin) > 2 && ~isempty(varargin{3})
%     these_colors = varargin{3};
% else
%     for i = 1:num_trials        
%         these_colors{i} = zeros(1,3);
%     end
% end

x_bins = experiment_setup.neighbourhood_params.x_bounds(1):spacing:experiment_setup.neighbourhood_params.x_bounds(2);
y_bins = experiment_setup.neighbourhood_params.y_bounds(1):spacing:experiment_setup.neighbourhood_params.y_bounds(2);
z_bins = 0;%stim_z_min:spacing:stim_z_max;

grid_dims = [length(x_bins) length(y_bins) length(z_bins)];

vclamp_map = cell(grid_dims);
psc_time_map = cell(grid_dims);
color_map = cell(grid_dims);

min_bin = [experiment_setup.neighbourhood_params.x_bounds(1) experiment_setup.neighbourhood_params.y_bounds(1) 0];
    

for j = 1:num_trials
    for k = 1:size(trials(j).locations,1)
        
        this_loc = round(trials(j).locations(k,:)/spacing)*spacing;
        map_index = (this_loc - min_bin + spacing)/spacing;
        map_index(3) = 1;
        if any(isnan(map_index))
            continue
        end

        vclamp_map{map_index(1),map_index(2),map_index(3)} = ...
            [vclamp_map{map_index(1),map_index(2),map_index(3)}; ...
            trials(j).voltage_clamp];

        switch trials(j).group_ID
            case 'undefined'
                this_color = rgb('SlateGray');
            case 'connected'
                this_color = rgb('Olive');
            case 'alive'
                this_color = rgb('DarkBlue');
        end
        color_map{map_index(1),map_index(2),map_index(3)} = ...
            [color_map{map_index(1),map_index(2),map_index(3)}; ...
            this_color];
        
        if isempty(psc_time_map{map_index(1),map_index(2),map_index(3)})
            psc_time_map{map_index(1),map_index(2),map_index(3)} = cell(0);
        end
        psc_time_map{map_index(1),map_index(2),map_index(3)}{end+1} = trials(j).event_times;

    end
end

