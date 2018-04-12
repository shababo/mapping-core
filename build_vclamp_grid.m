function [vclamp_map, psc_time_map, color_map, linewidth_map, cell_map] = build_vclamp_grid(experiment_setup,trials,spacing,varargin)

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

cell_map = zeros([grid_dims 3]);
vclamp_map = cell(grid_dims);
psc_time_map = cell(grid_dims);
color_map = cell(grid_dims);
linewidth_map = cell(grid_dims);

min_bin = [experiment_setup.neighbourhood_params.x_bounds(1) experiment_setup.neighbourhood_params.y_bounds(1) 0];
targeted_cell_IDs = [];
% cell_colors = zeros(length(experiment_setup.neurons),3);

for j = 1:num_trials

    targeted_cell_IDs = union(targeted_cell_IDs,trials(j).cell_IDs);
    
    if size(trials(j).locations,1) == 1
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


    %         for i = 1:length(trials(j).cell_IDs)
    %             if ~isnan(trials(j).cell_IDs(i))
    %                 cell_colors(trials(j).cell_IDs(i),:) = cell_color;
    %             end
    %         end
            switch trials(j).location_IDs(k)
                case 1
                    this_color = [0 0 0];
                otherwise
                    this_color = [.4 .4 .4];
            end
            color_map{map_index(1),map_index(2),map_index(3)} = ...
                [color_map{map_index(1),map_index(2),map_index(3)}; ...
                this_color];

            linewidth_map{map_index(1),map_index(2),map_index(3)} = ...
                [linewidth_map{map_index(1),map_index(2),map_index(3)}; ...
                trials(j).power_levels(k)/30];

            if isempty(psc_time_map{map_index(1),map_index(2),map_index(3)})
                psc_time_map{map_index(1),map_index(2),map_index(3)} = cell(0);
            end
            psc_time_map{map_index(1),map_index(2),map_index(3)}{end+1} = trials(j).event_times;

        end
    end
end

targeted_cell_IDs(isnan(targeted_cell_IDs)) = [];
assignin('base','targeted_cell_IDs',targeted_cell_IDs)

for j = 1:length(targeted_cell_IDs)
    
    i_neuron = [experiment_setup.neurons.cell_ID] == targeted_cell_IDs(j);
    this_loc = round(experiment_setup.neurons(i_neuron).location/spacing)*spacing;
    map_index = (this_loc - min_bin + spacing)/spacing;
    map_index(3) = 1;
    if any(isnan(map_index) | map_index < 1 | map_index > [300 300 100])
        continue
    end
    
    switch experiment_setup.neurons(i_neuron).group_ID
        case 'undefined'
            cell_color = rgb('DarkGray');
        case 'connected'
            cell_color = rgb('ForestGreen');
        case 'alive'
            cell_color = rgb('BlueViolet');
        case 'disconnected'
            cell_color = rgb('DarkRed');
        case 'secondary'
            cell_color = rbg('Snow');
        otherwise
            cell_color = zeros(1,3);
    end

    cell_map(map_index(1),map_index(2),map_index(3),:) = cell_color;


    
end
