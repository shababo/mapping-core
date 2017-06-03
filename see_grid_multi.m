function [maps,map_index,corr_maps,stddev_maps] = see_grid_multi(traces,sequence,stim_key,spacing,show_raw_data,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    do_std_map = varargin{1};
else
    do_std_map = 0;
end

if length(varargin) > 1 && ~isempty(varargin{2})
    do_corr_map = varargin{1};
else
    do_corr_map = 0;
end

[maps, map_index] = build_slm_maps_multi(...
    traces,sequence,stim_key,spacing,[-150 150],[-150 150]);

if show_raw_data
    figure
%     subplot(121)
    plot_trace_stack_grid(maps{1},Inf,5,0);
    title(['Power = ' num2str(sequence(1).target_power) ' mW'])
%     subplot(223)
%     plot_trace_stack_grid(maps{2},Inf,5,0);
end

if do_std_map
    stddev_maps{1} = get_stdev_map(trace_array,do_plot,neighborhood_size);
    stddev_maps{2} = get_stdev_map(trace_array,do_plot,neighborhood_size);
else
    stddev_maps = cell();
end

if do_corr_map
    corr_maps{1} = get_corr_image(maps{1},0,0);
    corr_maps{2} = get_corr_image(maps{2},0,0);
else
    corr_maps = cell();
end

if do_corr_map
figure
subplot(221); 
imagesc(corr_maps{1}); caxis([0 0.5])
title(['Cell 1 Corr Map: Power = ' num2str(sequence(1).target_power) ' mW'])
figure
subplot(221); 
imagesc(stddev_maps{1}); %caxis([0 1])
title(['Cell 1 Corr Map: Power = ' num2str(sequence(1).target_power) ' mW'])
subplot(221); 
imagesc(corr_maps{2}); caxis([0 0.5])
title(['Cell 2 Corr Map: Power = ' num2str(sequence(1).target_power) ' mW'])
figure
subplot(221); 
imagesc(stddev_maps{2}); %caxis([0 1])
title(['Cell 2 Corr Map: Power = ' num2str(sequence(1).target_power) ' mW'])

end

% assignin('base','maps',maps)
