function [maps,mpp_maps,map_index,corr_maps] = see_grid_multi(traces,mpp,sequence,stim_key,spacing,do_plot)



[maps, mpp_maps, map_index] = build_slm_maps_multi(...
    traces,mpp,sequence,stim_key,spacing);

if do_plot
    figure
%     subplot(121)
    plot_trace_stack_grid(maps{1},Inf,1,0);
    title(['Power = ' num2str(sequence(1).target_power) ' mW'])
%     subplot(223)
%     plot_trace_stack_grid(maps{2},Inf,5,0);
end

corr_maps{1} = get_corr_image(maps{1},0,0);
% corr_maps{2} = get_corr_image(maps{2},0,0);


if do_plot
figure
% subplot(122); 
% subplot(3,3,subplot_ind)
% size(corr_maps{1})
% assignin('base','corr_map',corr_maps{1});
imagesc(corr_maps{1}); caxis([0 .5])
title(['Power = ' num2str(sequence(1).target_power) ' mW'])
% subplot(224); imagesc(corr_maps{2}); %caxis([0 1])

end

maps = maps{1};
corr_maps = corr_maps{1};
mpp_maps = mpp_maps{1};
% assignin('base','maps',maps)
