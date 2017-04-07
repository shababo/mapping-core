function [maps,corr_maps] = see_grid_multi(traces,sequence,stim_key,spacing,do_plot)



maps = build_slm_maps_multi(...
    traces,sequence,stim_key,spacing);

if do_plot
    figure
    subplot(221)
    plot_trace_stack_grid(maps{1},Inf,1,0);
    subplot(223)
    plot_trace_stack_grid(maps{2},Inf,1,0);
end

corr_maps{1} = get_corr_image(maps{1},0,0);
corr_maps{2} = get_corr_image(maps{2},0,0);


if do_plot

subplot(222); imagesc(corr_maps{1}); %caxis([0 1])
subplot(224); imagesc(corr_maps{2}); %caxis([0 1])

end
