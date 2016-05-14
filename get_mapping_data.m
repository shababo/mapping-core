function [trace_grid, current_image, corr_image] = get_mapping_data(filename,run,min_or_max,plot_grid,plot_image,varargin)%,max_traces,trace_limits, plot_grid,min_or_max,plot_image)

% trace_grid = build_trace_grid(filename, run, 10, max_traces, trace_limits, plot_grid);
% 
% current_image = get_current_image(trace_grid,min_or_max,plot_image);
figure
subplot(221)
 varargin{:}

trace_grid_ch1 = build_trace_grid(filename, run, 1, 3, [], 1,[],varargin{:});
trace_grid = trace_grid_ch1;
title('Cell 1')
% subplot(132)
% current_image = get_current_image(trace_grid,min_or_max,1);
% title('amplitude')
% axis off
% surf1 = gca;
subplot(222)

corr_image_ch1 = get_corr_image(trace_grid_ch1,plot_image);
current_image = corr_image_ch1;
title('local correlations')
axis off

subplot(223)
 varargin{:}

trace_grid_ch2 = build_trace_grid(filename, run, 2, 5, [], 1,[],varargin{:});
title('Cell 2')
% subplot(132)
% current_image = get_current_image(trace_grid,min_or_max,1);
% title('amplitude')
% axis off
% surf1 = gca;
subplot(224)

corr_image_ch2 = get_corr_image(trace_grid_ch2,plot_image);
current_image = corr_image_ch2;
title('local correlations')
axis off


% caxis([0 1])
surf2 = gca;
% subplot(224)
% surf(1:21,21:-1:1,current_image .* (corr_image > 2.5*std(corr_image(:))))
% imagesc(current_image .* corr_image)
surf3 = gca;
% 
% color_base = sign(current_image .* corr_image);
% color_matrix = color_base .* (color_base > 10) .* imregionalmax(corr_image);%) .* (current_image > 2.0*std(current_image(:)));
% % grid_color.clims=get(surf2,'CLim');
% grid_color.clims = [min(color_matrix(color_matrix > 0)) max(color_matrix(color_matrix > 0))];
% grid_color.colormap= colormap(hot(64));
% grid_color.color_i = color_matrix;
% % grid_color.color_i = corr_image .* (corr_image > 2.5*std(corr_image(:)));
% 
% subplot(221)
% 
% trace_grid = build_trace_grid(filename, run, 3, [], plot_grid,grid_color);
% linkaxes([surf1 surf2 surf3])