function [trace_grid] = get_mapping_data(filename,run,min_or_max,plot_grid,plot_image,varargin)%,max_traces,trace_limits, plot_grid,min_or_max,plot_image)

%% corr images - 2 cells
% 
% trace_grid = build_trace_grid(filename, run, 10, max_traces, trace_limits, plot_grid);
% 
% current_image = get_current_image(trace_grid,min_or_max,plot_image);


figure
subplot(221)
varargin{:}

trace_grid_ch1 = build_trace_grid(filename, run, 1, 3, [], 1,[],0,varargin{:});
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

trace_grid_ch2 = build_trace_grid(filename, run, 2, 3, [], 1,[],0,varargin{:});
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

disp('corr of images: ')
corr([corr_image_ch1(:) corr_image_ch2(:)])


%% 1 cell multi z spiking
% figure
% varargin{:}
% 
% trace_grid = build_trace_grid(filename, run, 1, 3, [], 1,[],1,varargin{:});

