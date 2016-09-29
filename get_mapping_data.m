function [trace_grids, current_image] = get_mapping_data(filename,run,min_or_max,plot_grid,plot_image,varargin)%,max_traces,trace_limits, plot_grid,min_or_max,plot_image)

%% corr images - 2 cells

figure
subplot(121)
varargin{:}

trace_grid_ch1 = build_trace_grid(filename, run, 1, 20, [], 1,[],0,varargin{:});
title('Cell 1')
subplot(122)
current_image = get_current_image(trace_grid_ch1,min_or_max,1);
title('amplitude')
axis off
colorbar
% surf1 = gca;
% subplot(222)
% % % 
corr_image_ch1 = get_corr_image(trace_grid_ch1,plot_image);
current_image = corr_image_ch1;
title('local correlations')
axis off
% % % 
% % % subplot(223)
% % %  varargin{:}
% % % 
% % % trace_grid_ch2 = build_trace_grid(filename, run, 2, 20, [], 1,[],0,varargin{:});
% % % title('Cell 2')
% % % % subplot(132)
% % % % current_image = get_current_image(trace_grid,min_or_max,1);
% % % % title('amplitude')
% % % % axis off
% % % % surf1 = gca;
% % % subplot(224)
% % % 
% % % corr_image_ch2 = get_corr_image(trace_grid_ch2,plot_image);
% % % current_image = corr_image_ch2;
% % % title('local correlations')
% % % axis off
% % % 
% % % disp('corr of images: ')
% % % corr([corr_image_ch1(:) corr_image_ch2(:)])
% % % 
% % % trace_grids = {trace_grid_ch1, trace_grid_ch2};

trace_grids = trace_grid_ch1;

% compare_trace_stack_grid(trace_grids,3,1,[subplot(221) subplot(223)],0,{'cell1','cell2'},1)

% set(gcf,'units','inches')
% set(gcf,'Position',[3.0521   -0.0521   14.6979   11.4062])
% 
% set(gcf,'PaperSize',[15 12])
% set(gcf,'PaperPosition',[.2   .2   14.6979   11.4062])
% 
% 
% fig_save_str = ['quickfigs/' filename(1:end-4) '_' num2str(run) '.pdf'];
% export_fig(fig_save_str,'-opengl')

% saveas(gcf,['quickfigs/' filename(1:end-4) '_' num2str(run) '.pdf'])


% saveas(gcf,['quickfigs/' filename(1:end-4) '_' num2str(run) '.jpg'])

        
%% 1 cell multi z spiking
% figure
% varargin{:}
% 
% trace_grids = build_trace_grid(filename, run, 2, 3, [], 1,[],1,varargin{:});

