function [maps,mpp_maps,map_index,corr_maps,stddev_maps] = see_grid_multi(traces,mpp,sequence,stim_key,spacing,show_raw_data,varargin)


if ~isempty(varargin) && ~isempty(varargin{1})
    do_std_map = varargin{1};
else
    do_std_map = 0;
end

if length(varargin) > 1 && ~isempty(varargin{2})
    do_corr_map = varargin{2};
else
    do_corr_map = 0;
end

if length(varargin) > 2 && ~isempty(varargin{3})
    these_colors = varargin{3};
else
    for i = 1:length(traces)        
        these_colors{i} = zeros(size(traces{i},1),3);
    end
end

assignin('base','these_colors',these_colors)

[maps, mpp_maps, map_index, color_maps] = build_slm_maps_multi(...
    traces,mpp,sequence,stim_key,spacing,[-150 150],[-150 150],these_colors);

assignin('base','color_maps',color_maps)

loc_names = cell(size(maps{1}));
center = ceil((size(maps{1})-1)*spacing/2) + 1;

% for i = 1:size(maps{1},1)
%     for j = 1:size(maps{1},2)
%         
%         loc_names{i,j} = [num2str((i-1)*spacing - center(1)) ', ' ...
%                          num2str((j-1)*spacing - center(2)) ' um'];
%                      
%     end
% end

% assignin('base','loc_names',loc_names)

trials_per_stack = 15;
if show_raw_data
    figure
    subplot(121)

    plot_trace_stack_grid(maps{1},trials_per_stack,1,0,[],[],[],mpp_maps{1});%,loc_names);
%     title(['Power = ' num2str(sequence(1).target_power) ' mW'])
    subplot(122)
%     figure
	plot_trace_stack_grid(maps{2},trials_per_stack,1,0,[],[],[],mpp_maps{2},loc_names);


end


if do_std_map
    figure
    stddev_maps{1} = get_stdev_map(maps{1},1,0);
%     stddev_maps{2} = get_stdev_map(maps{2},0,0);
else
    stddev_maps = cell(2,1);
end

 
if do_corr_map
    corr_maps{1} = get_corr_image(maps{1},0,0);
    corr_maps{2} = get_corr_image(maps{2},0,0);
else
    corr_maps = cell(2,1);
end
% 
% if do_corr_map
%     
%     figure
%     
%     subplot(121); 
%     imagesc(corr_maps{1}); caxis([0 0.5])
%     title(['Cell 1 Corr Map: Power = ' num2str(sequence(1).target_power) ' mW'])
% 
%     subplot(122); 
%     imagesc(stddev_maps{1}); %caxis([0 1])
%     title(['Cell 1 Stddev Map: Power = ' num2str(sequence(1).target_power) ' mW'])
%     
% end
% if do_std_map
%     
%     figure 
%     
%     subplot(121);
%     imagesc(corr_maps{2}); caxis([0 0.5])
%     title(['Cell 2 Corr Map: Power = ' num2str(sequence(1).target_power) ' mW'])
% 
%     subplot(122); 
%     imagesc(stddev_maps{2}); %caxis([0 1])
%     title(['Cell 2 Stddev Map: Power = ' num2str(sequence(1).target_power) ' mW'])
% 
% 
% end

% maps = maps;
corr_maps = corr_maps{1};
mpp_maps = mpp_maps{1};
% assignin('base','maps',maps)
