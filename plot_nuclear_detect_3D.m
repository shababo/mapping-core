function detect_img = plot_nuclear_detect_3D(stack_filename,centers,varargin)

info = imfinfo(stack_filename);
if ~isempty(varargin) && ~isempty(varargin{1})
    colors = varargin{1};
else
%     colors = flipud(jet(100));
colors = zeros(size(centers,2),3);
%     color_inds = ceil(centers(3,:)/max(centers(3,:))*100);
end
if length(varargin) > 1 && ~isempty(varargin{2})
    sizes = varargin{2};
else
    sizes = 5*ones(size(centers,2),1);
end

sizes = 5*ones(size(centers,2),1);
centers(:,sizes <= 0) = [];
sizes(sizes <= 0) = [];



% h = figure;
% set(h,'position',[724 111 948 815])
plotted_cells = [];
M = 12500;
G = linspace(0,1,M)';
myGmap = horzcat(G, zeros(size(G)) , zeros(size(G)));
for i = 1:length(info)
    
    tmp=imread(stack_filename,i);
%     imagesc(tmp);
%     caxis([0 4000])
%     colormap(myGmap)
%     axis image;
%     hold on;
    if i == 1
        max_proj = tmp;%zeros(size(tmp));
    else
        max_proj(tmp(:) > max_proj(:)) = tmp(tmp(:) > max_proj(:));
    end
        
    
    
    these_cells = find(abs(centers(3,:) - i) < 7);
    plotted_cells = union(plotted_cells,these_cells);
%     scatter(centers(1,these_cells), centers(2,these_cells),sizes(these_cells),colors(these_cells,:)+1,'filled')
%     hold on
% %     scatter(5, 200,sizes(1),[0 0 0],'filled')
%     tmp2=getframe;
%     image(tmp2.cdata);
%     if i == 1
%         imwrite(tmp2.cdata,'tmp.tif');
%     else
%         imwrite(tmp2.cdata,'tmp.tif', 'writemode', 'append');
%     end
%     hold off
    
end
% return
% close(h)
h = figure;
% max_proj(:) = 0;
imagesc(max_proj);
caxis([0 4000])
colormap(myGmap)
axis image;
axis off
hold on;


% these_cells = centers;
scatter(centers(1,plotted_cells), centers(2,plotted_cells),sizes(plotted_cells),...
    colors(plotted_cells,:)+1,'filled')
hold on
% scatter(5, 200,sizes(1)+10,[0 0 0],'filled')
tmp2=getframe;
detect_img = tmp2.cdata;
% xlim([145 200]);
% ylim([110 165])
% close(h)

    
    