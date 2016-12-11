%Load your images

image_tags = {'100816_slice3','101216_slice4','101216_slice2','100716_slice3',...
    '103116_slice4','103116_slice5','103116_slice6','110116_slice1_L5p',...
    '110116_slice3','110116_slice4','110116_slice5_L45p','110116_slice6',...
    '110116_slice7','110116_slice8','110216_slice1','110216_slice2','110216_slice3',...
    '110216_slice4','101216_slice3'};

image_pair_type = [0 0 0 0 1 1 1 1 0 0 1 0 0 0 1 1 1 1 0];


%%

for i = 1:length(image_tags)
    
    i
    expname=image_tags{i};%'100816_slice3'
    dirname = '/media/shababo/data/aligned-maps/';
    layers_filename=[dirname,expname,'_alignment_layers_fix-01.jpg']
    try
        
        layers_image = imread(layers_filename);
    catch e
        continue
    end

    ch1_filename=[dirname,expname,'_alignment_q02_ch1-01.jpg'];
    % x=imfinfo(dapi_filename);
    ch1_image = imread(ch1_filename);

    ch2_filename=[dirname,expname,'_alignment_q02_ch2-01.jpg'];
    ch2_image = imread(ch2_filename);


    layers_image = squeeze(layers_image(:,:,1));
    ch1_image = squeeze(ch1_image(:,:,1));
    ch2_image = squeeze(ch2_image(:,:,1));
    % rotate images

    merge_image=[]
    merge_image(:,:,3)=double(layers_image)/max(double(layers_image(:)));
    merge_image(:,:,2)=double(ch2_image)/max(double(ch2_image(:)));
% 
% 
%     % global rtHandles
    moveon = 0;
    rotateToolNew(merge_image);
    while ~moveon
        drawnow
        pause(.5)
    end
    
    % rtHandles.figure1.UserData=merge_image

    %

    %     rot_rec=[];
    %     rot_ch1=[];
    %     rot_tdt=[];
        disp(theta)
        theta_all(i) = theta
        layers_image([1,2,3,4,5,end-4,end-3,end-2,end-1,end],:) = 0;
        ch1_image([1,2,3,4,5,6,7,8,end-4,end-3,end-2,end-1,end],:) = 0;
        ch2_image([1,2,3,4,5,6,7,8,end-4,end-3,end-2,end-1,end],:) = 0;
        layers_image = double(layers_image) + eps;
        ch1_image = double(ch1_image) + eps;
        ch2_image = double(ch2_image) + eps;
        rot_layers=imrotate(layers_image,theta_all(i),'bilinear','loose');
        rot_ch1=imrotate(ch1_image,theta_all(i),'bilinear','loose');
        rot_ch2=imrotate(ch2_image,theta_all(i),'bilinear','loose');
        rot_layers(rot_layers == 0) = NaN;
        rot_ch1(rot_ch1 == 0) = NaN;
        rot_ch2(rot_ch2 == 0) = NaN;
        rot_layers = double(rot_layers) - eps;
        rot_ch1 = double(rot_ch1) - eps;
        rot_ch2 = double(rot_ch2) - eps;
        ch1_maps{i} = double(rot_ch1 > 20);
        ch1_maps{i}(isnan(rot_ch1)) = NaN;
        ch2_maps{i} = double(rot_ch2 > 20);
        ch2_maps{i}(isnan(rot_ch2)) = NaN;
    %     rot_merge=[];
    %     rot_merge(:,:,3)=double(rot_dapi)/max(double(rot_dapi(:)));
    %     rot_merge(:,:,1)=double(rot_rec_max)/max(double(rot_rec_max(:)));




        % Label layer boundaries


        figure
        bigfig
        colormap gray
        hando(1)=subplot(1,2,1)
        imagesc(rot_layers)
        axis image
        hando(3)=subplot(1,2,2)
        hold on
        set(gca,'Ydir','reverse')
        plot(sydapinorm,1:length(sydapinorm),'b')
        plot(ydapinorm,1:length(ydapinorm),'r')
        ylim([0 size(crot_dapi,1)])
        set(hando(3),'PlotBoxAspectRatio',[size(crot_rec_max,2) size(crot_rec_max,1) 1]);
        set(hando(3),'Ytick',get(hando(1),'Ytick'))
        title(hando(3),'Then hit enter')

    %     title(hando(1),'Pia')
    %     temp =imline(hando(1), [size(crot_rec_max,2)*0.2 size(crot_rec_max,2)*0.8], [size(crot_rec_max,1)/2 size(crot_rec_max,1)/2]);
    %     pause
    %     pia=getPosition(temp)
    %     pia=round(mean(pia(:,2)));
    %     
    %     title(hando(1),'L1-2/3')
    %     temp =imline(hando(1), [size(crot_rec_max,2)*0.2 size(crot_rec_max,2)*0.8], [size(crot_rec_max,1)/2 size(crot_rec_max,1)/2]);
    %     pause
    %     L1L23border=getPosition(temp)
    %     L1L23=round(mean(L1L23border(:,2)));
    %     
    %     title(hando(1),'L2/3-L4')
    %     temp =imline(hando(1), [size(crot_rec_max,2)*0.2 size(crot_rec_max,2)*0.8], [size(crot_rec_max,1)/2 size(crot_rec_max,1)/2]);
    %     pause
    %     L23L4border=getPosition(temp)
    %     L23L4=round(mean(L23L4border(:,2)));

        title(hando(1),'L4-L5')
        temp =imline(hando(1), [size(crot_rec_max,2)*0.2 size(crot_rec_max,2)*0.8], [size(crot_rec_max,1)/2 size(crot_rec_max,1)/2]);
        pause
        L4L5border=getPosition(temp)
        L4L5(i)=round(mean(L4L5border(:,2)));

        title(hando(1),'L5-L6')
        temp =imline(hando(1), [size(crot_rec_max,2)*0.2 size(crot_rec_max,2)*0.8], [size(crot_rec_max,1)/2 size(crot_rec_max,1)/2]);
        pause
        L5L6border=getPosition(temp)
        L5L6(i)=round(mean(L5L6border(:,2)));

        % title(hando(1),'White Matter')
        % temp =imline(hando(1), [size(crot_rec_max,2)*0.2 size(crot_rec_max,2)*0.8], [size(crot_rec_max,1)/2 size(crot_rec_max,1)/2]);
        % pause
        % WM=getPosition(temp)
        % WM=round(mean(WM(:,2)));

        close(gcf)

end

%%
l45_border = max(L4L5);

avg_map_base = nan(max(cellfun(@(x,y) l45_border - L4L5(i) + size(x,1),ch1_maps,mat2cell(L4L5,1,ones(1,19)))),...
    max(cellfun(@(x) size(x,2),ch1_maps)));

% avg_map_l4 = zeros(size(avg_map_l5));
center_col = ceil(size(avg_map_base,2)/2);
l5_map_count = 1;
l4_map_count = 1;

aligned_l5_maps = [];
aligned_l4_maps = [];
for i = 1:length(image_tags)
    
    start_row = l45_border - L4L5(i) + 1;
    end_row = start_row + size(ch1_maps{i},1) - 1;
    start_col = center_col - ceil(size(ch1_maps{i},2)/2) + 1;
    end_col = start_col + size(ch1_maps{i},2) - 1;
    
    aligned_l5_maps(:,:,l5_map_count) = avg_map_base;
    aligned_l5_maps(start_row:end_row,start_col:end_col,l5_map_count) = ch1_maps{i};
    l5_map_count = l5_map_count + 1;
    
    if image_pair_type(i)
        aligned_l5_maps(:,:,l5_map_count) = avg_map_base;
        aligned_l5_maps(start_row:end_row,start_col:end_col,l5_map_count) = ch2_maps{i};
        l5_map_count = l5_map_count + 1; 
    else
        aligned_l4_maps(:,:,l4_map_count) = avg_map_base;
        aligned_l4_maps(start_row:end_row,start_col:end_col,l4_map_count) = ch2_maps{i};
        l4_map_count = l4_map_count + 1;
    end
    
end
    
%%

for i = 1:size(aligned_l4_maps,3)
    i
    aligned_l4_maps_smooth(:,:,i) = imgaussfilt(aligned_l4_maps(:,:,i),5*2.5);
end

for i = 1:size(aligned_l5_maps,3)
    i
    aligned_l5_maps_smooth(:,:,i) = imgaussfilt(aligned_l5_maps(:,:,i),5*2.5);
end

%%
figure; 
subplot(121)
imagesc(mean(~isnan(aligned_l4_maps),3) > .5)
hold on
plot([1 size(aligned_l4_maps,2)],[l45_border l45_border],'w')
subplot(122)
imagesc(mean(~isnan(aligned_l5_maps),3) > .5)%(:,:,randsample(1:size(aligned_l5_maps,3),6))
hold on
plot([1 size(aligned_l4_maps,2)],[l45_border l45_border],'w')

%%
occupancy = mean(~isnan(aligned_l4_maps),3) > .8;
xlimits = [min(find(occupancy(1000,:)))+50 max(find(occupancy(1000,:)))-50];
ylimits = [min(find(occupancy(:,600))) max(find(occupancy(:,600)))-50];


%%

aligned_l4_maps_trunc = aligned_l4_maps_smooth(ylimits(1):ylimits(2),xlimits(1):xlimits(2),:);
aligned_l5_maps_trunc = aligned_l5_maps_smooth(ylimits(1):ylimits(2),xlimits(1):xlimits(2),:);

figure; 
subplot(121)
imagesc(nanmean(aligned_l4_maps_trunc,3))
hold on
plot([1 size(aligned_l4_maps_trunc,2)],[l45_border l45_border]-ylimits(1),'w')
hold on
% plot([1 size(aligned_l4_maps_trunc,2)],[l45_border l45_border]-ylimits(1)+120*2.5,'w--')
hold on
% plot([1 size(aligned_l4_maps_trunc,2)],[l45_border l45_border]-ylimits(1)+300*2.5,'w')
hold on
plot(center_col + -100*2.5*[1 1]-xlimits(1),[1 size(aligned_l4_maps_trunc,1)],'w--')
hold on
plot(center_col + 100*2.5*[1 1]-xlimits(1),[1 size(aligned_l4_maps_trunc,1)],'w--')
axis image
caxis([0 .75])
% xlim(xlimits)
% ylim(ylimits)
axis off
subplot(122)
imagesc(nanmean(aligned_l5_maps_trunc,3))%(:,:,randsample(1:size(aligned_l5_maps,3),6))
hold on
% plot([1 size(aligned_l4_maps_trunc,2)],[l45_border l45_border]-ylimits(1)+120*2.5,'w--')
hold on
plot([1 size(aligned_l4_maps_trunc,2)],[l45_border l45_border]-ylimits(1),'w')
hold on
% plot([1 size(aligned_l4_maps_trunc,2)],[l45_border l45_border]-ylimits(1)+300*2.5,'w')
hold on
plot(center_col + -100*2.5*[1 1]-xlimits(1),[1 size(aligned_l4_maps_trunc,1)],'w--')
hold on
plot(center_col + 100*2.5*[1 1]-xlimits(1),[1 size(aligned_l4_maps_trunc,1)],'w--')
axis image
% xlim(xlimits)
% ylim(ylimits)
axis off
caxis([0 .75])
%%
window_size = 800;
smooth_std = 300;

figure; 
% subplot(211)
x = (1:size(aligned_l4_maps_trunc,2))';
% mean_vec = smoothts(squeeze(nanmean(nanmean(aligned_l4_maps_trunc,1),3)),'g',window_size,smooth_std)';
mean_vec = squeeze(nanmean(nanmean(aligned_l4_maps_trunc,1),3))';
% error = smoothts(squeeze(nanstd(nanmean(aligned_l4_maps_trunc,1),[],3)),'g',window_size,smooth_std)'/sqrt(size(aligned_l4_maps_trunc,3));  
error = squeeze(nanstd(nanmean(aligned_l4_maps_trunc,1),[],3))'/sqrt(size(aligned_l4_maps_trunc,3));
fill([x;flipud(x)],[(mean_vec-error);flipud(mean_vec+error)]/max(mean_vec),[246 146 31]*.25/255,'linestyle','none'); alpha(.5)
hold on
plot(x,mean_vec/max(mean_vec),'color',[246 146 31]/255)
hold on
plot(center_col + -100*2.5*[1 1]-xlimits(1),[0 1],'k--')
hold on
plot(center_col + 100*2.5*[1 1]-xlimits(1),[0 1],'k--')
axis off
figure
% subplot(212)
% mean_vec = smoothts(squeeze(nanmean(nanmean(aligned_l5_maps_trunc,1),3)),'g',window_size,smooth_std)';
mean_vec = squeeze(nanmean(nanmean(aligned_l5_maps_trunc,1),3))';
error = squeeze(nanstd(nanmean(aligned_l5_maps_trunc,1),[],3))'/sqrt(size(aligned_l5_maps_trunc,3));  
fill([x;flipud(x)],[(mean_vec-error);flipud(mean_vec+error)]/max(mean_vec),[39 168 224]*.25/255,'linestyle','none'); alpha(.5)
hold on
plot(x,mean_vec/max(mean_vec),'Color',[39 168 224]/255)
hold on
plot(center_col + -100*2.5*[1 1]-xlimits(1),[0 1],'k--')
hold on
plot(center_col + 100*2.5*[1 1]-xlimits(1),[0 1],'k--')
axis off
% xlim(xlimits)
figure
% subplot(121)
x = flipud((1:size(aligned_l4_maps_trunc,1))');
% mean_vec = smoothts(squeeze(nanmean(nanmean(aligned_l4_maps_trunc,2),3))','g',window_size,smooth_std)';
mean_vec = squeeze(nanmean(nanmean(aligned_l4_maps_trunc,2),3));
error = squeeze(nanstd(nanmean(aligned_l4_maps_trunc,2),[],3))/sqrt(size(aligned_l4_maps_trunc,3)); 
fill([x;flipud(x)],[(mean_vec-error);flipud(mean_vec+error)]/max(mean_vec),[246 146 31]/255*.25,'linestyle','none'); alpha(.5)
hold on
plot(x,mean_vec/max(mean_vec),'color',[246 146 31]/255)
hold on
hold on
plot(size(aligned_l4_maps_trunc,1) - [l45_border l45_border]+ylimits(1),[0 1],'k')
hold on
% plot(size(aligned_l4_maps_trunc,1) - [l45_border l45_border]-ylimits(1),[0 .1],'k')
axis off
view(-90,90) 
% xlim(size(aligned_l4_maps,1) - fliplr(ylimits))
figure
% subplot(122)
% mean_vec = smoothts(squeeze(nanmean(nanmean(aligned_l5_maps_trunc,2),3))','g',window_size,smooth_std)';
mean_vec = squeeze(nanmean(nanmean(aligned_l5_maps_trunc,2),3));
error = squeeze(nanstd(nanmean(aligned_l5_maps_trunc,2),[],3))/sqrt(size(aligned_l5_maps_trunc,3));  
fill([x;flipud(x)],[(mean_vec-error);flipud(mean_vec+error)]/max(mean_vec),[39 168 224]*.25/255,'linestyle','none'); alpha(.5)
hold on
plot(x,mean_vec/max(mean_vec),'color',[39 168 224]/255)
hold on
plot(size(aligned_l4_maps_trunc,1) - [l45_border l45_border]+ylimits(1),[0 1],'k')
hold on
% plot(size(aligned_l4_maps_trunc,1) - ([l45_border l45_border] + 300*2.5-ylimits(1)),[0 1],'k')
hold on
% plot(size(aligned_l4_maps_trunc,1) - ([l45_border l45_border] + 125*2.5-ylimits(1)),[0 1],'k--')
axis off
view(-90,90) 
% xlim(size(aligned_l4_maps_trunc,1) - fliplr(ylimits))
    