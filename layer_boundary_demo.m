%Load your images

clear global
clearvars -except SOMmorph


expname='161008_slice5'
dirname = '/media/shababo/data/som-map-imaging/projections/';
recon_filename=[dirname,expname,'_647.tif']
recon_image = imread(recon_filename);

dapi_filename=[dirname,expname,'_DAPI.tif']
x=imfinfo(dapi_filename);
dapi_image = imread(dapi_filename);

tdt_filename=[dirname,expname,'_TdT.tif']
tdt_image = imread(tdt_filename);

%% make sure all are the right size/resolution


imsize=max([length(recon_image) length(dapi_image) length(tdt_image)])

tdt_image=resizem(tdt_image,[imsize imsize]);
recon_image=resizem(recon_image,[imsize imsize]);
dapi_image=resizem(dapi_image,[imsize imsize]);

%% rotate images

merge_image=[]
merge_image(:,:,3)=double(dapi_image)/max(double(dapi_image(:)));
merge_image(:,:,2)=double(recon_image)/max(double(recon_image(:)));


% global rtHandles
rotateToolNew(merge_image);
% rtHandles.figure1.UserData=merge_image

%% 

    rot_rec=[];
    rot_dapi=[];
    rot_tdt=[];

    rot_rec_max=imrotate(recon_image,theta,'bilinear','loose');
    rot_dapi=imrotate(dapi_image,theta,'bilinear','loose');
    rot_tdt=imrotate(tdt_image,theta,'bilinear','loose');

    rot_merge=[];
    rot_merge(:,:,3)=double(rot_dapi)/max(double(rot_dapi(:)));
    rot_merge(:,:,1)=double(rot_rec_max)/max(double(rot_rec_max(:)));

    %% Crop images

  
     
    figure
    bigfig
    imagesc(rot_merge)
    axis image
    [x,cropbord]=imcrop;
    imagesc(x)
    axis image
    
    
    crot_dapi=rot_dapi(floor(cropbord(2)):floor(cropbord(2))+floor(cropbord(4)),...
        floor(cropbord(1)):floor(cropbord(1))+floor(cropbord(3)));
    crot_rec_max=rot_rec_max(floor(cropbord(2)):floor(cropbord(2))+floor(cropbord(4)),...
        floor(cropbord(1)):floor(cropbord(1))+floor(cropbord(3)));
    crot_tdt=rot_tdt(floor(cropbord(2)):floor(cropbord(2))+floor(cropbord(4)),...
        floor(cropbord(1)):floor(cropbord(1))+floor(cropbord(3)));
    
% %     optional: rescale axes to microns
%     rsfactor=720/imsize %for 20x images
%     crot_dapi=resizem(crot_dapi,rsfactor)
%     crot_rec_max=resizem(crot_rec_max,rsfactor)
%     crot_tdt=resizem(crot_tdt,rsfactor)

% %     optional: rescale axes to 0.5 micron units
    rsfactor=720/imsize*2 %for 20x images
    crot_dapi=resizem(crot_dapi,rsfactor);
    crot_rec_max=resizem(crot_rec_max,rsfactor);
    crot_tdt=resizem(crot_tdt,rsfactor);

    
    crot_dapi(ismembertol(crot_dapi,0))=nan;
    crot_rec_max(ismembertol(crot_rec_max,0))=nan;
    crot_tdt(ismembertol(crot_tdt,0))=nan;
    
    ydapi=nanmean(crot_dapi,2);
    yrec_max=nanmean(crot_rec_max,2);
    ydapinorm=ydapi/max(ydapi);
    yrecnorm_max=yrec_max/max(yrec_max);
    
    sydapi=smooth(ydapi,100);
    syrec_max=smooth(yrec_max,100);
    
    sydapinorm=sydapi/max(sydapi);
    syrec_maxnorm=syrec_max/max(syrec_max);

    
    
    %% Label layer boundaries
    
    
    figure
    bigfig
    colormap gray
    hando(1)=subplot(1,2,1)
    imagesc(crot_dapi)
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
    L4L5=round(mean(L4L5border(:,2)));
    
    title(hando(1),'L5-L6')
    temp =imline(hando(1), [size(crot_rec_max,2)*0.2 size(crot_rec_max,2)*0.8], [size(crot_rec_max,1)/2 size(crot_rec_max,1)/2]);
    pause
    L5L6border=getPosition(temp)
    L5L6=round(mean(L5L6border(:,2)));
    
    % title(hando(1),'White Matter')
    % temp =imline(hando(1), [size(crot_rec_max,2)*0.2 size(crot_rec_max,2)*0.8], [size(crot_rec_max,1)/2 size(crot_rec_max,1)/2]);
    % pause
    % WM=getPosition(temp)
    % WM=round(mean(WM(:,2)));
   
    close(gcf)
    
    %% 
    
    figure
    bigfig

    imagesc(crot_rec_max)
    axis image
    title('Cell1')
    soma1pos=ginput(1)
    title('Cell2')
    soma2pos=ginput(1)
