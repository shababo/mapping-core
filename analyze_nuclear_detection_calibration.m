% function analyze_nuclear_detection_calibration(images,nuclear_locs_image_coord,fluor_vals)
% 
% 
% 
% figure 
% linespecs ={'b.','ro'};
% for i = 1:length(nuclear_locs_image_coord)
%     subplot(121)
%     histogram(fluor_vals{i})
% %     histogram(fluor_vals{i}/max(fluor_vals{i}),(0:.1:30)/30)
%     hold on
%     subplot(122)
%     scatter3(nuclear_locs_image_coord{i}(:,1),nuclear_locs_image_coord{i}(:,2),...
%         nuclear_locs_image_coord{i}(:,3),linespecs{i})
%     hold on
% end

%%

median_size = [1 2 5];
radius_bg = [25 50 100];
blur_size = [2 5 10];
radius_dilation = [1 2 5];
min_voxel_remove = [3 7 15];
curvature_distance = [2 4 8];
curvature_numvoxels = [2 5 9];
sigma_default = [7.5 7.5 15]; 
thrdist = [1 3 8];
removefp_mindist = [3 6 10];

origfile = ['/media/shababo/data/nuc_detect_calibration/pockels_2_200_C0.tif'];
basefile = ['/media/shababo/data/nuc_detect_calibration/pockels_2_200_C0_%s.tif'];

exp_count = 1;
for i1 = 1%:length(median_size)
    for i2 = 2%1:length(radius_bg)
        for i3 = 1%:length(blur_size)
            for i4 = 2%1:length(radius_dilation)
                for i5 = 1:length(min_voxel_remove)
                    for i6 = 1:length(curvature_distance)
                        for i7 = 1:length(curvature_numvoxels)
%                             for i8= 1:length(sigma_default)
                                for i9 = 3%1:length(thrdist)
                                    for i10 = 2%1:length(removefp_mindist)
                                        
                                        tic
                                        params.prefilter_option = sprintf('x=%d y=%d z=%d',median_size(i1)*ones(1,3));
                                        params.radius_background = radius_bg(i2);
                                        params.blur_option =sprintf('x=%d y=%d z=%d',blur_size(i3)*ones(1,3));
                                        params.radius_dilation = [1,1,1]*radius_dilation(i4);
                                        params.min_voxel_remove = min_voxel_remove(i5);
                                        params.curvature_distance = [1,1,1]*curvature_distance(i6);
                                        params.curvature_numvoxels = curvature_numvoxels(i7);
                                        params.sigma_default = sigma_default;
                                        params.thrdist = thrdist(i9);
                                        params.removefp_min_dist = removefp_mindist(i10);
                                        
                                        exp_id = ['experiment_negcurve_' num2str(exp_count)];
                                        save([exp_id '_params.mat'],'params')
                                        newfilename = sprintf(basefile,exp_id);
                                        copyfile(origfile,newfilename)
                                        [nuclear_locs,fluor_vals,nuclear_locs_image_coord] = detect_nuclei(newfilename(1:end-4),...
                                            1.89,[145.924; 147.597],...
                                            2,1,0,0,params);
                                        exp_time = toc;
                                        save([exp_id '_result.mat'],'params',...
                                            'nuclear_locs','fluor_vals','nuclear_locs_image_coord','exp_time')
                                        exp_count = exp_count + 1;
                                        
                                    end
                                end
%                             end
                        end
                    end
                end
            end
        end
    end
end


