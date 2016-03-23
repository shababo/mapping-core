

norm_current_image_3_4_s1c2_r7 = current_image_3_4_s1c2_r7/max(max(current_image_3_4_s1c2_r7));
norm_current_image_3_4_s2c3_r3 = current_image_3_4_s2c3_r3/max(max(current_image_3_4_s2c3_r3));
norm_current_image_3_4_s3c1_r3 = current_image_3_4_s3c1_r3/max(max(current_image_3_4_s3c1_r3));
norm_current_image_3_5_s3c1_r3 = current_image_3_5_s3c1_r3/max(max(current_image_3_5_s3c1_r3));
norm_current_image_3_7_s2c2_r3 = current_image_3_7_s2c2_r3/max(max(current_image_3_7_s2c2_r3));
norm_current_image_3_7_s4c3_r3 = current_image_3_7_s4c3_r3/max(max(current_image_3_7_s4c3_r3));

mean_soma_map = (current_image_3_4_s1c2_r7 + ...
                current_image_3_4_s2c3_r3 + ...
                current_image_3_4_s3c1_r3 + ...
                current_image_3_5_s3c1_r3 + ...
                current_image_3_7_s2c2_r3 + ...
                current_image_3_7_s4c3_r3)/6;
            
mean_norm_soma_map = (norm_current_image_3_4_s1c2_r7 + ...
    norm_current_image_3_4_s2c3_r3 + ...
    norm_current_image_3_4_s3c1_r3 + ...
    norm_current_image_3_5_s3c1_r3 + ...
    norm_current_image_3_7_s2c2_r3 + ...
    norm_current_image_3_7_s4c3_r3)/6;

figure;
imagesc(mean_soma_map)
colormap hot
colorbar
axis off
title('Average soma-ChR2 Current Map (N = 6)')

figure;
imagesc(mean_norm_soma_map)
colormap hot
colorbar
caxis([0 1])
axis off
title('Average Normalized soma-ChR2 Current Map (N = 6)')



norm_current_image_3_8_s3c2_r3 = current_image_3_8_s3c2_r3/max(max(current_image_3_8_s3c2_r3));
norm_current_image_3_8_s4c1_r4 = current_image_3_8_s4c1_r4/max(max(current_image_3_8_s4c1_r4));
norm_current_image_3_8_s4c2_r4 = current_image_3_8_s4c2_r4/max(max(current_image_3_8_s4c2_r4));
norm_current_image_3_8_s5c1_r3 = current_image_3_8_s5c1_r3/max(max(current_image_3_8_s5c1_r3));
norm_current_image_3_8_s5c2_r1 = current_image_3_8_s5c2_r1/max(max(current_image_3_8_s5c2_r1));
norm_current_image_3_3_s4c1_r3 = current_image_3_3_s4c1_r3/max(max(current_image_3_3_s4c1_r3));


mean_chrimson_map = (current_image_3_8_s3c2_r3 + ...
                current_image_3_8_s4c1_r4 + ...
                current_image_3_8_s4c2_r4 + ...
                current_image_3_8_s5c1_r3 + ...
                current_image_3_3_s4c1_r3 + ...
                current_image_3_8_s5c2_r1)/6;
            
mean_norm_chrimson_map = (norm_current_image_3_8_s3c2_r3 + ...
    norm_current_image_3_8_s4c1_r4 + ...
    norm_current_image_3_8_s4c2_r4 + ...
    norm_current_image_3_8_s5c1_r3 + ...
    norm_current_image_3_3_s4c1_r3 + ...
    norm_current_image_3_8_s5c2_r1)/6;

figure;
imagesc(mean_chrimson_map)
colormap hot
colorbar
axis off
title('Average ChrimsonR Current Map (N = 6)')

figure;
imagesc(mean_norm_chrimson_map)
colormap hot
colorbar
caxis([0 1])
axis off
title('Average Normalized ChrimsonR Current Map (N = 6)')

%%

soma_maps = zeros(6,21,21);
soma_maps(1,:,:) = current_image_3_4_s1c2_r7;
soma_maps(2,:,:) = current_image_3_4_s2c3_r3;
soma_maps(3,:,:) = current_image_3_4_s3c1_r3;
soma_maps(4,:,:) = current_image_3_5_s3c1_r3;
soma_maps(5,:,:) = current_image_3_7_s2c2_r3;
soma_maps(6,:,:) = current_image_3_7_s4c3_r3;

chrimson_maps = zeros(6,21,21);
chrimson_maps(1,:,:) = current_imgiage_3_8_s3c2_r3;
chrimson_maps(2,:,:) = current_image_3_8_s4c1_r4;
chrimson_maps(3,:,:) = current_image_3_8_s4c2_r4;
chrimson_maps(4,:,:) = current_image_3_8_s5c1_r3;
chrimson_maps(5,:,:) = current_image_3_3_s4c1_r3;
chrimson_maps(6,:,:) = current_image_3_8_s5c2_r1;

figure
for i = 1:size(soma_maps,1)
    subplot(2,3,i)
    imagesc(squeeze(soma_maps(i,:,:)))
    colormap hot
    colorbar
    axis off
end

figure
for i = 1:size(chrimson_maps,1)
    subplot(2,3,i)
    imagesc(squeeze(chrimson_maps(i,:,:)))
    colormap hot
    colorbar
    axis off
end
                
            
%%    
  
figure; 
errorbar(peak_currents_cells_by_trial_mean_norm(:,1:switch_ind)',peak_currents_cells_by_trial_std_norm(:,1:switch_ind)');
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            