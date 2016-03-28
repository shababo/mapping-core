

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


soma_maps = zeros(21,21,6);
soma_maps(:,:,1) = current_image_3_4_s1c2_r7;
soma_maps(:,:,2) = current_image_3_4_s2c3_r3;
soma_maps(:,:,3) = current_image_3_4_s3c1_r3;
soma_maps(:,:,4) = current_image_3_5_s3c1_r3;
soma_maps(:,:,5) = current_image_3_7_s2c2_r3;
soma_maps(:,:,6) = current_image_3_7_s4c3_r3;

chrimson_maps = zeros(21,21,6);
chrimson_maps(:,:,1) = current_image_3_8_s3c2_r3;
chrimson_maps(:,:,2) = current_image_3_8_s4c1_r4;
chrimson_maps(:,:,3) = current_image_3_8_s4c2_r4;
chrimson_maps(:,:,4) = current_image_3_8_s5c1_r3;
chrimson_maps(:,:,5) = current_image_3_3_s4c1_r3;
chrimson_maps(:,:,6) = current_image_3_8_s5c2_r1;

norm_soma_maps = zeros(21,21,6);
norm_soma_maps(:,:,1) = norm_current_image_3_4_s1c2_r7;
norm_soma_maps(:,:,2) = norm_current_image_3_4_s2c3_r3;
norm_soma_maps(:,:,3) = norm_current_image_3_4_s3c1_r3;
norm_soma_maps(:,:,4) = norm_current_image_3_5_s3c1_r3;
norm_soma_maps(:,:,5) = norm_current_image_3_7_s2c2_r3;
norm_soma_maps(:,:,6) = norm_current_image_3_7_s4c3_r3;

norm_chrimson_maps = zeros(21,21,6);
norm_chrimson_maps(:,:,1) = norm_current_image_3_8_s3c2_r3;
norm_chrimson_maps(:,:,2) = norm_current_image_3_8_s4c1_r4;
norm_chrimson_maps(:,:,3) = norm_current_image_3_8_s4c2_r4;
norm_chrimson_maps(:,:,4) = norm_current_image_3_8_s5c1_r3;
norm_chrimson_maps(:,:,5) = norm_current_image_3_3_s4c1_r3;
norm_chrimson_maps(:,:,6) = norm_current_image_3_8_s5c2_r1;


%%

figure
for i = 1:size(soma_maps,3)
    subplot(2,3,i)
    imagesc(soma_maps(:,:,i))
    colormap hot
    colorbar
    axis off
end

figure
for i = 1:size(chrimson_maps,3)
    subplot(2,3,i)
    imagesc(chrimson_maps(:,:,i))
    colormap hot
    colorbar
    axis off
end
                

%%

mean_soma = mean(soma_maps,3);
mean_chrimson = mean(chrimson_maps,3);
std_soma = std(soma_maps,[],3);
std_chrimson = std(chrimson_maps,[],3);

norm_mean_soma = mean(norm_soma_maps,3);
norm_mean_chrimson = mean(norm_chrimson_maps,3);
norm_std_soma = std(norm_soma_maps,[],3);
norm_std_chrimson = std(norm_chrimson_maps,[],3);
%%
figure;          
subplot(2,2,1)
imagesc(mean_soma)
subplot(2,2,2)
imagesc(std_soma)
subplot(2,2,3)
imagesc(mean_chrimson)
subplot(2,2,4)
imagesc(std_chrimson)
            
%%    
  
figure; 
errorbar(-100:10:100,mean_chrimson(11,:)',std_chrimson(11,:)','r');
hold on
errorbar(-100:10:100,mean_soma(11,:)',std_soma(11,:)','b');
xlim([-100 100])
%%

figure; 
errorbar(-100:10:100,norm_mean_chrimson(11,:)',norm_std_chrimson(11,:)');
hold on
errorbar(-100:10:100,norm_mean_soma(11,:)',norm_std_soma(11,:)');

%%

gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
startPoints = [1 0 100 .25];
            
f1 = fit((-100:10:100)',norm_mean_chrimson(11,:)',gaussEqn,'Start', startPoints)
f2 = fit((-100:10:100)',norm_mean_soma(11,:)',gaussEqn,'Start', startPoints)
figure
plot(f1,(-100:10:100)',norm_mean_chrimson(11,:)')
hold on
plot(f2,(-100:10:100)',norm_mean_soma(11,:)')

%%

figure; 
errorbar(-100:10:100,norm_mean_chrimson(11,:)',norm_std_chrimson(11,:)','.r');
hold on
errorbar(-100:10:100,norm_mean_soma(11,:)',norm_std_soma(11,:)','.b');
hold on
plot(f1,'r-')
hold on
plot(f2,'b-')
xlim([-100 100])


            

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            