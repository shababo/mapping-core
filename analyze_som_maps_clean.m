detection_grid_5_13_s2c2_r4_2000 = cell(size(posteriors_grid_5_13_s2c2_r4_2000));



all_points = [];
for i = 1:size(posteriors_grid_5_13_s2c2_r4_2000,1)
    for j = 1:size(posteriors_grid_5_13_s2c2_r4_2000,2)
        for k = 1:size(posteriors_grid_5_13_s2c2_r4_2000{i,j},1)
            all_points = [all_points posteriors_grid_5_13_s2c2_r4_2000{i,j}(k,:)];
        end
    end
end

std_all_points = std(all_points);
delete all_points

for i = 1:size(posteriors_grid_5_13_s2c2_r4_2000,1)
    for j = 1:size(posteriors_grid_5_13_s2c2_r4_2000,2)
        
        detection_grid_5_13_s2c2_r4_2000{i,j} = zeros(size(posteriors_grid_5_13_s2c2_r4_2000{i,j}));
        
        for k = 1:size(posteriors_grid_5_13_s2c2_r4_2000{i,j},1)
            [~, event_times] = findpeaks(posteriors_grid_5_13_s2c2_r4_2000{i,j}(k,:),'MinPeakHeight',1.0*std_all_points,'MinPeakDistance',60);
            detection_grid_5_13_s2c2_r4_2000{i,j}(k,event_times) = 50;
        end
    end
end

figure; compare_trace_stack_grid_overlap({detection_grid_5_13_s2c2_r4_2000,traces_5_13_s2c2_r4_tracegrid},3,1,[],0,{'L4','L5'},1)
%%

max_xcorr = cell(size(detection_grid_5_13_s2c1_r4_2000));
mad_xcorr_lag = cell(size(detection_grid_5_13_s2c1_r4_2000));

max_xcorr_img_shuff = zeros(size(detection_grid_5_13_s2c1_r4_2000));
mad_xcorr_lag_img = zeros(size(detection_grid_5_13_s2c1_r4_2000));

shuffle_rows = randperm(size(detection_grid_5_13_s2c1_r4_2000,1));
shuffle_cols = randperm(size(detection_grid_5_13_s2c1_r4_2000,2));
shuffle_trials = randperm(size(detection_grid_5_13_s2c1_r4_2000{1,1},1));
while isequal(1:size(detection_grid_5_13_s2c1_r4_2000{1,1},1),shuffle_trials)
    shuffle_trials = randperm(size(detection_grid_5_13_s2c1_r4_2000{1,1},1));
end

for i = 1:size(detection_grid_5_13_s2c1_r4_2000,1)
    for j = 1:size(detection_grid_5_13_s2c1_r4_2000,2)
        max_xcorr{i,j} = zeros(size(detection_grid_5_13_s2c1_r4_2000{i,j},1),1);
        mad_xcorr_lag{i,j} = zeros(size(detection_grid_5_13_s2c1_r4_2000{i,j},1),1);
        trace1_tmp = [];
        trace2_tmp = [];
        for k = 1:size(detection_grid_5_13_s2c1_r4_2000{i,j},1)
            trace1_tmp = [trace1_tmp smoothts(detection_grid_5_13_s2c1_r4_2000{i,j}(k,:),'g',100,10)];
            trace2_tmp = [trace2_tmp smoothts(detection_grid_5_13_s2c2_r4_2000{i,j}(shuffle_trials(k),:),'g',100,10)];
        end
        full_xcorr = xcorr(trace1_tmp,trace2_tmp,'coeff');
        [max_xcorr_img_shuff(i,j), mad_xcorr_lag_img(i,j)] = max(full_xcorr(9900:10100));
    end
end

%%
max_xcorr = cell(size(detection_grid_5_13_s2c1_r4_2000));
mad_xcorr_lag = cell(size(detection_grid_5_13_s2c1_r4_2000));

max_xcorr_img = zeros(size(detection_grid_5_13_s2c1_r4_2000));
mad_xcorr_lag_img = zeros(size(detection_grid_5_13_s2c1_r4_2000));

for i = 1:size(detection_grid_5_13_s2c1_r4_2000,1)
    for j = 1:size(detection_grid_5_13_s2c1_r4_2000,2)
        max_xcorr{i,j} = zeros(size(detection_grid_5_13_s2c1_r4_2000{i,j},1),1);
        mad_xcorr_lag{i,j} = zeros(size(detection_grid_5_13_s2c1_r4_2000{i,j},1),1);
        trace1_tmp = [];
        trace2_tmp = [];
        for k = 1:size(detection_grid_5_13_s2c1_r4_2000{i,j},1)
            trace1_tmp = [trace1_tmp smoothts(detection_grid_5_13_s2c1_r4_2000{i,j}(k,:),'g',100,10)];
            trace2_tmp = [trace2_tmp smoothts(detection_grid_5_13_s2c2_r4_2000{i,j}(k,:),'g',100,10)];
        end
        full_xcorr = xcorr(trace1_tmp,trace2_tmp,'coeff');
        [max_xcorr_img(i,j), mad_xcorr_lag_img(i,j)] = max(full_xcorr(9900:10100));
    end
end

%%

max_xcorr_img_sig = max_xcorr_img;
mad_xcorr_lag_img_sig = mad_xcorr_lag_img;
max_xcorr_img_sig(max_xcorr_img_sig < quantile(max_xcorr_img_shuff(:),.975)) = 0;
mad_xcorr_lag_img_sig(max_xcorr_img_sig < quantile(max_xcorr_img_shuff(:),.975)) = 0;

figure;
subplot(121)
imagesc(max_xcorr_img_sig)
colorbar
subplot(122)
% mad_xcorr_lag_img(mad_xcorr_lag_img >= 5060) = 0;
imagesc(mad_xcorr_lag_img_sig - 100)
% caxis([4940 5060])
colorbar

colormap hot