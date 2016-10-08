function see_grid(

trial_ind = 7; traces = [data.sweeps{trial_ind}(:,1)'; data.sweeps{trial_ind}(:,2)'];
stim = data.sweeps{trial_ind}(:,4)' > .35; sum(diff(stim) == 1)
figure; plot(stim)
stim(1:384000) = 0;
maps_1007_s5c1_trial7 = build_slm_maps(traces,stim,map_index,.1*20000);
trial_ind = 8; traces = [data.sweeps{trial_ind}(:,1)'; data.sweeps{trial_ind}(:,2)'];
stim = data.sweeps{trial_ind}(:,4)' > .35; sum(diff(stim) == 1)
maps_1007_s5c1_trial8 = build_slm_maps(traces,stim,map_index,.1*20000);
map_l5_s5c1_trials12 = stack_grids({maps_1007_s5c1_trial1{1},maps_1007_s5c1_trial2{1}});
map_l5_s5c1_trials78 = stack_grids({maps_1007_s5c1_trial7{1},maps_1007_s5c1_trial8{1}});
map_l4_s5c1_trials78 = stack_grids({maps_1007_s5c1_trial7{2},maps_1007_s5c1_trial8{2}});
figure;compare_trace_stack_grid({map_l4_s5c1_trials78,map_l5_s5c1_trials78},6,1,[],0,{'raw','detected events'})
figure; corr_image_l5_x = get_corr_image(map_l5_s5c1_trials78,1);
figure; corr_image_l4_x = get_corr_image(map_l4_s5c1_trials78,1);
figure; corr_image_l5_x = get_corr_image(map_l5_s5c1_trials12,1);
figure; corr_image_l5_x = get_corr_image(map_l4_s5c1_trials12,1);
trial_ind = 4; traces = [data.sweeps{trial_ind}(:,1)'; data.sweeps{trial_ind}(:,2)'];
stim = data.sweeps{trial_ind}(:,4)' > .35; sum(diff(stim) == 1)
figure; plot(stim)
stim(1:233600) = 0;
maps_1007_s5c1_trial4 = build_slm_maps(traces,stim,map_index,.1*20000);
trial_ind = 5; traces = [data.sweeps{trial_ind}(:,1)'; data.sweeps{trial_ind}(:,2)'];
stim = data.sweeps{trial_ind}(:,4)' > .35; sum(diff(stim) == 1)
maps_1007_s5c1_trial5 = build_slm_maps(traces,stim,map_index,.1*20000);
map_l4_s5c1_trials45 = stack_grids({maps_1007_s5c1_trial4{2},maps_1007_s5c1_trial5{2}});
map_l5_s5c1_trials45 = stack_grids({maps_1007_s5c1_trial4{1},maps_1007_s5c1_trial5{1}});
figure; corr_image_l4_x = get_corr_image(map_l4_s5c1_trials45,1);
figure; corr_image_l4_x = get_corr_image(map_l5_s5c1_trials45,1);