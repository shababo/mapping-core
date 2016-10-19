function [map_ch1,map_ch2,corr_ch1,corr_ch2] = see_grid(data,trials,map_index)

traces_ch1 = cell(1,length(trials));
traces_ch2 = cell(1,length(trials));
for i = 1:length(trials)
    
    trial_ind = trials(i); 
    traces = [data.sweeps{trial_ind}(:,1)'; data.sweeps{trial_ind}(:,2)'];
    stim = data.sweeps{trial_ind}(:,4)' > .35; sum(diff(stim) == 1)
    stim_starts = find(diff(stim) == 1);
    if length(stim_starts) ~= 1323
        figure; plot(stim)
        stim_starts(1:length(stim_starts)-1323) = [];
    end
    
    maps = build_slm_maps(traces,stim_starts,map_index,.1*20000);
    traces_ch1{i} = maps{1};
    traces_ch2{i} = maps{2};
end
map_ch1 = stack_grids(traces_ch1);
map_ch2 = stack_grids(traces_ch2);

figure;compare_trace_stack_grid({map_ch1,map_ch2},Inf,1,[],0,{'raw','detected events'})
figure;compare_trace_stack_grid({map_ch1,map_ch2},Inf,1,[],1,{'raw','detected events'})
figure;compare_trace_stack_grid_overlap({map_ch1,map_ch2},Inf,1,[],0,{'L4','L5'},1)

corr_ch1 = get_corr_image(map_ch1,0,0);
corr_ch2 = get_corr_image(map_ch2,0,0);

figure; 
subplot(121); imagesc(corr_ch1); %caxis([0 1])
subplot(122); imagesc(corr_ch2); %caxis([0 1])

colormap hot

charge_map_ch1 = get_charge_map(map_ch1);
charge_map_ch2 = get_charge_map(map_ch2);

figure; 
subplot(121); imagesc(charge_map_ch1); %caxis([0 1])
subplot(122); imagesc(charge_map_ch2); %caxis([0 1])

colormap hot