%%
    
         
stacknames = {'/media/shababo/data/01232018images/s2c1-pre - 1_C0'
              '/media/shababo/data/01232018images/s2c1-pre - 2_C0'
              '/media/shababo/data/01232018images/s2c1-pre - 5_C0'
              '/media/shababo/data/01252018images/s2c1-pre_C0'
              '/media/shababo/data/01252018images/s2c1-pre - 1_C0'
              '/media/shababo/data/01252018images/s2c1-pre - 3_C0'
              '/media/shababo/data/01252018images/s2c1-pre - 4_C0'
              '/media/shababo/data/0126272018images/s2c1-pre_C0'
              '/media/shababo/data/0126272018images/s2c1-pre - 1_C0'
              '/media/shababo/data/0126272018images/s2c1-pre - 2_C0'
              '/media/shababo/data/0126272018images/s2c1-pre - 10_C0'
              '/media/shababo/data/0126272018images/s2c1-pre - 11_C0'
              '/media/shababo/data/0126272018images/s2c1-pre - 12_C0'};

filenames = {'2_24_slice1_cell1.mat', '2_24_13_57_data.mat'
             '2_24_slice1_cell2.mat', '2_24_14_20_data.mat'
             '2_24_slice2_cell1.mat', '2_24_14_47_data.mat'
             '2_26_slice1_cell1.mat', '2_26_13_20_data.mat'
             '2_26_slice1_cell2.mat', '2_26_13_48_data.mat'
             '2_26_slice2_cell1.mat', '2_26_14_11_data.mat'
             '2_26_slice2_cell2.mat', '2_26_14_27_data.mat'
             '2_26_slice3_cell1.mat', '2_26_15_10_data.mat'
             '2_26_slice3_cell2.mat', '2_26_15_32_data.mat'
             '2_26_slice3_cell3.mat', '2_26_15_53_data.mat'
             '2_26_slice3_cell4.mat', '2_26_16_15_data.mat'};     
         
spike_trials = {[],[],[],[],1,1,1,1,1,1,1};
current_trials = {3:8,3:8,3,3:8,4:9,4:6,4:9,4:9,4:9,4:9,4:9};
intrinsic_trials = {2,2,2,2,3,3,3,3,3,3,3};
z_center = z_pos(i)*ones(size(filenames,1),1);

colors = jet(4);
% colors(:) = 0;
% colors(size(filenames,1)+1:size(filenames,1)*2) = 1; %green
% z_pos = [-8 0 8];


%%


clear result_shape

for j = 1:size(filenames,1)
    
    load(filenames{j,2});
    load(filenames{j,1}); 
    
    result_shape(j).quadrant = experiment_setup.quadrant;
    
    if ~isempty(spike_trials{j})
        [result_shape(j).spike_traces, ~, this_seq, result_shape(j).spike_stim_traces] = get_traces(data,spike_trials{j});
        result_shape(j).spike_targ_pos = bsxfun(@minus,data.trial_metadata(spike_trials{j}).stim_key([this_seq.precomputed_target_index],:),experiment_setup.center_pos_um);
    
        result_shape(j).spike_times = detect_peaks(-bsxfun(@minus,result_shape(j).spike_traces,result_shape(j).spike_traces(:,1)),20,30,0,Inf,-Inf,0,0,1);
        cell_spike_times = result_shape(j).spike_times;
        result_shape(j).spike_times = zeros(size(result_shape(j).spike_times));

        for i = 1:length(cell_spike_times)
            if ~isempty(cell_spike_times{i})
                result_shape(j).spike_times(i) = cell_spike_times{i};
            else
                result_shape(j).spike_times(i) = NaN;
            end
        end
    end
    
    [result_shape(j).current_traces, ~,this_seq, result_shape(j).current_stim_traces] = get_traces(data,current_trials{j});
    result_shape(j).current_targ_pos = ...
        [bsxfun(@minus,data.trial_metadata(current_trials{j}(1)).stim_key([this_seq.precomputed_target_index],:),experiment_setup.center_pos_um) ...
            [this_seq.piezo_z]'];
    result_shape(j).max_curr = result_shape(j).current_traces(:,1) - min(result_shape(j).current_traces,[],2);

    
%     
%     tested_pos = unique(result_shape(j).spike_targ_pos(:,1));
%     if ~isempty(setdiff(tested_pos,unique(result_shape(j).current_targ_pos(:,1))))
%         disp([filenames{j,1} ' spike and current pos sets differ'])
%         continue
%     elseif ~isempty(setdiff(tested_pos,unique(result_shape(j).spike_targ_pos(:,2))))
%         disp([filenames{j,1} ' y and x pos sets differ'])
%         continue
%     elseif ~isempty(setdiff(tested_pos,unique(result_shape(j).current_targ_pos(:,2))))
%         disp([filenames{j,1} ' y and x pos sets differ across spike curr'])
%         continue
%     end
%     
%     these_x_power = x_measurements(:,result_shape(j).quadrant);
%     these_y_power = y_measurements(:,result_shape(j).quadrant);
%     result_shape(j).spike_targ_power = zeros(size(result_shape(j).spike_times));
%     result_shape(j).curr_targ_power = zeros(size(result_shape(j).max_curr));
%     for i = 1:length(tested_pos)
%         these_trials = result_shape(j).spike_targ_pos(:,1) == 0 & result_shape(j).spike_targ_pos(:,2) == tested_pos(i);
%         result_shape(j).x_spike_time_means(i) = nanmean(result_shape(j).spike_times(these_trials));
%         result_shape(j).x_spike_time_jitter(i) = nanstd(result_shape(j).spike_times(these_trials));
%         result_shape(j).spike_targ_power(these_trials) = these_x_power(i);
%         these_trials = result_shape(j).current_targ_pos(:,1) == 0 & result_shape(j).current_targ_pos(:,2) == tested_pos(i);
%         result_shape(j).x_max_curr_means(i) = nanmean(result_shape(j).max_curr(these_trials));
%         result_shape(j).curr_targ_power(these_trials) = these_x_power(i);
%         
%         these_trials = result_shape(j).spike_targ_pos(:,2) == 0 & result_shape(j).spike_targ_pos(:,1) == tested_pos(i);
%         result_shape(j).y_spike_time_means(i) = nanmean(result_shape(j).spike_times(these_trials));
%         result_shape(j).y_spike_time_jitter(i) = nanstd(result_shape(j).spike_times(these_trials));
%         result_shape(j).spike_targ_power(these_trials) = these_y_power(i);
%         these_trials = result_shape(j).current_targ_pos(:,2) == 0 & result_shape(j).current_targ_pos(:,1) == tested_pos(i);
%         result_shape(j).y_max_curr_means(i) = nanmean(result_shape(j).max_curr(these_trials));
%         result_shape(j).curr_targ_power(these_trials) = these_y_power(i);
%     end
%     
%     figure(pow_vs_cur_and_spike_time)
%     subplot(221)
% %     these_trials = result_shape(j).spike_targ_pos(:,1) == 0;
%     scatter(result_shape(j).spike_targ_power,result_shape(j).spike_times/20,[],colors(experiment_setup.quadrant,:));
%     hold on
% %     plot(tested_pos,result_shape(j).x_spike_time_means/20,'color',colors(experiment_setup.quadrant,:))
% %     xlim([-20 20])
% %     xlabel('Horizontal Distance (um)')
% %     ylabel('Spike Time (msec)')
% %     title('Horizontal Distance vs. Spike Time')
%     subplot(122)
% %     these_trials = result_shape(j).current_targ_pos(:,1) == 0;
%     scatter(result_shape(j).curr_targ_power,result_shape(j).max_curr,[],colors(experiment_setup.quadrant,:));
% %     plot(tested_pos,result_shape(j).x_max_curr_means,'color',colors(experiment_setup.quadrant,:))
%     hold on
% %     xlim([-20 20])
% %     xlabel('Horizontal Distance (um)')
% %     ylabel('Peak Current (pA)')
% %     title('Horizontal Distance vs. Peak Current')
%     
%     figure(pos_vs_cur_and_spike_time)
%     subplot(221)
%     these_trials = result_shape(j).spike_targ_pos(:,1) == 0;
%     scatter(result_shape(j).spike_targ_pos(these_trials,2),result_shape(j).spike_times(these_trials)/20,[],colors(experiment_setup.quadrant,:));
%     hold on
%     plot(tested_pos,result_shape(j).x_spike_time_means/20,'color',colors(experiment_setup.quadrant,:))
%     xlim([-20 20])
%     xlabel('Horizontal Distance (um)')
%     ylabel('Spike Time (msec)')
%     title('Horizontal Distance vs. Spike Time')
%     subplot(222)
%     these_trials = result_shape(j).current_targ_pos(:,1) == 0;
%     scatter(result_shape(j).current_targ_pos(these_trials,2),result_shape(j).max_curr(these_trials),[],colors(experiment_setup.quadrant,:));
%     plot(tested_pos,result_shape(j).x_max_curr_means,'color',colors(experiment_setup.quadrant,:))
%     hold on
%     xlim([-20 20])
%     xlabel('Horizontal Distance (um)')
%     ylabel('Peak Current (pA)')
%     title('Horizontal Distance vs. Peak Current')
%     subplot(223)
%     these_trials = result_shape(j).spike_targ_pos(:,2) == 0;
%     scatter(result_shape(j).spike_targ_pos(these_trials,1),result_shape(j).spike_times(these_trials)/20,[],colors(experiment_setup.quadrant,:));
%     hold on
%     plot(tested_pos,result_shape(j).y_spike_time_means/20,'color',colors(experiment_setup.quadrant,:))
%     xlim([-20 20])
%     xlabel('Vertical Distance (um)')
%     ylabel('Spike Time (msec)')
%     title('Vertical Distance vs. Spike Time')
%     subplot(224)
%     these_trials = result_shape(j).current_targ_pos(:,2) == 0;
%     scatter(result_shape(j).current_targ_pos(these_trials,1),result_shape(j).max_curr(these_trials),[],colors(experiment_setup.quadrant,:));
%     hold on
%     plot(tested_pos,result_shape(j).y_max_curr_means,'color',colors(experiment_setup.quadrant,:))
%     xlim([-20 20])
%     xlabel('Vertical Distance (um)')
%     ylabel('Peak Current (pA)')
%     title('Vertical Distance vs. Peak Current')
%     
%     figure(curr_vs_time)
%     subplot(121)
%     [ordered_curr, curr_order_x] = sort(result_shape(j).x_max_curr_means);
%     plot(result_shape(j).x_max_curr_means(curr_order_x),result_shape(j).x_spike_time_means(curr_order_x)/20,'-r','Linewidth',2);
%     hold on
%     [ordered_curr, curr_order_y] = sort(result_shape(j).y_max_curr_means);
%     plot(result_shape(j).y_max_curr_means(curr_order_y),result_shape(j).y_spike_time_means(curr_order_y)/20,'-r','Linewidth',2);
%     ylim([0 15])
%     if new_fig
%         xlabel('Peak Current Mean (pA)')
%         ylabel('Spike Time Mean (msec)')
%         title('Peak Current vs. Spike Time')
%     end
%     subplot(122)
%     semilogy(result_shape(j).x_max_curr_means(curr_order_x),result_shape(j).x_spike_time_jitter(curr_order_x)/20,'-r','Linewidth',2);
%     hold on
%     semilogy(result_shape(j).y_max_curr_means(curr_order_y),result_shape(j).y_spike_time_jitter(curr_order_y)/20,'-r','Linewidth',2);
%     if new_fig
%         xlabel('Peak Current Mean (pA)')
%         ylabel('Spike Time Std. Dev. (msec)')
%         title('Peak Current vs. Spike Time Jitter')
%     end

    
    
    
    
end

%%

z_pos = [160 200 240];
figure
count = 1;
for i = 1:length(z_pos)
    

    for j = 1:size(filenames,1)
        count
        subplot(length(z_pos),size(filenames,1),count)
        scatter3(result_shape(j).current_targ_pos(result_shape(j).current_targ_pos(:,3) == z_pos(i) & result_shape(j).max_curr < 2000,1), ...
                 result_shape(j).current_targ_pos(result_shape(j).current_targ_pos(:,3) == z_pos(i) & result_shape(j).max_curr < 2000,2), ...
                 result_shape(j).max_curr(result_shape(j).current_targ_pos(:,3) == z_pos(i) & result_shape(j).max_curr < 2000),[], ...
                 result_shape(j).max_curr(result_shape(j).current_targ_pos(:,3) == z_pos(i) & result_shape(j).max_curr < 2000)/...
                 max(result_shape(j).max_curr(result_shape(j).max_curr < 2000)))
             zlim([0 2000])
         count = count + 1;

    end
end
