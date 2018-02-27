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

filenames = {'2_13_slice3_cell1.mat', '2_13_15_26_data.mat'
             '2_13_slice3_cell2.mat', '2_13_15_37_data.mat'
             '2_13_slice4_cell3.mat', '2_13_16_29_data.mat'
             '2_14_slice2_cell1.mat', '2_14_16_42_data.mat'
             '2_14_slice2_cell3.mat', '2_14_16_58_data.mat'
             '2_14_slice2_cell4.mat', '2_14_17_13_data.mat'
             '2_14_slice3_cell1.mat', '2_14_17_38_data.mat'};     
         
spike_trials = {1,1,1,1,1,1,1};
current_trials = {4,4,4,4,4,4,4};

colors = jet(4);
% colors(:) = 0;
% colors(size(filenames,1)+1:size(filenames,1)*2) = 1; %green
% z_pos = [-8 0 8];


%%


clear result_xy

pos_vs_cur_and_spike_time = figure;
pow_vs_cur_and_spike_time = figure;
new_fig = 0;
if ~exist('curr_vs_time','var')
    new_fig = 1;
    curr_vs_time = figure;
end
for j = 1:size(filenames,1)
    
    load(filenames{j,2});
    load(filenames{j,1}); 
    
    result_xy(j).quadrant = experiment_setup.quadrant;
    
    [result_xy(j).spike_traces, ~, this_seq, result_xy(j).spike_stim_traces] = get_traces(data,spike_trials{j});
    result_xy(j).spike_targ_pos = bsxfun(@minus,data.trial_metadata(spike_trials{j}).stim_key([this_seq.precomputed_target_index],:),experiment_setup.center_pos_um);
    
    [result_xy(j).current_traces, ~,this_seq, result_xy(j).current_stim_traces] = get_traces(data,current_trials{j});
    result_xy(j).current_targ_pos = bsxfun(@minus,data.trial_metadata(current_trials{j}).stim_key([this_seq.precomputed_target_index],:),experiment_setup.center_pos_um);
    result_xy(j).max_curr = result_xy(j).current_traces(:,1) - min(result_xy(j).current_traces,[],2);
    result_xy(j).spike_times = detect_peaks(-bsxfun(@minus,result_xy(j).spike_traces,result_xy(j).spike_traces(:,1)),20,30,0,Inf,-Inf,0,0,1);
    cell_spike_times = result_xy(j).spike_times;
    result_xy(j).spike_times = zeros(size(result_xy(j).spike_times));
    
    for i = 1:length(cell_spike_times)
        if ~isempty(cell_spike_times{i})
            result_xy(j).spike_times(i) = cell_spike_times{i};
        else
            result_xy(j).spike_times(i) = NaN;
        end
    end

    
    
    tested_pos = unique(result_xy(j).spike_targ_pos(:,1));
    if ~isempty(setdiff(tested_pos,unique(result_xy(j).current_targ_pos(:,1))))
        disp([filenames{j,1} ' spike and current pos sets differ'])
        continue
    elseif ~isempty(setdiff(tested_pos,unique(result_xy(j).spike_targ_pos(:,2))))
        disp([filenames{j,1} ' y and x pos sets differ'])
        continue
    elseif ~isempty(setdiff(tested_pos,unique(result_xy(j).current_targ_pos(:,2))))
        disp([filenames{j,1} ' y and x pos sets differ across spike curr'])
        continue
    end
    
    these_x_power = x_measurements(:,result_xy(j).quadrant);
    these_y_power = y_measurements(:,result_xy(j).quadrant);
    result_xy(j).spike_targ_power = zeros(size(result_xy(j).spike_times));
    result_xy(j).curr_targ_power = zeros(size(result_xy(j).max_curr));
    for i = 1:length(tested_pos)
        these_trials = result_xy(j).spike_targ_pos(:,1) == 0 & result_xy(j).spike_targ_pos(:,2) == tested_pos(i);
        result_xy(j).x_spike_time_means(i) = nanmean(result_xy(j).spike_times(these_trials));
        result_xy(j).x_spike_time_jitter(i) = nanstd(result_xy(j).spike_times(these_trials));
        result_xy(j).spike_targ_power(these_trials) = these_x_power(i);
        these_trials = result_xy(j).current_targ_pos(:,1) == 0 & result_xy(j).current_targ_pos(:,2) == tested_pos(i);
        result_xy(j).x_max_curr_means(i) = nanmean(result_xy(j).max_curr(these_trials));
        result_xy(j).curr_targ_power(these_trials) = these_x_power(i);
        
        these_trials = result_xy(j).spike_targ_pos(:,2) == 0 & result_xy(j).spike_targ_pos(:,1) == tested_pos(i);
        result_xy(j).y_spike_time_means(i) = nanmean(result_xy(j).spike_times(these_trials));
        result_xy(j).y_spike_time_jitter(i) = nanstd(result_xy(j).spike_times(these_trials));
        result_xy(j).spike_targ_power(these_trials) = these_y_power(i);
        these_trials = result_xy(j).current_targ_pos(:,2) == 0 & result_xy(j).current_targ_pos(:,1) == tested_pos(i);
        result_xy(j).y_max_curr_means(i) = nanmean(result_xy(j).max_curr(these_trials));
        result_xy(j).curr_targ_power(these_trials) = these_y_power(i);
    end
    
    figure(pow_vs_cur_and_spike_time)
    subplot(221)
%     these_trials = result_xy(j).spike_targ_pos(:,1) == 0;
    scatter(result_xy(j).spike_targ_power,result_xy(j).spike_times/20,[],colors(experiment_setup.quadrant,:));
    hold on
%     plot(tested_pos,result_xy(j).x_spike_time_means/20,'color',colors(experiment_setup.quadrant,:))
%     xlim([-20 20])
%     xlabel('Horizontal Distance (um)')
%     ylabel('Spike Time (msec)')
%     title('Horizontal Distance vs. Spike Time')
    subplot(122)
%     these_trials = result_xy(j).current_targ_pos(:,1) == 0;
    scatter(result_xy(j).curr_targ_power,result_xy(j).max_curr,[],colors(experiment_setup.quadrant,:));
%     plot(tested_pos,result_xy(j).x_max_curr_means,'color',colors(experiment_setup.quadrant,:))
    hold on
%     xlim([-20 20])
%     xlabel('Horizontal Distance (um)')
%     ylabel('Peak Current (pA)')
%     title('Horizontal Distance vs. Peak Current')
    
    figure(pos_vs_cur_and_spike_time)
    subplot(221)
    these_trials = result_xy(j).spike_targ_pos(:,1) == 0;
    scatter(result_xy(j).spike_targ_pos(these_trials,2),result_xy(j).spike_times(these_trials)/20,[],colors(experiment_setup.quadrant,:));
    hold on
    plot(tested_pos,result_xy(j).x_spike_time_means/20,'color',colors(experiment_setup.quadrant,:))
    xlim([-20 20])
    xlabel('Horizontal Distance (um)')
    ylabel('Spike Time (msec)')
    title('Horizontal Distance vs. Spike Time')
    subplot(222)
    these_trials = result_xy(j).current_targ_pos(:,1) == 0;
    scatter(result_xy(j).current_targ_pos(these_trials,2),result_xy(j).max_curr(these_trials),[],colors(experiment_setup.quadrant,:));
    plot(tested_pos,result_xy(j).x_max_curr_means,'color',colors(experiment_setup.quadrant,:))
    hold on
    xlim([-20 20])
    xlabel('Horizontal Distance (um)')
    ylabel('Peak Current (pA)')
    title('Horizontal Distance vs. Peak Current')
    subplot(223)
    these_trials = result_xy(j).spike_targ_pos(:,2) == 0;
    scatter(result_xy(j).spike_targ_pos(these_trials,1),result_xy(j).spike_times(these_trials)/20,[],colors(experiment_setup.quadrant,:));
    hold on
    plot(tested_pos,result_xy(j).y_spike_time_means/20,'color',colors(experiment_setup.quadrant,:))
    xlim([-20 20])
    xlabel('Vertical Distance (um)')
    ylabel('Spike Time (msec)')
    title('Vertical Distance vs. Spike Time')
    subplot(224)
    these_trials = result_xy(j).current_targ_pos(:,2) == 0;
    scatter(result_xy(j).current_targ_pos(these_trials,1),result_xy(j).max_curr(these_trials),[],colors(experiment_setup.quadrant,:));
    hold on
    plot(tested_pos,result_xy(j).y_max_curr_means,'color',colors(experiment_setup.quadrant,:))
    xlim([-20 20])
    xlabel('Vertical Distance (um)')
    ylabel('Peak Current (pA)')
    title('Vertical Distance vs. Peak Current')
    
    figure(curr_vs_time)
    subplot(121)
    [ordered_curr, curr_order_x] = sort(result_xy(j).x_max_curr_means);
    plot(result_xy(j).x_max_curr_means(curr_order_x),result_xy(j).x_spike_time_means(curr_order_x)/20,'-r','Linewidth',2);
    hold on
    [ordered_curr, curr_order_y] = sort(result_xy(j).y_max_curr_means);
    plot(result_xy(j).y_max_curr_means(curr_order_y),result_xy(j).y_spike_time_means(curr_order_y)/20,'-r','Linewidth',2);
    ylim([0 15])
    if new_fig
        xlabel('Peak Current Mean (pA)')
        ylabel('Spike Time Mean (msec)')
        title('Peak Current vs. Spike Time')
    end
    subplot(122)
    semilogy(result_xy(j).x_max_curr_means(curr_order_x),result_xy(j).x_spike_time_jitter(curr_order_x)/20,'-r','Linewidth',2);
    hold on
    semilogy(result_xy(j).y_max_curr_means(curr_order_y),result_xy(j).y_spike_time_jitter(curr_order_y)/20,'-r','Linewidth',2);
    if new_fig
        xlabel('Peak Current Mean (pA)')
        ylabel('Spike Time Std. Dev. (msec)')
        title('Peak Current vs. Spike Time Jitter')
    end

    
    
    
    
end

