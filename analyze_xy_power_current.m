%%
datapath = '/media/shababo/data/'

filenames = {'2_13_slice3_cell1.mat', '2_13_15_26_data.mat', '/media/shababo/data/02132018images/s2c1-pre - 3_C0' % good %P19
             '2_13_slice3_cell2.mat', '2_13_15_37_data.mat', '/media/shababo/data/02132018images/s2c1-pre - 4_C0' % good
             '2_13_slice4_cell3.mat', '2_13_16_29_data.mat', '/media/shababo/data/02132018images/s2c1-pre - 8_C0' % good
             '2_14_slice2_cell1.mat', '2_14_16_42_data.mat', '/media/shababo/data/02142018images/s2c1-pre - 3_C0' % good %P20
             '2_14_slice2_cell3.mat', '2_14_16_58_data.mat', '/media/shababo/data/02142018images/s2c1-pre - 5_C0' % good
             '2_14_slice2_cell4.mat', '2_14_17_13_data.mat', '/media/shababo/data/02142018images/s2c1-pre - 6_C0' % good
             '2_14_slice3_cell1.mat', '2_14_17_38_data.mat', '/media/shababo/data/02142018images/s2c1-pre - 7_C0'};  % meh   
         
spike_trials = {1,1,1,1,1,1,1};
current_trials = {4,4,4,4,4,4,4};

colors = jet(4);
% colors(:) = 0;
% colors(size(filenames,1)+1:size(filenames,1)*2) = 1; %green
% z_pos = [-8 0 8];
colors = jet(7);
figure_for_cosyne = figure;
%%


clear result_xy

pos_vs_cur_and_spike_time = figure;
new_fig = 0;
if ~exist('curr_vs_time','var')
    new_fig = 1;
    curr_vs_time = figure;
end
do_detect = 0;

for j = 1:size(filenames,1)%find([result_xy_bu.quadrant] == 1)%

    
    load([datapath filenames{j,2})];
    load([datapath filenames{j,1})]; 
    
    result_xy(j).quadrant = experiment_setup.quadrant;
    
    if ~isempty(filenames{j,3})
        
        [nuclear_locs,fluor_vals,nuclear_locs_image_coord] = detect_nuclei(filenames{j,3},[],[],[],do_detect,[],0);
        offsets = nuclear_locs - [experiment_setup.center_pos_um(1:2) 20];

        [targ_error, index] = min(sqrt(sum(offsets.^2,2)));
        result_xy(j).fluor_val = fluor_vals(index);
        result_xy(j).cell_pos = nuclear_locs(index,:);
    else
        result_xy(j).fluor_val = NaN;
        result_xy(j).cell_pos = [experiment_setup.center_pos_um(1:2) experiment_setup.piezo_center];

    end
    
    [result_xy(j).spike_traces, ~, this_seq, result_xy(j).spike_stim_traces] = get_traces(data,spike_trials{j});
    result_xy(j).spike_targ_pos = bsxfun(@minus,data.trial_metadata(spike_trials{j}).stim_key([this_seq.precomputed_target_index],:),result_xy(j).cell_pos);%experiment_setup.center_pos_um);
    
    [result_xy(j).current_traces, ~,this_seq, result_xy(j).current_stim_traces] = get_traces(data,current_trials{j});
    result_xy(j).current_targ_pos = bsxfun(@minus,data.trial_metadata(current_trials{j}).stim_key([this_seq.precomputed_target_index],:),result_xy(j).cell_pos);%experiment_setup.center_pos_um);
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

    
    
    tested_pos_y = unique(result_xy(j).spike_targ_pos(:,1));
    tested_pos_x = unique(result_xy(j).spike_targ_pos(:,2));
%     if ~isempty(setdiff(tested_pos,unique(result_xy(j).current_targ_pos(:,1))))
%         disp([filenames{j,1} ' spike and current pos sets differ'])
%         continue
%     elseif ~isempty(setdiff(tested_pos,unique(result_xy(j).spike_targ_pos(:,2))))
%         disp([filenames{j,1} ' y and x pos sets differ'])
%         continue
%     elseif ~isempty(setdiff(tested_pos,unique(result_xy(j).current_targ_pos(:,2))))
%         disp([filenames{j,1} ' y and x pos sets differ across spike curr'])
%         continue
%     end
    
    
    
%     tested_pos = [-15 -10 -7 -4 -2 0 2 4 7 10 15];
    result_xy(j).these_x_power = x_measurements(:,result_xy(j).quadrant);%.*shape_template(sub2ind(size(shape_template),round(36 + tested_pos_y(6))*ones(size(tested_pos_x)),round(tested_pos_x)+36));
    result_xy(j).these_y_power = y_measurements(:,result_xy(j).quadrant);%.*shape_template(sub2ind(size(shape_template),round(tested_pos_y)+36,round(36 + tested_pos_x(6))*ones(size(tested_pos_y))));
    result_xy(j).spike_targ_power = zeros(size(result_xy(j).spike_times));
    result_xy(j).curr_targ_power = zeros(size(result_xy(j).max_curr));
    for i = 1:length(tested_pos_x)
        these_trials = result_xy(j).spike_targ_pos(:,1) == tested_pos_y(6) & result_xy(j).spike_targ_pos(:,2) == tested_pos_x(i);
        result_xy(j).x_spike_time_means(i) = nanmean(result_xy(j).spike_times(these_trials));
        result_xy(j).x_spike_time_jitter(i) = nanstd(result_xy(j).spike_times(these_trials));
        result_xy(j).spike_targ_power(these_trials) = result_xy(j).these_x_power(i);
        these_trials = result_xy(j).current_targ_pos(:,1) == tested_pos_y(6) & result_xy(j).current_targ_pos(:,2) == tested_pos_x(i);
        result_xy(j).x_max_curr_means(i) = nanmean(result_xy(j).max_curr(these_trials));
%         result_xy(j).x_max_curr_means_pownorm(i) = nanmean(result_xy(j).max_curr(these_trials))./result_xy(j).curr_targ_power(these_trials);
        result_xy(j).curr_targ_power(these_trials) = result_xy(j).these_x_power(i);
    end
    for i = 1:length(tested_pos_y)
        these_trials = result_xy(j).spike_targ_pos(:,2) == tested_pos_x(6) & result_xy(j).spike_targ_pos(:,1) == tested_pos_y(i);
        result_xy(j).y_spike_time_means(i) = nanmean(result_xy(j).spike_times(these_trials));
        result_xy(j).y_spike_time_jitter(i) = nanstd(result_xy(j).spike_times(these_trials));
        result_xy(j).spike_targ_power(these_trials) = result_xy(j).these_y_power(i);
        these_trials = result_xy(j).current_targ_pos(:,2) == tested_pos_x(6) & result_xy(j).current_targ_pos(:,1) == tested_pos_y(i);
        result_xy(j).y_max_curr_means(i) = nanmean(result_xy(j).max_curr(these_trials));
%         result_xy(j).y_max_curr_means_pownorm(i) = nanmean(result_xy(j).max_curr(these_trials))./result_xy(j).curr_targ_power(these_trials);
        result_xy(j).curr_targ_power(these_trials) = result_xy(j).these_y_power(i);
    end
    
    trial_hologram_power = result_xy(j).spike_targ_power;
    trial_hologram_position = result_xy(j).spike_targ_pos;
    trial_shape_gain = shape_template(sub2ind(size(shape_template),round(trial_hologram_position(:,1)+36),round(trial_hologram_position(:,1)+36)));
    result_xy(j).spatial_adj_power = trial_hologram_power .* trial_shape_gain;
    
%     figure(pow_vs_cur_and_spike_time)
%     subplot(221)
% %     these_trials = result_xy(j).spike_targ_pos(:,1) == 0;
%     scatter(result_xy(j).spike_targ_power,result_xy(j).spike_times/20,[],colors(j,:));
%     hold on
% %     plot(tested_pos,result_xy(j).x_spike_time_means/20,'color',colors(j,:))
% %     xlim([-20 20])
% %     xlabel('Horizontal Distance (um)')
% %     ylabel('Spike Time (msec)')
% %     title('Horizontal Distance vs. Spike Time')
%     subplot(122)
% %     these_trials = result_xy(j).current_targ_pos(:,1) == 0;
%     scatter(result_xy(j).curr_targ_power,result_xy(j).max_curr,[],colors(j,:));
% %     plot(tested_pos,result_xy(j).x_max_curr_means,'color',colors(j,:))
%     hold on
% %     xlim([-20 20])
% %     xlabel('Horizontal Distance (um)')
% %     ylabel('Peak Current (pA)')
% %     title('Horizontal Distance vs. Peak Current')
    
    figure(figure_for_cosyne)
%     subplot(221)
    subplot(2,size(filenames,1),j)
    yyaxis left
    these_trials = result_xy(j).spike_targ_pos(:,1) == tested_pos_y(6);
%     scatter(result_xy(j).spike_targ_pos(these_trials,2),result_xy(j).spike_times(these_trials)/20,[],colors(j,:));
    hold on
    plot(tested_pos_x,result_xy(j).x_spike_time_means/20,'color',colors(2,:),'linewidth',2)
    ylim([0 14])
    xlim([-20 20])
    xlabel('Horizontal Distance (um)')
    ylabel('Spike Time (msec)')
    title('Horizontal Distance vs. Spike Time')
%     subplot(222)
    yyaxis right
    these_trials = result_xy(j).current_targ_pos(:,1) == tested_pos_y(6);
%     scatter(result_xy(j).current_targ_pos(these_trials,2),result_xy(j).max_curr(these_trials),[],colors(j,:));
    hold on
    [~,this_zero_pos] = min(abs(tested_pos_x));
    scaling = result_xy(j).x_max_curr_means(this_zero_pos);
    [~,this_zero_pos] = min(abs(result_xy(j).current_targ_pos(these_trials,2)));
    these_powers = result_xy(j).spatial_adj_power(these_trials);

%     scaling = these_powers(this_zero_pos)*gain_mle(28+j)*1000;
%     plot(-20:20,scaling*shape_template(sub2ind(size(shape_template),36*ones(size(-20:20)),(-20:20)+36)),'color',colors(j,:),'linewidth',1);
    plot(tested_pos_x,result_xy(j).x_max_curr_means,'color',colors(6,:),'linewidth',2)
    ylim([0 1600])
    xlim([-20 20])
    xlabel('Horizontal Distance (um)')
    ylabel('Peak Current (pA)')
    title('Horizontal Distance vs. Peak Current')
%     subplot(223)
    subplot(2,size(filenames,1),size(filenames,1) + j)
    yyaxis left
    these_trials = result_xy(j).spike_targ_pos(:,2) == tested_pos_x(6);
%     scatter(result_xy(j).spike_targ_pos(these_trials,1),result_xy(j).spike_times(these_trials)/20,[],colors(j,:));
    hold on
    plot(tested_pos_y,result_xy(j).y_spike_time_means/20,'color',colors(2,:),'linewidth',2)
    ylim([0 14])
    xlim([-20 20])
    xlabel('Vertical Distance (um)')
    ylabel('Spike Time (msec)')
    title('Vertical Distance vs. Spike Time')
%     subplot(224)
    yyaxis right
    these_trials = result_xy(j).current_targ_pos(:,2) == tested_pos_x(6);
%     scatter(result_xy(j).current_targ_pos(these_trials,1),result_xy(j).max_curr(these_trials),[],colors(j,:));
    hold on
    [~,this_zero_pos] = min(abs(tested_pos_y))
    scaling = result_xy(j).y_max_curr_means(this_zero_pos);
    [~,this_zero_pos] = min(abs(result_xy(j).current_targ_pos(these_trials,1)));
    these_powers = result_xy(j).spatial_adj_power(these_trials);

%     scaling = these_powers(this_zero_pos)*gain_mle(28+j)*1000;
%     plot(-20:20,scaling*shape_template(sub2ind(size(shape_template),(-20:20)+36,36*ones(size(-20:20)))),'color',colors(j,:),'linewidth',1)
    plot(tested_pos_y,result_xy(j).y_max_curr_means,'color',colors(6,:),'linewidth',2)
    ylim([0 1600])
    xlim([-20 20])
    xlabel('Vertical Distance (um)')
    ylabel('Peak Current (pA)')
    title('Vertical Distance vs. Peak Current')
%     
%     figure(curr_vs_time)
%     subplot(121)
%     [ordered_curr, curr_order_x] = sort(result_xy(j).x_max_curr_means);
%     plot(result_xy(j).x_max_curr_means(curr_order_x),result_xy(j).x_spike_time_means(curr_order_x)/20,'-r','Linewidth',2);
%     hold on
%     [ordered_curr, curr_order_y] = sort(result_xy(j).y_max_curr_means);
%     plot(result_xy(j).y_max_curr_means(curr_order_y),result_xy(j).y_spike_time_means(curr_order_y)/20,'-r','Linewidth',2);
%     ylim([0 15])
%     if new_fig
%         xlabel('Peak Current Mean (pA)')
%         ylabel('Spike Time Mean (msec)')
%         title('Peak Current vs. Spike Time')
%     end
%     subplot(122)
%     semilogy(result_xy(j).x_max_curr_means(curr_order_x),result_xy(j).x_spike_time_jitter(curr_order_x)/20,'-r','Linewidth',2);
%     hold on
%     semilogy(result_xy(j).y_max_curr_means(curr_order_y),result_xy(j).y_spike_time_jitter(curr_order_y)/20,'-r','Linewidth',2);
%     if new_fig
%         xlabel('Peak Current Mean (pA)')
%         ylabel('Spike Time Std. Dev. (msec)')
%         title('Peak Current vs. Spike Time Jitter')
%     end

    
    
    
    
end


%% adjust powers


% shape_template = mean_cell(:,:,2)';
% 
% for i_cell = 1:length(result_xy)
%     
%     trial_hologram_power = result_xy(i_cell).spike_targ_power;
%     trial_hologram_position = result_xy(i_cell).spike_targ_pos;
%     trial_shape_gain = shape_template(sub2ind(size(shape_template),round(trial_hologram_position(:,1)+36),round(trial_hologram_position(:,1)+36)));
%     result_xy(i_cell).spatial_adj_power = trial_hologram_power .* trial_shape_gain;
%     
% end

%% nuc detect

% do_detect = 0;
% for j = 1:7
%     if ~isempty(filenames{j,3})
%         
%         load(filenames{j,2})
%         [nuclear_locs,fluor_vals,nuclear_locs_image_coord] = detect_nuclei(filenames{j,3},[],[],[],do_detect,[],0);
%         offsets = nuclear_locs - [experiment_setup.center_pos_um(1:2) 20];
% 
%         [targ_error, index] = min(sqrt(sum(offsets.^2,2)));
%         result_xy(j).fluor_val = fluor_vals(index);
%         result_xy(j).cell_pos = nuclear_locs(index,:);
%         result_xy(j).exp_cell_pos = [experiment_setup.center_pos_um(1:2) 20];
%         result_xy(j).err_cell_pos = result_xy(j).exp_cell_pos - result_xy(j).cell_pos;
%         result_xy(j).err_cell_pos_norm = norm(result_xy(j).exp_cell_pos - result_xy(j).cell_pos);
%         
%     else
%         
%         result_xy(j).fluor_val = NaN;
%         result_xy(j).cell_pos = NaN;
% 
%     end
% %     title(['Cell: ' num2str(j)])
% end

%%

figure

for j = 1:size(filenames,1)
    subplot(211)
    these_trials = result_xy(j).current_targ_pos(:,1) == tested_pos_y(6);
    %     scatter(result_xy(j).current_targ_pos(these_trials,2),result_xy(j).max_curr(these_trials),[],colors(j,:));
    hold on
    [~,this_zero_pos] = min(abs(tested_pos_x));
    scaling = result_xy(j).x_max_curr_means(this_zero_pos);
    [~,this_zero_pos] = min(abs(result_xy(j).current_targ_pos(these_trials,2)));
    these_powers = result_xy(j).spatial_adj_power(these_trials);
    %     scaling = these_powers(this_zero_pos)*gain_mle(28+j)*1000;
    %     plot(-20:20,scaling*shape_template(sub2ind(size(shape_template),36*ones(size(-20:20)),(-20:20)+36)),'color',colors(j,:),'linewidth',1);
    plot(tested_pos_x,result_xy(j).x_max_curr_means/max(result_xy(j).x_max_curr_means),'color',colors(result_xy(j).quadrant*2-1,:),'linewidth',2)
    ylim([0 1.2])
    xlim([-20 20])
    xlabel('Horizontal Distance (um)')
    ylabel('Peak Current (pA)')
    title('Horizontal Distance vs. Peak Current')
    subplot(212)
    these_trials = result_xy(j).current_targ_pos(:,2) == tested_pos_x(6);
    %     scatter(result_xy(j).current_targ_pos(these_trials,1),result_xy(j).max_curr(these_trials),[],colors(j,:));
    hold on
    [~,this_zero_pos] = min(abs(tested_pos_y))
    scaling = result_xy(j).y_max_curr_means(this_zero_pos);
    [~,this_zero_pos] = min(abs(result_xy(j).current_targ_pos(these_trials,1)));
    these_powers = result_xy(j).spatial_adj_power(these_trials);
    %     scaling = these_powers(this_zero_pos)*gain_mle(28+j)*1000;
    %     plot(-20:20,scaling*shape_template(sub2ind(size(shape_template),(-20:20)+36,36*ones(size(-20:20)))),'color',colors(j,:),'linewidth',1)
    plot(tested_pos_y,result_xy(j).y_max_curr_means/max(result_xy(j).y_max_curr_means),'color',colors(result_xy(j).quadrant*2-1,:),'linewidth',2)
    ylim([0 1.2])
    xlim([-20 20])
    xlabel('Vertical Distance (um)')
    ylabel('Peak Current (pA)')
    title('Vertical Distance vs. Peak Current')
end






















