%%
     

% filenames = {'1_31_slice1_cell2.mat', '1_31_15_3_data.mat', '/media/shababo/data/01312018images/s2c1-pre - 2_C0'
%              '1_31_slice1_cell3.mat', '1_31_15_13_data.mat', '/media/shababo/data/01312018images/s2c1-pre - 3_C0'
%              '2_21_slice1_cell1.mat', '2_21_16_57_data.mat', '/media/shababo/data/0221232018_zdata/s2c1-pre - 1_C0'
%              '2_22_slice1_cell1.mat', '2_22_15_17_data.mat', '/media/shababo/data/0221232018_zdata/s2c1-pre - 2_C0'
%              '2_23_slice1_cell1.mat', '2_23_12_15_data.mat', '/media/shababo/data/0221232018_zdata/s2c1-pre - 4_C0'
%              '2_23_slice1_cell2.mat', '2_23_12_27_data.mat', '/media/shababo/data/0221232018_zdata/s2c1-pre - 5_C0'
%              '2_23_slice2_cell1.mat', '2_23_12_59_data.mat', '/media/shababo/data/0221232018_zdata/s2c1-pre - 6_C0'
%              '2_23_slice2_cell2.mat', '2_23_13_24_data.mat', '/media/shababo/data/0221232018_zdata/s2c1-pre - 7_C0'
%              '2_23_slice3_cell1.mat', '2_23_14_5_data.mat', '/media/shababo/data/0221232018_zdata/s2c1-pre - 8_C0'
%              '2_23_slice3_cell2.mat', '2_23_14_28_data.mat', '/media/shababo/data/0221232018_zdata/s2c1-pre - 9_C0'
%              '2_23_slice3_cell3.mat', '2_23_14_40_data.mat', '/media/shababo/data/0221232018_zdata/s2c1-pre - 10_C0'}; 
%          
%          
% spike_trials = {1,1,1,1,1,1,1,1,1,1,1};
% current_trials = {4,4,4,4,4,4,4,4,4,4,4};
% z_center = [30 30 200 200 200 200 200 200 200 200 200];


filenames = {'7_17_slice1_cell1.mat', '7_17_17_59_data.mat' % good
             '7_17_slice1_cell2.mat', '7_17_18_12_data.mat' % good
             '7_17_slice1_cell3.mat', '7_17_18_28_data.mat' % good
             '7_17_slice2_cell1.mat', '7_17_18_50_data.mat' % crappy cell - EXCLUDE
             '7_17_slice2_cell3.mat', '7_17_19_13_data.mat' % good
             '7_17_slice1_cell4.mat', '7_17_19_32_data.mat' % good
             } 
         
 spike_trials = {NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN};
current_trials = {4, 4, 4, 3, 4, 4, 2, 6, 4, 3, 3};
z_center = [70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70];
         
% filenames = {'7_13_slice1_cell4.mat', '7_13_16_37_data.mat'
%              '7_13_slice1_cell4.mat', '7_13_16_37_data.mat'
%              '7_13_slice1_cell4.mat', '7_13_16_37_data.mat'
%              '7_13_slice1_cell4.mat', '7_13_16_37_data.mat'
%              '7_13_slice1_cell4.mat', '7_13_16_37_data.mat'} 
%          
% spike_trials = {NaN, NaN, NaN, NaN, NaN};
% current_trials = {4, 5, 6, 7, 8};
% z_center = [70, 70, 70, 70, 70];

colors = jet(length(filenames));
% colors(:) = 0;
% colors(1:size(filenames,1)) = 1; %red
% figure_for_cosyne = figure

%%

% clear result_x
% z_pos_vs_cur_and_spike_time = figure;
% new_fig = 0;
% if ~exist('curr_vs_time','var')
%     new_fig = 1;
%     curr_vs_time = figure;
% end
cell_choice = setdiff(1:size(filenames,1),[]);
% cell_choice = [5];
for jj = 1:length(cell_choice)%setsetdiff(1:size(filenames,1),6)
    jj
    j = cell_choice(jj)
%     clear result_x(result_ind)
    load(filenames{j,2});
    load(filenames{j,1});
    experiment_setup = exp_data.experiment_setup;
    
    result_ind = jj;
    
    do_detect = 0;


%     if ~isempty(filenames{j,3})
%         disp('...')
% %         load(filenames{j,2});
%         [nuclear_locs,fluor_vals,nuclear_locs_image_coord] = detect_nuclei(filenames{j,3},[],[],[],do_detect,[],0);
%         offsets = nuclear_locs - [experiment_setup.center_pos_um(1:2) 30];
% 
%         [targ_error, index] = min(sqrt(sum(offsets.^2,2)));
%         result_x(result_ind).fluor_val = fluor_vals(index);
%         result_x(result_ind).cell_pos = nuclear_locs(index,:);
%         result_x(result_ind).exp_cell_pos = [experiment_setup.center_pos_um(1:2) 30];
%         result_x(result_ind).err_cell_pos = result_x(result_ind).exp_cell_pos - result_x(result_ind).cell_pos;
%         result_x(result_ind).err_cell_pos = result_x(result_ind).exp_cell_pos - result_x(result_ind).cell_pos;
%     else
%         result_x(result_ind).fluor_val = NaN;
%         result_x(result_ind).cell_pos = [experiment_setup.center_pos_um(1:2) experiment_setup.piezo_center];

        result_x(result_ind).fluor_val = experiment_setup.fluor_vals(experiment_setup.select_cell_index_full);
        result_x(result_ind).cell_pos = [experiment_setup.patched_cell_loc];

%     end
    
    
    
    if ~isnan(current_trials{j})
        disp('wahat')
        [result_x(result_ind).current_traces, ~,~, result_x(result_ind).current_stim_traces] = get_traces(data,current_trials{j});
        target_pos = data.trial_metadata(current_trials{j}).stim_key([data.trial_metadata(current_trials{j}).sequence.precomputed_target_index],:);
        result_x(result_ind).current_x_pos = target_pos(:,1) - result_x(result_ind).cell_pos(1);
        result_x(result_ind).max_curr = result_x(result_ind).current_traces(:,1) - min(result_x(result_ind).current_traces,[],2);
    end 
    
%     assignin('base','result_x',result_x)
%     return
    
%     if ~isnan(spike_trials{j}) && ~isnan(current_trials{j})
%         tested_pos = unique(result_x(result_ind).spike_z_pos);
%         if ~isempty(setdiff(tested_pos,unique(result_x(result_ind).current_z_pos)))
%             disp([filenames{j,1} ' spike and current pos sets differ'])
%             continue
%         end
% 
% 
%         for i = 1:length(tested_pos)
%             these_trials = result_x(result_ind).spike_z_pos == tested_pos(i);
%             result_x(result_ind).spike_time_means(i) = nanmean(result_x(result_ind).spike_times(these_trials));
%             result_x(result_ind).spike_time_jitter(i) = nanstd(result_x(result_ind).spike_times(these_trials));
%             these_trials = result_x(result_ind).current_z_pos == tested_pos(i);
%             result_x(result_ind).max_curr_means(i) = nanmean(result_x(result_ind).max_curr(these_trials));
%         end
% 
% %         figure(curr_vs_time)
% %         subplot(121)
% %         [ordered_curr, curr_order] = sort(result_x(result_ind).max_curr_means);
% %         plot(result_x(result_ind).max_curr_means(curr_order),result_x(result_ind).spike_time_means(curr_order)/20,'-k','LineWidth',2);
% %         hold on
% %         if new_fig
% %             xlabel('Peak Current Mean (pA)')
% %             ylabel('Spike Time Mean (msec)')
% %             title('Peak Current vs. Spike Time')
% %         end
% %         subplot(122)
% %         semilogy(result_x(result_ind).max_curr_means(curr_order),result_x(result_ind).spike_time_jitter(curr_order)/20,'-k','LineWidth',2);
% %         hold on
% %         if new_fig
% %             xlabel('Peak Current Mean (pA)')
% %             ylabel('Spike Time Std. Dev. (msec)')
% %             title('Peak Current vs. Spike Time Jitter')
% %         end
%     end
%     assignin('base','result_x',result_x)
%     
%     
%     figure(figure_for_cosyne)
%     if ~isnan(spike_trials{j})
%         subplot(1,length(cell_choice),jj)
%         yyaxis left
% %         scatter(result_x(result_ind).spike_z_pos,result_x(result_ind).spike_times/20,[],colors(j,:));
%         hold on
%         plot(tested_pos,result_x(result_ind).spike_time_means/20,'color',colors(2,:))
%         xlabel('z distance (um)')
%         ylabel('spike time (msec)')
%         title('Z Distance vs. Spike Times')
%         
%     end
%     if ~isnan(current_trials{j})
%         subplot(1,length(cell_choice),jj)
%         yyaxis right
% %         scatter(result_x(result_ind).current_z_pos,result_x(result_ind).max_curr,[],colors(j,:));
%         hold on
%         plot(tested_pos,result_x(result_ind).max_curr_means,'color',colors(1,:))
%         ylim([0 2500])
%         xlabel('z distance (um)')
%         ylabel('peak current (pA)')
%         title('Z Distance vs. Peak Current')
%     end
    
end

%%
figure
for jj = 1:length(result_x)
    
    scatter(result_x(jj).current_x_pos,result_x(jj).max_curr,10,colors(jj,:))
    hold on
    
end

xlabel('x (vertical) dist from nuc detect location (um)')
ylabel('current (pA)')
title('6 cells, x (vertical) modulation current response')

%%
figure
for jj = 1:length(result_x)
    
    scatter(result_x(jj).fluor_val,max(result_x(jj).max_curr),20,colors(mod(jj,size(colors,1))+1,:),'linewidth',4)
    hold on
    
end

for jj = 1:length(result_z)
    
    scatter(result_z(jj).fluor_val,max(result_z(jj).max_curr),20,colors(mod(jj,size(colors,1))+1,:),'linewidth',4)
    hold on
    
end

xlim([0 700])
%%
figure
cells_to_plot = [3:1:11];
[ha, pos] = tight_subplot(1,length(cells_to_plot),.01,.25,.1); 

for jj = 1:length(cells_to_plot)
    j = cells_to_plot(jj);
    axes(ha(jj))
    yyaxis left
    tested_pos = unique([result_x(result_ind).spike_z_pos]);
%     these_trials = result_x(result_ind).current_targ_pos(:,1) == tested_pos_y(6);
    %     scatter(result_x(result_ind).current_targ_pos(these_trials,2),result_x(result_ind).max_curr(these_trials),[],colors(j,:));
    hold on
%     [~,this_zero_pos] = min(abs(tested_pos_x));
%     scaling = result_x(result_ind).x_max_curr_means(this_zero_pos);
%     [~,this_zero_pos] = min(abs(result_x(result_ind).current_targ_pos(these_trials,2)));
%     these_powers = result_x(result_ind).spatial_adj_power(these_trials);
    %     scaling = these_powers(this_zero_pos)*gain_mle(28+j)*1000;
    %     plot(-20:20,scaling*shape_template(sub2ind(size(shape_template),36*ones(size(-20:20)),(-20:20)+36)),'color',colors(j,:),'linewidth',1);
    plot(tested_pos,result_x(result_ind).max_curr_means,'-x','color',colors(1,:),'linewidth',2)
    ylim([0 1600])
    xlim([-60 60])
    xlabel('Z Distance (um)')
    if jj == 1
        ylabel('Peak Current (pA)')
    else
%         set(ha(jj),'YTickLabel','')
    end
    title(['Cell ' num2str(j)])
    
    
%     subplot(2,size(filenames,1),j+size(filenames,1))
    yyaxis right
    tested_pos = unique([result_x(result_ind).current_z_pos]);
    %     scatter(result_x(result_ind).current_targ_pos(these_trials,1),result_x(result_ind).max_curr(these_trials),[],colors(j,:));
    hold on
%     [~,this_zero_pos] = min(abs(tested_pos_y))
%     scaling = result_x(result_ind).y_max_curr_means(this_zero_pos);
%     [~,this_zero_pos] = min(abs(result_x(result_ind).current_targ_pos(these_trials,1)));
%     these_powers = result_x(result_ind).spatial_adj_power(these_trials);
    %     scaling = these_powers(this_zero_pos)*gain_mle(28+j)*1000;
    %     plot(-20:20,scaling*shape_template(sub2ind(size(shape_template),(-20:20)+36,36*ones(size(-20:20)))),'color',colors(j,:),'linewidth',1)
    plot(tested_pos,result_x(result_ind).spike_time_means/20,'-x','color',colors(4,:),'linewidth',2)
    ylim([0 15])
    xlim([-60 60])
%     xlabel('Vertical Distance (um)')
    if jj == length(cells_to_plot)
        ylabel('Spike Time (msec)')
    else
%         set(ha(jj),'YTickLabel','')
    end
%     title('Vertical Distance vs. Peak Current')
end

%% detect nucs
do_detect = 0;
for jj = 1:length(cell_choice)
    
    j = cell_choice(jj)
    if ~isempty(filenames{j,3})
        disp('...')
        load(filenames{j,2});
        [nuclear_locs,fluor_vals,nuclear_locs_image_coord] = detect_nuclei(filenames{j,3},[],[],[],do_detect,[],0);
        offsets = nuclear_locs - [experiment_setup.center_pos_um(1:2) 30];

        [targ_error, index] = min(sqrt(sum(offsets.^2,2)));
        result_x(result_ind).fluor_val = fluor_vals(index);
        result_x(result_ind).cell_pos = nuclear_locs(index,:);
        result_x(result_ind).exp_cell_pos = [experiment_setup.center_pos_um(1:2) 30];
        result_x(result_ind).err_cell_pos = result_x(result_ind).exp_cell_pos - result_x(result_ind).cell_pos;
        result_x(result_ind).err_cell_pos = result_x(result_ind).exp_cell_pos - result_x(result_ind).cell_pos;
    else
        result_x(result_ind).fluor_val = NaN;
        result_x(result_ind).cell_pos = [experiment_setup.center_pos_um(1:2) experiment_setup.piezo_center];

    end

end