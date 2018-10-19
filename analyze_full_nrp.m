%%

filenames = {'9_27_slice3_cell1_2.mat', '9_27_16_0_data.mat'
             '9_28_slice1_cell1_2.mat', '9_28_13_50_data.mat'
             '9_28_slice1_cell3_4.mat', '9_28_14_8_data.mat'
             '9_28_slice1_cell5_6.mat', '9_28_14_21_data.mat'
             '9_28_slice2_cell3_4.mat', '9_28_15_19_data.mat'
             '9_28_slice2_cell5_6.mat', '9_28_15_34_data.mat'
             };      
                                                    
ch1_cell = [1,1,1,1,1,1];
ch2_cell = [1,1,1,1,1,1];
                                                    
spike_trials = {2:4,2:7,2:4,2:4,2:4,3:5};
current_trials = {[]};


%%

% clear result_full_nrp

for j = 6:size(filenames,1)
    
    load(filenames{j,2});
    load(filenames{j,1}); 
    experiment_setup = exp_data.experiment_setup;
    
    if ~isempty(spike_trials{j})
        
        [result_full_nrp(j).spike_traces_c1, result_full_nrp(j).spike_traces_c2, this_seq, result_full_nrp(j).spike_stim_traces, full_stim_key] = get_traces(data,spike_trials{j});
        result_full_nrp(j).spike_traces_c2 = result_full_nrp(j).spike_traces_c2*1000;
        result_full_nrp(j).spike_targ_pos_c1 = bsxfun(@minus,full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc);
        result_full_nrp(j).c1_pos = experiment_setup.patched_cell_loc;
        result_full_nrp(j).c2_pos = experiment_setup.patched_cell_loc_2;
        result_full_nrp(j).spike_targ_pos_c2 = bsxfun(@minus,full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc_2);
        result_full_nrp(j).spike_targ_pos = full_stim_key([this_seq.precomputed_target_index],:);
        result_full_nrp(j).spike_targ_power = [this_seq.target_power];
    
        cell_spike_times_c1 = detect_peaks(-bsxfun(@minus,result_full_nrp(j).spike_traces_c1,result_full_nrp(j).spike_traces_c1(:,1)),20,30,0,Inf,-Inf,0,0,1);
        cell_spike_times_c2 = detect_peaks(-bsxfun(@minus,result_full_nrp(j).spike_traces_c2,result_full_nrp(j).spike_traces_c2(:,1)),20,30,0,Inf,-Inf,0,0,1);
        result_full_nrp(j).spike_times_c1 = zeros(size(cell_spike_times_c1));
        result_full_nrp(j).spike_times_c2 = zeros(size(cell_spike_times_c2));

        for i = 1:length(cell_spike_times_c1)
            if ~isempty(cell_spike_times_c1{i})
                result_full_nrp(j).spike_times_c1(i) = cell_spike_times_c1{i};
            else
                result_full_nrp(j).spike_times_c1(i) = NaN;
            end
        end
        for i = 1:length(cell_spike_times_c2)
            if ~isempty(cell_spike_times_c2{i})
                result_full_nrp(j).spike_times_c2(i) = cell_spike_times_c2{i};
            else
                result_full_nrp(j).spike_times_c2(i) = NaN;
            end
        end
    end
    
end
%%
j = 5
offset = .05;
unique_powers = unique(result_full_nrp(j).spike_targ_power);

figure; 
for i = 1:length(unique_powers)
    these_trials = result_full_nrp(j).spike_targ_power == unique_powers(i);
    these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c1),:);
    subplot(1,length(unique_powers),i)

    scatter3(result_full_nrp(j).spike_targ_pos(:,1),result_full_nrp(j).spike_targ_pos(:,2),result_full_nrp(j).spike_targ_pos(:,3),5,'k','filled')
    hold on
    scatter3(these_locs(:,1),these_locs(:,2),these_locs(:,3),[],[1 0 0],'filled','MarkerFaceAlpha',.5)
    hold on
    scatter3(result_full_nrp(j).c1_pos(1),result_full_nrp(j).c1_pos(2),result_full_nrp(j).c1_pos(3),[],[1 0 0])
    
    these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c2),:);
    scatter3(these_locs(:,1)+.25,these_locs(:,2)+offset,these_locs(:,3)+offset,[],[0 0 1],'filled','MarkerFaceAlpha',.5)
    hold on
    scatter3(result_full_nrp(j).c2_pos(1)+offset,result_full_nrp(j).c2_pos(2)+offset,result_full_nrp(j).c2_pos(3),[],[0 0 1])
    xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
    ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
    zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
end
% subplot(1,length(unique_powers)+1,length(unique_powers)+1)
% scatter3(result_full_nrp(j).spike_targ_pos(:,1),result_full_nrp(j).spike_targ_pos(:,2),result_full_nrp(j).spike_targ_pos(:,3),'filled')

% figure; 
% for i = 1:length(unique_powers)
%     these_trials = result_full_nrp(j).spike_targ_power == unique_powers(i);
%     these_locs = result_full_nrp.spike_targ_pos(these_trials,:);
%     subplot(1,length(unique_powers)+1,i)
%     scatter3(these_locs(:,1),these_locs(:,2),these_locs(:,3),[],(500-result_full_nrp.spike_times_c2(these_trials))/500,'filled','MarkerFaceAlpha',.33)
%     hold on
%     scatter3(result_full_nrp(j).c2_pos(1),result_full_nrp(j).c2_pos(2),result_full_nrp(j).c2_pos(3))
% end
% subplot(1,length(unique_powers)+1,length(unique_powers)+1)
% scatter3(these_locs(:,1),these_locs(:,2),these_locs(:,3),'filled')