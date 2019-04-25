
%%
realign = 1;

    
         

if ~realign  
    
    filenames_nrp = {'9_27_slice3_cell1_2.mat', '9_27_16_0_data.mat'
             '9_28_slice1_cell1_2.mat', '9_28_13_50_data.mat'
             '9_28_slice1_cell3_4.mat', '9_28_14_8_data.mat'
             '9_28_slice1_cell5_6.mat', '9_28_14_21_data.mat'
             '9_28_slice2_cell3_4.mat', '9_28_15_19_data.mat'
             '9_28_slice2_cell5_6.mat', '9_28_15_34_data.mat'
             '10_13_slice1_cell5_6.mat', '10_13_15_0_data.mat'
             '10_13_slice1_cell1_2.mat', '10_13_14_21_data.mat'
             '10_16_slice1_cell3_4.mat', '10_16_15_14_data.mat'
             '10_23_slice1_cell1_2.mat', '10_23_17_29_data.mat'
             '10_23_slice1_cell3_4.mat', '10_23_17_47_data.mat'
%              '10_24_slice1_cell1_2.mat', '10_24_16_5_data.mat'
             '10_24_slice1_cell3_4.mat', '10_24_16_21_data.mat'
             '10_24_slice2_cell1_2.mat', '10_24_16_37_data.mat'
             };  
         
    ch1_cell = [1,1,1,1,1,1,1,1,1,1,1,1,0];
    ch2_cell = [1,1,1,1,1,1,0,0,0,0,1,0,1];

    spike_trials = {2:4,2:7,2:4,2:4,2:4,3:5,[1 3:5],[1 3:5],[1 3:5],5:7,2:5,2:5,3:5};
    current_trials = {[]};
    
    set_color = {'r','r','r','r','r','r','b','b','b','b','b','b','b'};
    
    spike_thresh_ch1 = 30*ones(length(ch1_cell),1);
    spike_thresh_ch2 = 30*ones(length(ch2_cell),1);
    do_hp = 0;
    
else
    
    filenames_nrp = {'1_14_slice1_cell2next.mat', '1_14_15_10_data.mat'
                     '1_15_slice1_cell2.mat', '1_15_15_48_data.mat'
                     '1_15_slice2_cell3_4.mat', '1_15_16_26_data.mat'
                     '1_15_slice2_cell5_6.mat', '1_15_16_40_data.mat'
                     '4_19_slice1_cell1.mat', '4_19_14_25_data.mat'};
                      
    ch1_cell = [1 1 1 0 1];
    ch2_cell = [0 0 1 1 0];
    
    spike_trials = {2:3, 2:5, 3:6, 2:5, 3:6}; %
    
    
    set_color = {'r','r','r','r','r','r','b','b','b','b','b','b','b'};
    
    spike_thresh_ch1 = [8 20 30 12 20];
    spike_thresh_ch2 = [8 20 25 12 20];
    
    do_hp_ch1 = [1 0 0 1 0];
    do_hp_ch2 = [1 0 0 1 0];
    
    total_cells = sum(ch1_cell) + sum(ch2_cell);
end

[-2 -6 0
 -3 3 -30];
%%

% clear result_full_nrp

for j = 1:size(filenames_nrp,1)
    
    j
    
    load(filenames_nrp{j,2});
    load(filenames_nrp{j,1}); 
    experiment_setup = exp_data.experiment_setup;
    
    if ~isempty(spike_trials{j})
        result_full_nrp(j).nuclear_locs = experiment_setup.nuclear_locs;
        [result_full_nrp(j).spike_traces_c1, result_full_nrp(j).spike_traces_c2, this_seq, result_full_nrp(j).spike_stim_traces, full_stim_key] = get_traces(data,spike_trials{j});
        if ch1_cell(j)
            disp('cell 1')
            
        
            result_full_nrp(j).spike_targ_pos_c1 = bsxfun(@minus,full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc);
            result_full_nrp(j).c1_pos = experiment_setup.patched_cell_loc;
            result_full_nrp(j).spike_targ_pos = full_stim_key([this_seq.precomputed_target_index],:);
            result_full_nrp(j).spike_targ_power = [this_seq.target_power];
    
            cell_spike_times_c1 = detect_peaks(-bsxfun(@minus,result_full_nrp(j).spike_traces_c1,...
                median(result_full_nrp(j).spike_traces_c1,2)),...
                spike_thresh_ch1(j),30,0,Inf,-Inf,0,do_hp_ch1(j),1);
            result_full_nrp(j).spike_times_c1 = zeros(size(cell_spike_times_c1));
            for i = 1:length(cell_spike_times_c1)
                if ~isempty(cell_spike_times_c1{i})
                    result_full_nrp(j).spike_times_c1(i) = cell_spike_times_c1{i};
                else
                    result_full_nrp(j).spike_times_c1(i) = NaN;
                end
            end
        end
        if ch2_cell(j)
            disp('cell 2')
%             result_full_nrp(j).spike_traces_c2 = result_full_nrp(j).spike_traces_c2;
            result_full_nrp(j).c2_pos = experiment_setup.patched_cell_loc_2;
            result_full_nrp(j).spike_targ_pos_c2 = bsxfun(@minus,full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc_2);
            result_full_nrp(j).spike_targ_pos = full_stim_key([this_seq.precomputed_target_index],:);
            result_full_nrp(j).spike_targ_power = [this_seq.target_power];
            cell_spike_times_c2 = detect_peaks(-bsxfun(@minus,result_full_nrp(j).spike_traces_c2,result_full_nrp(j).spike_traces_c2(:,1)),...
                spike_thresh_ch2(j),30,0,Inf,-Inf,0,do_hp_ch2(j),1);
            result_full_nrp(j).spike_times_c2 = zeros(size(cell_spike_times_c2));
            for i = 1:length(cell_spike_times_c2)
                if ~isempty(cell_spike_times_c2{i})
                    result_full_nrp(j).spike_times_c2(i) = cell_spike_times_c2{i};
                else
                    result_full_nrp(j).spike_times_c2(i) = NaN;
                end
            end
        end
        
    end
    
end

%% plot spiking locs

offset = .00;

spike_time_max = 100;

views = [0 90; 0 0; 90 0];
colors = jet(total_cells);
count = 1;
for j = 1:size(filenames_nrp,1)
    
    figure
    unique_powers = unique(result_full_nrp(j).spike_targ_power);%[15 35 55]; %
    
    
    
    for k = 1:3
        for i = 1:length(unique_powers)
            
            subplot(3,length(unique_powers),i + length(unique_powers)*(k-1))
            these_trials = result_full_nrp(j).spike_targ_power == unique_powers(i);
            scatter3(result_full_nrp(j).spike_targ_pos(these_trials,1),result_full_nrp(j).spike_targ_pos(these_trials,2),-result_full_nrp(j).spike_targ_pos(these_trials,3),5,'k','filled')
            hold on
            scatter3(result_full_nrp(j).nuclear_locs(:,1),result_full_nrp(j).nuclear_locs(:,2),-result_full_nrp(j).nuclear_locs(:,3),10,'r','filled')
            hold on
            if ch1_cell(j)
                
                these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c1) &...
                    result_full_nrp(j).spike_times_c1 < spike_time_max,:);                
                scatter3(these_locs(:,1),these_locs(:,2),-these_locs(:,3),[],colors(count,:),'filled','MarkerFaceAlpha',.33)
                hold on
                scatter3(result_full_nrp(j).c1_pos(1),result_full_nrp(j).c1_pos(2),-result_full_nrp(j).c1_pos(3),60,colors(count,:))
            end
            if ch2_cell(j)
                these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c2) &...
                    result_full_nrp(j).spike_times_c2 < spike_time_max,:);
                scatter3(these_locs(:,1)+offset,these_locs(:,2)+offset,-these_locs(:,3)+offset,[],colors(count+ch1_cell(j),:),'filled','MarkerFaceAlpha',.33)
                hold on
                scatter3(result_full_nrp(j).c2_pos(1)+offset,result_full_nrp(j).c2_pos(2)+offset,-result_full_nrp(j).c2_pos(3),60,colors(count+ch1_cell(j),:))


            end
            
            
            xlabel('vertical')
            ylabel('horizontal')
            zlabel('axial/horizontal')
            title(sprintf('Pair %d, Power: %d',j,unique_powers(i)))
            view(views(k,1),views(k,2))
            axis image
            xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
            ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
            zlim([min(-result_full_nrp(j).spike_targ_pos(:,3)) max(-result_full_nrp(j).spike_targ_pos(:,3))])
        end
%         subplot(3,length(unique_powers)+1,length(unique_powers)+1+(length(unique_powers)+1)*(k-1))
%         scatter3(result_full_nrp(j).spike_targ_pos(:,1),result_full_nrp(j).spike_targ_pos(:,2),result_full_nrp(j).spike_targ_pos(:,3),'filled')
%         xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
%         ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
%         zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
%         xlabel('vertical')
%         ylabel('horizontal')
%         zlabel('axial/horizontal')
%         title(sprintf('Pair %d, All Targets',j))
%         view(views(k,1),views(k,2))
        
    end
    count = count + ch1_cell(j) + ch2_cell(j);
end
%% plot spiking locs with overlap

offset = .05;

spike_time_max = 200;

views = [0 90; 0 0; 90 0];

for j = 10
    unique_powers = unique(result_full_nrp(j).spike_targ_power);
    figure; 
    for k = 1:3
        for i = 1:length(unique_powers)
            
            subplot(3,length(unique_powers),i + length(unique_powers)*(k-1))
            
            if ch1_cell(j) && ch2_cell(j)
                these_trials = result_full_nrp(j).spike_targ_power == unique_powers(i);
                these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c1) &...
                    result_full_nrp(j).spike_times_c1 < spike_time_max & (isnan(result_full_nrp(j).spike_times_c2) | result_full_nrp(j).spike_times_c2 > spike_time_max),:);

                scatter3(result_full_nrp(j).spike_targ_pos(:,1),result_full_nrp(j).spike_targ_pos(:,2),result_full_nrp(j).spike_targ_pos(:,3),5,'k','filled')
                hold on
                scatter3(these_locs(:,1),these_locs(:,2),these_locs(:,3),[],[1 0 0],'filled','MarkerFaceAlpha',.33)
                hold on
                scatter3(result_full_nrp(j).c1_pos(1),result_full_nrp(j).c1_pos(2),result_full_nrp(j).c1_pos(3),60,[1 0 0])
%             end
%             if ch2_cell(j)
                these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c2) &...
                    result_full_nrp(j).spike_times_c2 < spike_time_max & (isnan(result_full_nrp(j).spike_times_c1) | result_full_nrp(j).spike_times_c1 > spike_time_max),:);
                scatter3(these_locs(:,1)+offset,these_locs(:,2)+offset,these_locs(:,3)+offset,[],[0 0 1],'filled','MarkerFaceAlpha',.33)
                hold on
                scatter3(result_full_nrp(j).c2_pos(1)+offset,result_full_nrp(j).c2_pos(2)+offset,result_full_nrp(j).c2_pos(3),60,[0 0 1])

                these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c2) &...
                    result_full_nrp(j).spike_times_c2 < spike_time_max & ~isnan(result_full_nrp(j).spike_times_c1) & result_full_nrp(j).spike_times_c1 < spike_time_max,:);
                scatter3(these_locs(:,1)+offset,these_locs(:,2)+offset,these_locs(:,3)+offset,[],[0 1 0],'filled','MarkerFaceAlpha',.33)
            end
            
            xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
            ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
            zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
            xlabel('vertical')
            ylabel('horizontal')
            zlabel('axial/horizontal')
            title(sprintf('Pair %d, Power: %d',j,unique_powers(i)))
            view(views(k,1),views(k,2))
            axis image
        end
%         subplot(3,length(unique_powers)+1,length(unique_powers)+1+(length(unique_powers)+1)*(k-1))
%         scatter3(result_full_nrp(j).spike_targ_pos(:,1),result_full_nrp(j).spike_targ_pos(:,2),result_full_nrp(j).spike_targ_pos(:,3),'filled')
%         xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
%         ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
%         zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
%         xlabel('vertical')
%         ylabel('horizontal')
%         zlabel('axial/horizontal')
%         title(sprintf('Pair %d, All Targets',j))
%         view(views(k,1),views(k,2))
    end
end


%% plot first spiker

offset = .05;
spike_time_max = 200;

views = [0 90; 0 0; 90 0];

for j = 1:6
    unique_powers = unique(result_full_nrp(j).spike_targ_power);
    figure; 
    for k = 1:3
        for i = 1:length(unique_powers)
            
            subplot(3,length(unique_powers)+1,i + (length(unique_powers)+1)*(k-1))
            
            these_trials = result_full_nrp(j).spike_targ_power' == unique_powers(i) & ~isnan(result_full_nrp(j).spike_times_c1) & ...
                result_full_nrp(j).spike_times_c1 < spike_time_max & (result_full_nrp(j).spike_times_c1 < result_full_nrp(j).spike_times_c2 | isnan(result_full_nrp(j).spike_times_c2));
            these_locs = result_full_nrp(j).spike_targ_pos(these_trials',:);
            
            scatter3(these_locs(:,1),these_locs(:,2),these_locs(:,3),40,[1 0 0],'filled','MarkerFaceAlpha',.33)
            hold on
            scatter3(result_full_nrp(j).c1_pos(1),result_full_nrp(j).c1_pos(2),result_full_nrp(j).c1_pos(3),50,[1 0 0])
            these_trials = result_full_nrp(j).spike_targ_power' == unique_powers(i) & ~isnan(result_full_nrp(j).spike_times_c2) & ...
                result_full_nrp(j).spike_times_c2 < spike_time_max & (result_full_nrp(j).spike_times_c2 < result_full_nrp(j).spike_times_c1 | isnan(result_full_nrp(j).spike_times_c1));
            these_locs = result_full_nrp(j).spike_targ_pos(these_trials',:);
            scatter3(these_locs(:,1)+.25,these_locs(:,2)+offset,these_locs(:,3)+offset,[],[0 0 1],'filled','MarkerFaceAlpha',.33)
            hold on
            scatter3(result_full_nrp(j).c2_pos(1)+offset,result_full_nrp(j).c2_pos(2)+offset,result_full_nrp(j).c2_pos(3),50,[0 0 1])
            xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
            ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
            zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
            xlabel('vertical')
            ylabel('horizontal')
            zlabel('axial/horizontal')
            title(sprintf('Pair %d, Power: %d',j,unique_powers(i)))
            view(views(k,1),views(k,2))
        end
        subplot(3,length(unique_powers)+1,length(unique_powers)+1+(length(unique_powers)+1)*(k-1))
        scatter3(result_full_nrp(j).spike_targ_pos(:,1),result_full_nrp(j).spike_targ_pos(:,2),result_full_nrp(j).spike_targ_pos(:,3),'filled')
        xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
        ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
        zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
        xlabel('vertical')
        ylabel('horizontal')
        zlabel('axial/horizontal')
        title(sprintf('Pair %d, All Targets',j))
        view(views(k,1),views(k,2))
    end
end

%% plot spike time diffs

spike_time_max = 200;



for j = 1:6
    figure
    unique_powers = unique(result_full_nrp(j).spike_targ_power);
        for i = 1:length(unique_powers)
            
            subplot(3,length(unique_powers),i + (length(unique_powers))*(1-1))
            
            these_trials = result_full_nrp(j).spike_targ_power' == unique_powers(i) & ~isnan(result_full_nrp(j).spike_times_c1) & ...
                result_full_nrp(j).spike_times_c1 < spike_time_max;
            histogram(result_full_nrp(j).spike_times_c1(these_trials),0:4:200)
            title(sprintf('Spike Times: Pair %d, Cell 1, power: %d',j,unique_powers(i)))
            xlim([0 200])
            ylim([0 10])
            
            subplot(3,length(unique_powers),i + (length(unique_powers))*(2-1))
            these_trials = result_full_nrp(j).spike_targ_power' == unique_powers(i) & ~isnan(result_full_nrp(j).spike_times_c2) & ...
                result_full_nrp(j).spike_times_c2 < spike_time_max;
            histogram(result_full_nrp(j).spike_times_c2(these_trials),0:4:200)
            title(sprintf('Spike Times: Pair %d, Cell 2, power: %d',j,unique_powers(i)))
            xlim([0 200])
            ylim([0 10])
            
            subplot(3,length(unique_powers),i + (length(unique_powers))*(3-1))
            these_trials = result_full_nrp(j).spike_targ_power' == unique_powers(i) & ~isnan(result_full_nrp(j).spike_times_c2) & ...
                ~isnan(result_full_nrp(j).spike_times_c1) & result_full_nrp(j).spike_times_c1 < spike_time_max & result_full_nrp(j).spike_times_c2 < spike_time_max;
            histogram(result_full_nrp(j).spike_times_c1(these_trials) - result_full_nrp(j).spike_times_c2(these_trials),-200:8:200)
            xlim([-200 200])
            title(sprintf('Spike Time Diffs (c1 - c2): Pair %d power: %d',j,unique_powers(i)))
            ylim([0 10])
        end

end


%% plot spiking locs centered

offset = .05;

spike_time_max = 140;

views = [0 90; 0 0; 90 0];
    figure; 

xy_bound = 8;
clear max_dist
clear max_dist_loc
z_test_locs = cell(5,1);
z_spike_locs = cell(5,1);
for k = 1:3
    
    unique_powers = unique(result_full_nrp(1).spike_targ_power);
    for i = 1:length(unique_powers)
        subplot(3,length(unique_powers),i + length(unique_powers)*(k-1))
        for j = 1:4
        
            unique_powers = unique(result_full_nrp(j).spike_targ_power);
            
            
            these_trials = result_full_nrp(j).spike_targ_power == unique_powers(i);
%             these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c1) &...
%                 result_full_nrp(j).spike_times_c1 < spike_time_max,:);
            
            if ch1_cell(j)
                c1_targs = result_full_nrp(j).spike_targ_pos_c1(these_trials,:);
%                 z_aligned = c1_targs(:,1) > -xy_bound & c1_targs(:,1) < xy_bound & c1_targs(:,2) > -xy_bound & c1_targs(:,2) < xy_bound;
                these_locs_cell = c1_targs(:,:);
                scatter3(c1_targs(:,1),c1_targs(:,2),-c1_targs(:,3),15,'k','filled','MarkerFaceAlpha',1./20.)

%                 z_test_locs{i} = [z_test_locs{i}; abs(these_locs_cell(:,3))];
                hold on
            end
%             scatter3(these_locs(:,1)-result_full_nrp(j).c1_pos(1),these_locs(:,2)-result_full_nrp(j).c1_pos(2),these_locs(:,3)-result_full_nrp(j).c1_pos(3),[],[0 0 1],'filled','MarkerFaceAlpha',.33)
%             hold on
%             scatter3(result_full_nrp(j).c1_pos(1),result_full_nrp(j).c1_pos(2),result_full_nrp(j).c1_pos(3),60,[1 0 0])
            
            if ch2_cell(j)
                c2_targs = result_full_nrp(j).spike_targ_pos_c2(these_trials,:);
%                 scatter3(c2_targs(:,1),c2_targs(:,2),-c2_targs(:,3),15,'k','filled','MarkerFaceAlpha',1./20.)
                z_aligned = c2_targs(:,1) > -xy_bound & c2_targs(:,1) < xy_bound & c2_targs(:,2) > -xy_bound & c2_targs(:,2) < xy_bound;
                these_locs_cell = c2_targs(z_aligned,:);
                z_test_locs{i} = [z_test_locs{i}; abs(these_locs_cell(:,3))];
                hold on
            end
%             these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c2) &...
%                 result_full_nrp(j).spike_times_c2 < spike_time_max,:);
%             scatter3(these_locs(:,1)-result_full_nrp(j).c2_pos(1),these_locs(:,2)-result_full_nrp(j).c2_pos(2),these_locs(:,3)-result_full_nrp(j).c2_pos(3),[],[0 0 1],'filled','MarkerFaceAlpha',.33) 
%             
%             xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
%             ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
%             zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
%             xlabel('vertical')
%             ylabel('horizontal')
%             zlabel('axial/horizontal')
%             title(sprintf('Pair %d, Power: %d',j,unique_powers(i)))
%             view(views(k,1),views(k,2))
        end
        subplot(3,length(unique_powers),i + length(unique_powers)*(k-1))
        cell_count = 1
        for j = 1:4
%             subplot(3,length(unique_powers),i + length(unique_powers)*(k-1))
            unique_powers = unique(result_full_nrp(j).spike_targ_power);
            these_trials = result_full_nrp(j).spike_targ_power == unique_powers(i);
            if ~any(these_trials)
                bad_pow(j,i) = 1;
            end
            if ch1_cell(j)
                these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c1) &...
                    result_full_nrp(j).spike_times_c1 < spike_time_max,:);
                these_locs_cell = result_full_nrp(j).spike_targ_pos_c1(these_trials' & ~isnan(result_full_nrp(j).spike_times_c1) &...
                    result_full_nrp(j).spike_times_c1 < spike_time_max,:);
%                 these_locs_dist = vecnorm(these_locs_cell');
%                 z_aligned = these_locs_cell(:,1) > -xy_bound & these_locs_cell(:,1) < xy_bound & these_locs_cell(:,2) > -xy_bound & these_locs_cell(:,2) < xy_bound;
%                 these_locs_cell = these_locs_cell(z_aligned,:);
                these_locs_dist = sqrt(sum(these_locs_cell.^2,2));
%                 if ~isempty(these_locs_dist)
%                     if k == 1
%                         [max_dist(i,cell_count),max_i] = max(these_locs_dist);
%                         max_dist_loc(i,cell_count,:) = these_locs_cell(max_i,:);
%                     end
% %                     scatter3(max_dist_loc(i,cell_count,1),max_dist_loc(i,cell_count,2),-max_dist_loc(i,cell_count,3),[],[0 0 1],'filled')
%                 else
%                     max_dist(i,cell_count) = 0;
%                     max_dist_loc(i,cell_count,:) = NaN;
%                 end
%                 cell_count = cell_count + 1;
    %             c1_targs = bsxfun(@minus,result_full_nrp(j).spike_targ_pos,result_full_nrp(j).c1_pos);
    %             scatter3(c1_targs(:,1),c1_targs(:,2),c1_targs(:,3),5,'k','filled')
                hold on
%                 z_spike_locs{i} = [z_spike_locs{i}; abs(these_locs_cell(:,3))];
                scatter3(these_locs(:,1)-result_full_nrp(j).c1_pos(1),these_locs(:,2)-result_full_nrp(j).c1_pos(2),these_locs(:,3)-result_full_nrp(j).c1_pos(3),[],...
                    [0 0 1],'filled','MarkerFaceAlpha',1./2.)
                hold on
            end
%             scatter3(result_full_nrp(j).c1_pos(1),result_full_nrp(j).c1_pos(2),result_full_nrp(j).c1_pos(3),60,[1 0 0])
            
%             c2_targs = bsxfun(@minus,result_full_nrp(j).spike_targ_pos,result_full_nrp(j).c2_pos);
%             scatter3(c2_targs(:,1),c2_targs(:,2),c2_targs(:,3),5,'k','filled')
            if ch2_cell(j)
                these_locs = result_full_nrp(j).spike_targ_pos(these_trials' & ~isnan(result_full_nrp(j).spike_times_c2) &...
                    result_full_nrp(j).spike_times_c2 < spike_time_max,:);
                these_locs_cell = result_full_nrp(j).spike_targ_pos_c2(these_trials' & ~isnan(result_full_nrp(j).spike_times_c2) &...
                    result_full_nrp(j).spike_times_c2 < spike_time_max,:);
                z_aligned = these_locs_cell(:,1) > -xy_bound & these_locs_cell(:,1) < xy_bound & these_locs_cell(:,2) > -xy_bound & these_locs_cell(:,2) < xy_bound;
                these_locs_cell = these_locs_cell(z_aligned,:);
                these_locs_dist = sqrt(sum(these_locs_cell.^2,2));
%                 if ~isempty(these_locs_dist)
%                     if k == 1
%                         [max_dist(i,cell_count),max_i] = max(these_locs_dist);
%                         max_dist_loc(i,cell_count,:) = these_locs_cell(max_i,:);
%                     end
% %                     scatter3(max_dist_loc(i,cell_count,1),max_dist_loc(i,cell_count,2),-max_dist_loc(i,cell_count,3),[],[0 0 1],'filled','MarkerFaceAlpha',1./4.)
%                 else
%                     max_dist(i,cell_count) = 0;
%                     max_dist_loc(i,cell_count,:) = NaN;
%                 end
                cell_count = cell_count + 1;
                z_spike_locs{i} = [z_spike_locs{i}; abs(these_locs_cell(:,3))];
                scatter3(these_locs(:,1)-result_full_nrp(j).c2_pos(1),these_locs(:,2)-result_full_nrp(j).c2_pos(2),these_locs(:,3)-result_full_nrp(j).c2_pos(3),[],...
                    [0 0 1],'filled','MarkerFaceAlpha',1./2.)
                hold on
            end
            xlim([-20 20]); ylim([-20 20]); zlim([-50 50])
%             zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
            xlabel('vertical')
            ylabel('horizontal')
            zlabel('axial/horizontal')
            title(sprintf('All Cells, Power: %d',unique_powers(i)))
            view(views(k,1),views(k,2))
            axis image
        end
%         subplot(3,length(unique_powers)+1,length(unique_powers)+1+(length(unique_powers)+1)*(k-1))
%         scatter3(result_full_nrp(j).spike_targ_pos(:,1),result_full_nrp(j).spike_targ_pos(:,2),result_full_nrp(j).spike_targ_pos(:,3),'filled')
%         xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
%         ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
%         zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
%         xlabel('vertical')
%         ylabel('horizontal')
%         zlabel('axial/horizontal')
%         title(sprintf('Pair %d, All Targets',j))
%         view(views(k,1),views(k,2))
    end
end

%%
test_locs_bin = zeros(3,10);
spike_locs_bin = zeros(3,10);
for i = 1:3
    test_locs_bin(i,:) = histcounts(z_test_locs{i},0:5:50);
    spike_locs_bin(i,:) = histcounts(z_spike_locs{i},0:5:50);
end
figure; plot(5:5:50,spike_locs_bin'./test_locs_bin')

%%
max_dist(:,[13:15 18]) = [];

max_dist_loc(:,[13:15 18],:) = [];


figure; cdfplot(max_dist(1,:)); hold on; cdfplot(max_dist(2,:)); cdfplot(max_dist(3,:)); legend({'25 mW','50 mW','75 mW'})
ylabel('Fraction of Cells')
xlabel('Distance (um)')

distance_x = 0:1:50;
percent_below = zeros(2,length(distance_x));

for pow_i = 1:length(unique_powers)
    for i = 1:length(distance_x)
        percent_below(pow_i,i) = sum(abs(max_dist_loc(pow_i,:,3)) > distance_x(i))/size(max_dist_loc,2);
    end
end
figure;
plot(distance_x,bsxfun(@times,percent_below,[.75/.9 .83 .98]'))
legend({'25 mW','50 mW','75 mW'})
xlabel('Axial Distance (um)')
ylabel('Fraction Cells Spiking')


%%
% z_bounds = 0:5:50;

num_followers = zeros(size(percent_below,1),length(distance_x));
num_within_range_mat = [num_within_range20; num_within_range20; num_within_range25; num_within_range25];

bounds = [20 25 30];

density_cells = 2.4889e-04;

for j = 1:size(percent_below,1)
    for i = 1:length(distance_x)
        num_followers(j,i) = percent_below(j,i)*density_cells*bounds(j)^2*2;
    end
end
figure;
plot(distance_x,num_followers')

% figure
% plot(z_bounds(2:end),num_within_range_mat')

figure;
bar([.7273 .8309 0.98;sum(num_followers(1:3,8:end),2)']')
hold on; plot([-1 4],[1.0 1.0],'k--');
xlim([0.5 3.5])
yticks([0:14])
sum(num_followers(1:3,8:end),2)

xticklabels({'25 mW','50 mW','75 mW'})
xlabel('target power')
ylabel('Spiking Cells')
% legend({'On Target','Off Target'})

%% plot min spike time loc

spike_time_max = 140;
spike_time_min = 0;
colors = jet(total_cells);

views = [0 90; 0 0; 90 0];
        for j = 1:size(filenames_nrp,1) 

    figure; 

     1
    unique_powers = unique(result_full_nrp(j).spike_targ_power);

for k = 1:3
    for i = 1:length(unique_powers)

        count = 1;
            subplot(3,length(unique_powers),i + length(unique_powers)*(k-1))
            
            if ch1_cell(j)
                if ~any(result_full_nrp(j).spike_times_c1 < spike_time_min)
                    these_trials = find(result_full_nrp(j).spike_targ_power' == unique_powers(i) & ~isnan(result_full_nrp(j).spike_times_c1) &...
                        result_full_nrp(j).spike_times_c1 < spike_time_max);
    %                 these_locs = result_full_nrp(j).spike_targ_pos(these_trials,:);
                    these_times = result_full_nrp(j).spike_times_c1(these_trials);
                    if ~isempty(these_times)
        %             c1_targs = bsxfun(@minus,result_full_nrp(j).spike_targ_pos,result_full_nrp(j).c1_pos);
        %             scatter3(c1_targs(:,1),c1_targs(:,2),c1_targs(:,3),5,'k','filled')
                        [this_time, this_ind] = min(these_times);
                        this_loc = result_full_nrp(j).spike_targ_pos(these_trials(this_ind),:) - result_full_nrp(j).c1_pos + rand(1,3)-.5;
    %                     scatter3(these_locs(:,1),these_locs(:,2),these_locs(:,3),5,[0 0 0],'filled','jitter','on','jitteramount',0.)
    %                     hold on
                        scatter3(this_loc(:,1),this_loc(:,2),this_loc(:,3),[],colors(count,:),'filled','jitter','on','jitteramount',0.)
                        hold on
        %             scatter3(result_full_nrp(j).c1_pos(1),result_full_nrp(j).c1_pos(2),result_full_nrp(j).c1_pos(3),60,[1 0 0])

        %             c2_targs = bsxfun(@minus,result_full_nrp(j).spike_targ_pos,result_full_nrp(j).c2_pos);
        %             scatter3(c2_targs(:,1),c2_targs(:,2),c2_targs(:,3),5,'k','filled')
                    end
                end
                count = count + 1;
            end
            if ch2_cell(j)
                if ~any(result_full_nrp(j).spike_times_c2 < spike_time_min)

                    these_trials = find(result_full_nrp(j).spike_targ_power' == unique_powers(i) & ~isnan(result_full_nrp(j).spike_times_c2) &...
                        result_full_nrp(j).spike_times_c2 < spike_time_max);
                    these_locs = result_full_nrp(j).spike_targ_pos(these_trials,:);
                    these_times = result_full_nrp(j).spike_times_c2(these_trials);
                    if ~isempty(these_times)
        %             c1_targs = bsxfun(@minus,result_full_nrp(j).spike_targ_pos,result_full_nrp(j).c1_pos);
        %             scatter3(c1_targs(:,1),c1_targs(:,2),c1_targs(:,3),5,'k','filled')
                        [this_time, this_ind] = min(these_times);
                        this_loc = result_full_nrp(j).spike_targ_pos(these_trials(this_ind),:) - result_full_nrp(j).c2_pos + rand(1,3)-.5;
    %                     scatter3(these_locs(:,1),these_locs(:,2),these_locs(:,3),5,[0 0 0],'filled','jitter','on','jitteramount',0.)
    %                     hold on
                        scatter3(this_loc(:,1),this_loc(:,2),this_loc(:,3),[],colors(count,:),'filled')
                        hold on
        %             scatter3(result_full_nrp(j).c1_pos(1),result_full_nrp(j).c1_pos(2),result_full_nrp(j).c1_pos(3),60,[1 0 0])

        %             c2_targs = bsxfun(@minus,result_full_nrp(j).spike_targ_pos,result_full_nrp(j).c2_pos);
        %             scatter3(c2_targs(:,1),c2_targs(:,2),c2_targs(:,3),5,'k','filled')

                    end
                end
                count = count + 1;
            end
            xlabel('vertical')
            ylabel('horizontal')
            zlabel('axial/horizontal')
            title(sprintf('All Cells, Power: %d',unique_powers(i)))
            view(views(k,1),views(k,2))
            axis image
            xlim([-20 20])
            ylim([-20 20])
            zlim([-45 45])
        end
%         subplot(3,length(unique_powers)+1,length(unique_powers)+1+(length(unique_powers)+1)*(k-1))
%         scatter3(result_full_nrp(j).spike_targ_pos(:,1),result_full_nrp(j).spike_targ_pos(:,2),result_full_nrp(j).spike_targ_pos(:,3),'filled')
%         xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
%         ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
%         zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
%         xlabel('vertical')
%         ylabel('horizontal')
%         zlabel('axial/horizontal')
%         title(sprintf('Pair %d, All Targets',j))
%         view(views(k,1),views(k,2))
    end
end

%% z- sections


offset = .00;

spike_time_max = 140;

views = [90 90; 90 -0];
subplot_order = [1 2 3];
colors = [.35 0 0; .68 0 0; 1 0 0];
colors = [.2 .2 1; 1 .2 .2; .2 1 .2; .2 .2 .2]*.75;

xlimits = [-30 30];
ylimits = [-30 30];
zlimits = [-50 50];

time_thresh = 140;

figure
n_choice = [5];
count = 1;
for jj = 1:length(n_choice)
    j = n_choice(jj);
%     figure
    j
    unique_powers = unique(result_full_nrp(j).spike_targ_power);%[15 35 55]; %
    unique_powers(unique_powers > 55) = [];
%     [unique_targs,~,unique_targs_trial_idx] = unique(result_full_nrp(j).spike_targ_pos_c1,'rows');
    
%     prob_spike_c2 = zeros(length(unique_powers),size(unique_targs,1));
    for i = length(unique_powers)
        i
        for ii = 1:5
            z_range = -50 + [(ii-1)*20 ii*20]
        for k = 1:size(views,1)
            k
            subplot(2,5,ii+5*(k-1))
            disp(['subplot : ' num2str(ii+5*(k-1))])
            hold on
            if k == 1 && jj == 1
                line(xlimits,[0 0],[0 0],'color','k','linewidth',0.5)
                hold on
                line([0 0],ylimits,[0 0],'color','k','linewidth',0.5)
                hold on
            elseif k == 2 && jj == 1
                line([0 0],[0 0],zlimits,'color','k','linewidth',0.5)
                hold on
                line([0 0],ylimits,[0 0],'color','k','linewidth',0.5)
                hold on
            end
%             subplot(3,6,count+6*(j-1))
%             subplot(2,2,subplot_order(k))
%             these_trials = result_full_nrp(j).targ_power == unique_powers(i);
%             these_trials = these_trials & result_full_nrp(j).spike_targ_pos_c1(:,3)' >= z_range(1) & result_full_nrp(j).spike_targ_pos_c1(:,3)' < z_range(2);
%             [unique_targs,~,unique_targs_trial_idx] = unique(result_full_nrp(j).spike_targ_pos_c1(these_trials,:),'rows');
            
%             hold on
            
            
            num_cells = 0;
            
            if ch2_cell(j) == 1
                
                these_trials = result_full_nrp(j).spike_targ_power == unique_powers(i);
            these_trials = these_trials & result_full_nrp(j).spike_targ_pos_c2(:,3)' >= z_range(1) & result_full_nrp(j).spike_targ_pos_c2(:,3)' < z_range(2);
            [unique_targs,~,unique_targs_trial_idx] = unique(result_full_nrp(j).spike_targ_pos_c2(these_trials,:),'rows');
%                 prob_spike_c2 = zeros(length(unique_powers),size(unique_targs,1));
%                 scatter3(result_full_nrp(j).c2_pos(1),result_full_nrp(j).c2_pos(2),result_full_nrp(j).c2_pos(3),100,[0 0 1],'filled')
                these_times = result_full_nrp(j).spike_times_c2(these_trials);
                for this_loc_ind = 1:size(unique_targs,1)
                    loc_times = these_times(unique_targs_trial_idx == this_loc_ind);
                    this_prob_spike = sum(~isnan(loc_times) & loc_times < time_thresh)/length(loc_times);
                    if this_prob_spike
                        scatter3(unique_targs(this_loc_ind,1),unique_targs(this_loc_ind,2),-unique_targs(this_loc_ind,3),...
                            40,colors(count,:),'filled','markerfacealpha',this_prob_spike^1.2-.1)
                        hold on
                    else
                        scatter3(unique_targs(this_loc_ind,1),unique_targs(this_loc_ind,2),-unique_targs(this_loc_ind,3),...
                            10,colors(count,:),'filled','MarkerFaceAlpha',.2)
                        hold on
                    end
%                     prob_spike_c2(i,this_loc_ind) = this_prob_spike;
                end
                num_cells = num_cells + 1;
            end
            
            if ch1_cell(j) == 1
                these_trials = result_full_nrp(j).spike_targ_power == unique_powers(i);
            these_trials = these_trials & result_full_nrp(j).spike_targ_pos_c1(:,3)' >= z_range(1) & result_full_nrp(j).spike_targ_pos_c1(:,3)' < z_range(2);
            [unique_targs,~,unique_targs_trial_idx] = unique(result_full_nrp(j).spike_targ_pos_c1(these_trials,:),'rows');
%             prob_spike_c1 = zeros(1,size(unique_targs,1));    
%             scatter3(result_full_nrp(j).c1_pos(1),result_full_nrp(j).c1_pos(2),result_full_nrp(j).c1_pos(3),40,[1 0 0])
                these_times = result_full_nrp(j).spike_times_c1(these_trials);
                
                for this_loc_ind = 1:size(unique_targs,1)
                    loc_times = these_times(unique_targs_trial_idx == this_loc_ind);
                    this_prob_spike = sum(~isnan(loc_times) & loc_times < time_thresh)/length(loc_times);
                    if this_prob_spike
                        scatter3(unique_targs(this_loc_ind,1),unique_targs(this_loc_ind,2),-unique_targs(this_loc_ind,3),40,colors(count+1,:),'filled','MarkerFaceAlpha',this_prob_spike^1.2-.1)
                        hold on
                    else
                        scatter3(unique_targs(this_loc_ind,1),unique_targs(this_loc_ind,2),-unique_targs(this_loc_ind,3),...
                            10,colors(count+1,:),'filled','MarkerFaceAlpha',.2)
                        hold on
                    end
%                     prob_spike_c1(i,this_loc_ind) = this_prob_spike;
                end
                num_cells = num_cells + 1;
            end
%             if i == length(unique_powers)
%                 spiking_targs = unique_targs(prob_spike_c2(i,:) > 0,:);
%             end
              if jj == length(n_choice)
                if k == 1
                    title(['Power: ' num2str(unique_powers(i)) 'mW'])
                end
                if ii == 1 && k == 1
                    xlabel('vertical')
                    set(gca,'Yticklabel',[]) 
                end
                if k == 2
                    ylabel('horizontal')
                end

                if ii == 1 && k == 2
                    zlabel('axial/horizontal')
                end
                if ii > 1 && ii < 5 && k == 1
                    set(gca,'Xticklabel',[]) 
                    set(gca,'Yticklabel',[]) 
                end
                if ii > 1 && ii < 5 && k == 2
                    set(gca,'Zticklabel',[]) 
                end
                if ii == 5 && k == 1
                    set(gca,'Yticklabel',[]) 
                    set(gca,'Xticklabel',[]) 
                end
                if ii == 5 && k == 2
                    set(gca,'Zticklabel',[]) 
                end
%                 title(sprintf('Pair %d, Power: %d',j,unique_powers(i)))
                view(views(k,1),views(k,2))
                axis equal
                xlim(xlimits)
                ylim(ylimits)
                zlim(zlimits)
%                 lim = axis;
                if k == 1
                    axis_pos_radial = get(gca, 'Position')
%                     lim_axial = axis;
                elseif k == 2
                    axis_pos_axial = get(gca, 'Position')   
%                     lim_horiz = axis;
                    set(gca, 'Position', [axis_pos_axial(1:3) diff(zlimits)/diff(xlimits)*axis_pos_radial(4)]) 
%                 else
%                     axis_pos3 = get(gca, 'Position') ;                    
%                     set(gca, 'Position', [axis_pos3(1:2)  (lim(6)-lim(5))/(lim(2)-lim(1))*axis_pos2(3)]) 
                end
%             end
%             axis image
              end
                end
%             count = count - 1;
        end
%         subplot(3,length(unique_powers)+1,length(unique_powers)+1+(length(unique_powers)+1)*(k-1))
%         scatter3(result_full_nrp(j).spike_targ_pos(:,1),result_full_nrp(j).spike_targ_pos(:,2),result_full_nrp(j).spike_targ_pos(:,3),'filled')
%         xlim([min(result_full_nrp(j).spike_targ_pos(:,1)) max(result_full_nrp(j).spike_targ_pos(:,1))])
%         ylim([min(result_full_nrp(j).spike_targ_pos(:,2)) max(result_full_nrp(j).spike_targ_pos(:,2))])
%         zlim([min(result_full_nrp(j).spike_targ_pos(:,3)) max(result_full_nrp(j).spike_targ_pos(:,3))])
%         xlabel('vertical')
%         ylabel('horizontal')
%         zlabel('axial/horizontal')
%         title(sprintf('Pair %d, All Targets',j))
%         view(views(k,1),views(k,2))
        
    end
    count = count + num_cells;
end

