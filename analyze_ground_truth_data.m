%%

filenames = {'10_25_slice1_cell1_2.mat','10_25_15_33_data.mat'
             '10_25_slice1_cell2_4.mat','10_25_15_54_data.mat'
             '10_25_slice2_cell1_2.mat','10_25_17_21_data.mat'
             };      
                                                    
ch1_cell_type = [2,2,1]; %0-no cell, 1-spike cell, 2-psc cell
ch2_cell_type = [1,1,1];

thresh = [30 30 30;
          25 30 7];
      
hpf = [0 0 0;
       0 0 1];
               
map_trials = {3:5,3:5,2:5};
connection_check_trials = {1,1,1};
post_aspiration_trials = {6,6,6};


%%
% thisdir = '~/projects/mapping/data/';
thisdir = '/media/shababo/data/';
%%

for j = 1:size(filenames,1)
    
    j
    
    % load data
    load([thisdir filenames{j,2}]);
    load([thisdir filenames{j,1}]); 
    experiment_setup = exp_data.experiment_setup;
    experiment_setup.trials.min_time = 30;
    experiment_setup.trials.max_time = 200;
    data_trials = map_trials{j};
    [result_ground_truth_set(j).traces_c1, result_ground_truth_set(j).traces_c2, this_seq, result_ground_truth_set(j).stim_traces, result_ground_truth_set(j).full_stim_key] = get_traces(data,data_trials);
    result_ground_truth_set(j).targ_pos = result_ground_truth_set(j).full_stim_key([this_seq.precomputed_target_index],:);
    result_ground_truth_set(j).targ_power = [this_seq.target_power];
    
    filename_base = ['/media/shababo/data/' experiment_setup.exp_id];
    % detect spikes
    disp('cell 1')
    if ch1_cell_type(j) == 1
        

        result_ground_truth_set(j).c1_targ_pos = bsxfun(@minus,result_ground_truth_set(j).full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc);
        result_ground_truth_set(j).c1_pos = experiment_setup.patched_cell_loc;
        
        cell_spike_times_c1 = detect_peaks(-bsxfun(@minus,result_ground_truth_set(j).traces_c1(:,30:end),median(result_ground_truth_set(j).traces_c1(:,30:end),2)),thresh(1,j),30,0,Inf,-Inf,0,hpf(1,j),1);
        result_ground_truth_set(j).spike_times_c1 = zeros(size(cell_spike_times_c1));
        for i = 1:length(cell_spike_times_c1)
            if ~isempty(cell_spike_times_c1{i})
                result_ground_truth_set(j).spike_times_c1(i) = cell_spike_times_c1{i} + 29;
            else
                result_ground_truth_set(j).spike_times_c1(i) = NaN;
            end
        end
    elseif ch1_cell_type(j) == 2
        
        %result_ground_truth_set(j).c1_targ_pos = bsxfun(@minus,result_ground_truth_set(j).full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc);
        result_ground_truth_set(j).c1_pos = experiment_setup.patched_cell_loc;
        
        fullsavepath = [filename_base '_traces.mat'];
        oasis_out_path = [filename_base '_traces_detect.mat'];
        traces = result_ground_truth_set(j).traces_c1;
        save(fullsavepath,'traces')
        cmd = 'python /home/shababo/projects/mapping/code/OASIS/run_oasis_online.py ';
        cmd = [cmd ' ' fullsavepath];
        error = system(cmd);
        if ~error
            % Wait for file to be created.
            maxSecondsToWait = 60*5; % Wait five minutes...?
            secondsWaitedSoFar  = 0;
            while secondsWaitedSoFar < maxSecondsToWait 
              if exist(oasis_out_path, 'file')
                break;
              end
              pause(1); % Wait 1 second.
              secondsWaitedSoFar = secondsWaitedSoFar + 1;
            end
            if exist(oasis_out_path, 'file')
              load(oasis_out_path)
              oasis_data = reshape(event_process,size(traces'))';

            else
              fprintf('Warning: x.log never got created after waiting %d seconds', secondsWaitedSoFar);
            %               uiwait(warndlg(warningMessage));
                oasis_data = zeros(size(traces));
            end
        else
          oasis_data = zeros(size(traces));
        end
          
        for jj = 1:size(traces,1)
            if ~isempty(find(oasis_data(jj,...
                        experiment_setup.trials.min_time:experiment_setup.trials.max_time),1))
                result_ground_truth_set(j).event_times_c1(jj) = ...
                    find(oasis_data(jj,...
                        experiment_setup.trials.min_time:experiment_setup.trials.max_time),1) + ...
                        experiment_setup.trials.min_time - 1;
            else
                result_ground_truth_set(j).event_times_c1(jj) = NaN;
            end
        end
    end 
        disp('cell 2')
    if ch2_cell_type(j) == 1
        

        result_ground_truth_set(j).c2_targ_pos = bsxfun(@minus,result_ground_truth_set(j).full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc_2);
        result_ground_truth_set(j).c2_pos = experiment_setup.patched_cell_loc_2;
        
        cell_spike_times_c2 = detect_peaks(-bsxfun(@minus,result_ground_truth_set(j).traces_c2(:,30:end),median(result_ground_truth_set(j).traces_c2(:,30:end),2)),thresh(2,j),30,0,Inf,-Inf,0,hpf(2,j),1);
        result_ground_truth_set(j).spike_times_c2 = zeros(size(cell_spike_times_c2));
        for i = 1:length(cell_spike_times_c2)
            if ~isempty(cell_spike_times_c2{i})
                result_ground_truth_set(j).spike_times_c2(i) = cell_spike_times_c2{i} + 29;
            else
                result_ground_truth_set(j).spike_times_c2(i) = NaN;
            end
        end
    elseif ch2_cell_type(j) == 2
        
        result_ground_truth_set(j).c2_targ_pos = bsxfun(@minus,result_ground_truth_set(j).full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc_2);
        result_ground_truth_set(j).c2_pos = experiment_setup.patched_cell_loc_2;
        
        fullsavepath = [filename_base '_traces.mat'];
        oasis_out_path = [filename_base '_traces_detect.mat'];
        traces = result_ground_truth_set(j).traces_c2;
        save(fullsavepath,'traces')
        cmd = 'python /home/shababo/projects/mapping/code/OASIS/run_oasis_online.py ';
        cmd = [cmd ' ' fullsavepath];
        error = system(cmd);
        if ~error
            % Wait for file to be created.
            maxSecondsToWait = 60*5; % Wait five minutes...?
            secondsWaitedSoFar  = 0;
            while secondsWaitedSoFar < maxSecondsToWait 
              if exist(oasis_out_path, 'file')
                break;
              end
              pause(1); % Wait 1 second.
              secondsWaitedSoFar = secondsWaitedSoFar + 1;
            end
            if exist(oasis_out_path, 'file')
              load(oasis_out_path)
              oasis_data = reshape(event_process,size(traces'))';

            else
              fprintf('Warning: x.log never got created after waiting %d seconds', secondsWaitedSoFar);
            %               uiwait(warndlg(warningMessage));
                oasis_data = zeros(size(traces));
            end
        else
          oasis_data = zeros(size(traces));
        end
          
        for jj = 1:size(traces,1)
            result_ground_truth_set(j).event_times_c2(jj) = ...
                find(oasis_data(jj,...
                        experiment_setup.trials.min_time:experiment_setup.trials.max_time),1) + ...
                        experiment_setup.trials.min_time - 1;
        end
        
    end
    
    % detect pscs
    
    
    
end   
    
    
 %%   
    
    
load([thisdir filenames{1,2}])
% load([thisdir filenames{9,2}])
experiment_setup = exp_data.experiment_setup;

local_nucs = experiment_setup.local_nuc_locs;
all_nucs = experiment_setup.nuclear_locs;

presyn_cell_pos = experiment_setup.presyn_cell_pos;
% if experiment_setup.ch1_pre
%     postsyn_cell_pos = experiment_setup.patched_cell_loc_2;
% else
%     postsyn_cell_pos = experiment_setup.patched_cell_loc_1;
% end

figure; scatter3(all_nucs(:,1), all_nucs(:,2), -all_nucs(:,3));
hold on
scatter3(local_nucs(:,1), local_nucs(:,2), -local_nucs(:,3)); 
scatter3(presyn_cell_pos(:,1), presyn_cell_pos(:,2), -presyn_cell_pos(:,3)); axis image
% scatter3(presyn_cell_pos(:,1), presyn_cell_pos(:,2), -presyn_cell_pos(:,3)); axis image

%%

% load(['~/projects/mapping/data/' filenames_nrp{9,2}])
% experiment_setup = exp_data.experiment_setup;

% local_nucs = experiment_setup.local_nuc_locs;
all_nucs = experiment_setup.nuclear_locs;
num_cells = size(all_nucs,1);
xy_bounds = [15];
z_bounds = 0:5:50;

total_cells = sum(all_nucs(:,1) > -100 & all_nucs(:,1) < 100 & all_nucs(:,2) > -100 & all_nucs(:,2) < 150 & all_nucs(:,3) > 20 & all_nucs(:,3) < 100);
density_cell = total_cells/(200 * 250 * 80)


num_within_range = zeros(num_cells,length(z_bounds)-1);
for i = 1:length(z_bounds)-1
    for j = 1:num_cells
        all_nucs_rel = abs(bsxfun(@minus,all_nucs,all_nucs(j,:)));
        num_within_range(j,i) = sum(all_nucs_rel(:,1) < xy_bounds & all_nucs_rel(:,2) < xy_bounds & all_nucs_rel(:,3) < z_bounds(i+1) & all_nucs_rel(:,3) >= z_bounds(i));
    end
    
end
num_within_range(:,1) = num_within_range(:,1)-1;
% figureall_nucs_rel
% boxplot(num_within_range)


num_within_range_all = mean(num_within_range,1);

figure;
plot(z_bounds(2:end),num_within_range_all)

num_within_range25 = num_within_range_all;

%%
paired_trials = ~isnan(result_ground_truth_set(2).spike_times_c2') & ~isnan(result_ground_truth_set(2).event_times_c1) & result_ground_truth_set(2).spike_times_c2' < result_ground_truth_set(2).event_times_c1;


figure;
subplot(121)
scatter(result_ground_truth_set(2).spike_times_c2(paired_trials),...
    result_ground_truth_set(2).event_times_c1(paired_trials));
hold on
% plot(0:1:200,0:1:200)
% subplot(132); histogram(result_ground_truth_set(1).event_times_c1(paired_trials))
% hold on
% histogram(result_ground_truth_set(1).spike_times_c2(paired_trials))

subplot(122);
histogram(result_ground_truth_set(2).event_times_c1(paired_trials)'...
    - result_ground_truth_set(2).spike_times_c2(paired_trials))
%%

figure;
unique_powers = unique(result_ground_truth_set(1).targ_power);
jitter_amt = .5;
for i = 1:length(unique_powers)

    subplot(1,length(unique_powers),i)
    trials = result_ground_truth_set(1).targ_power == unique_powers(i) & ~isnan(result_ground_truth_set(1).event_times_c1) & ...
        result_ground_truth_set(1).event_times_c1 > 50  & result_ground_truth_set(1).event_times_c1 < 160;
    scatter3(result_ground_truth_set(1).c2_targ_pos(trials,1)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        result_ground_truth_set(1).c2_targ_pos(trials,2)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        -result_ground_truth_set(1).c2_targ_pos(trials,3)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
    80,200 - result_ground_truth_set(1).event_times_c1(trials),'filled','markerfacealpha',.75)

end
hold on
scatter3(experiment_setup.local_nuc_locs(:,1),experiment_setup.local_nuc_locs(:,2),-experiment_setup.local_nuc_locs(:,3))
%%

figure;

unique_powers = unique(result_ground_truth_set(2).targ_power);
jitter_amt = .5;
for i = 1:length(unique_powers)
    subplot(1,length(unique_powers),i)
%     subplot(1,length(unique_powers),i)
    trials = result_ground_truth_set(2).targ_power' == unique_powers(i);
%     scatter3(experiment_setup.local_nuc_locs(:,1),experiment_setup.local_nuc_locs(:,2)-6,-(experiment_setup.local_nuc_locs(:,3)-6),100,[1 0 0],'filled','markerfacealpha',0.15)
    hold on
    scatter3(result_ground_truth_set(2).c2_targ_pos(trials,1),...
        result_ground_truth_set(2).c2_targ_pos(trials,2),...
        -result_ground_truth_set(2).c2_targ_pos(trials,3),...
        5,[0 0 1],'filled','markerfacealpha',.25)
    hold on
    scatter3(0,...
        0,...
        0,20,...
        'x')
    trials = result_ground_truth_set(2).targ_power' == unique_powers(i) & ~isnan(result_ground_truth_set(2).spike_times_c2);
    scatter3(result_ground_truth_set(2).c2_targ_pos(trials,1)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        result_ground_truth_set(2).c2_targ_pos(trials,2)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        -result_ground_truth_set(2).c2_targ_pos(trials,3)+rand(sum(trials),1)*jitter_amt+jitter_amt/2,...
        50,200 - result_ground_truth_set(2).spike_times_c2(trials),'filled','markerfacealpha',.75)
    xlabel('vertical')
    ylabel('horiz')
    zlabel('axial')
    caxis([0 150])
    axis image
end

%%

% figure;

unique_powers = unique(result_ground_truth_set(1).targ_power);
jitter_amt = .5;
for i = 1:length(unique_powers)
    subplot(1,length(unique_powers),i)%+length(unique_powers))
%     subplot(1,length(unique_powers),i)
    trials = result_ground_truth_set(1).targ_power' == unique_powers(i) & ~isnan(result_ground_truth_set(1).spike_times_c2);
%     scatter3(experiment_setup.local_nuc_locs(:,1),experiment_setup.local_nuc_locs(:,2),-experiment_setup.local_nuc_locs(:,3),100,[1 0 0],'filled','markerfacealpha',0.15)
    hold on
    scatter3(result_ground_truth_set(1).c2_targ_pos(:,1),...
        result_ground_truth_set(1).c2_targ_pos(:,2),...
        -result_ground_truth_set(1).c2_targ_pos(:,3),...
        5,[0 0 1],'filled','markerfacealpha',.25)
    hold on
    scatter3(result_ground_truth_set(1).c2_targ_pos(trials,1)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        result_ground_truth_set(1).c2_targ_pos(trials,2)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        -result_ground_truth_set(1).c2_targ_pos(trials,3)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        50,250 - result_ground_truth_set(1).event_times_c1(trials),'filled','markerfacealpha',.75)
    caxis([0 200])
    axis image
end


%%
figure
unique_powers = unique(result_ground_truth_set(1).targ_power);
jitter_amt = .5;
for i = 1:length(unique_powers)
    subplot(1,length(unique_powers),i)%+2*length(unique_powers)
%     subplot(1,length(unique_powers),i)
    trials = result_ground_truth_set(1).targ_power == unique_powers(i) & result_ground_truth_set(1).event_times_c1 < 200  & result_ground_truth_set(1).event_times_c1 > 50;
    scatter3(experiment_setup.local_nuc_locs(:,1),experiment_setup.local_nuc_locs(:,2),-experiment_setup.local_nuc_locs(:,3),125,[1 0 0])
    hold on
%     scatter3(result_ground_truth_set(1).c2_targ_pos(:,1),...
%         result_ground_truth_set(1).c2_targ_pos(:,2),...
%         -result_ground_truth_set(1).c2_targ_pos(:,3),...
%         5,[0 0 1],'filled','markerfacealpha',.25)
%     hold on
    scatter3(result_ground_truth_set(1).targ_pos(trials,1)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        result_ground_truth_set(1).targ_pos(trials,2)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        -result_ground_truth_set(1).targ_pos(trials,3)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        50,200 - result_ground_truth_set(1).event_times_c1(trials),'filled','markerfacealpha',.75)
    caxis([0 150])
    axis image
end


%% plot spiking locs

offset = .00;

spike_time_max = 140;

views = [90 90; 90 -0];
subplot_order = [1 2 3];
colors = [.35 0 0; .68 0 0; 1 0 0];
colors = [.2 .2 1; .2 .2 1; .2 .2 1]*.75;

xlimits = [-40 40];
ylimits = [-30 30];
zlimits = [-50 50];

figure
n_choice = [2];
for jj = 1:length(n_choice)
    j = n_choice(jj);
%     figure
    j
    unique_powers = unique(result_ground_truth_set(j).targ_power);%[15 35 55]; %
    unique_powers(unique_powers > 55) = [];
    [unique_targs,~,unique_targs_trial_idx] = unique(result_ground_truth_set(j).c2_targ_pos,'rows');
    prob_spike_c1 = zeros(length(unique_powers),size(unique_targs,1));
    prob_spike_c2 = zeros(length(unique_powers),size(unique_targs,1));
    for i = 1:length(unique_powers) 
        i
        for k = 1:size(views,1)
            k
            subplot(2,3,i+3*(k-1))
            disp(['subplot : ' num2str(i+3*(k-1))])
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
            these_trials = result_ground_truth_set(j).targ_power == unique_powers(i);
            [unique_targs,~,unique_targs_trial_idx] = unique(result_ground_truth_set(j).c2_targ_pos(these_trials,:),'rows');
            
%             hold on
            
            
            
            
            if ch2_cell_type(j) == 1
                
%                 scatter3(result_ground_truth_set(j).c2_pos(1),result_ground_truth_set(j).c2_pos(2),result_ground_truth_set(j).c2_pos(3),100,[0 0 1],'filled')
                these_times = result_ground_truth_set(j).spike_times_c2(these_trials);
                for this_loc_ind = 1:size(unique_targs,1)
                    loc_times = these_times(unique_targs_trial_idx == this_loc_ind);
                    this_prob_spike = sum(~isnan(loc_times) & loc_times < spike_time_max)/length(loc_times);
                    if this_prob_spike
                        scatter3(unique_targs(this_loc_ind,1),unique_targs(this_loc_ind,2),-unique_targs(this_loc_ind,3),...
                            40,colors(j,:),'filled','markerfacealpha',this_prob_spike^1.2)
                        hold on
                    else
                        scatter3(unique_targs(this_loc_ind,1),unique_targs(this_loc_ind,2),-unique_targs(this_loc_ind,3),...
                            10,colors(j,:),'filled','MarkerFaceAlpha',.2)
                    end
                    prob_spike_c2(i,this_loc_ind) = this_prob_spike;
                end
            end
            
%             if ch1_cell_type(j) == 1
%                 scatter3(result_ground_truth_set(j).c1_pos(1),result_ground_truth_set(j).c1_pos(2),result_ground_truth_set(j).c1_pos(3),40,[1 0 0])
%                 these_times = result_ground_truth_set(j).spike_times_c1(these_trials);
%                 
%                 for this_loc_ind = 1:size(unique_targs,1)
%                     loc_times = these_times(unique_targs_trial_idx == this_loc_ind);
%                     this_prob_spike = sum(~isnan(loc_times) & loc_times < 140)/length(loc_times);
%                     if this_prob_spike
%                         scatter3(unique_targs(this_loc_ind,1),unique_targs(this_loc_ind,2),unique_targs(this_loc_ind,3),30,[1 0 0],'filled','MarkerFaceAlpha',1)
%                         hold on
%                     end
%                     prob_spike_c1(i,this_loc_ind) = this_prob_spike;
%                 end
%             end
%             if i == length(unique_powers)
%                 spiking_targs = unique_targs(prob_spike_c2(i,:) > 0,:);
%             end
              if jj == length(n_choice)
                if k == 1
                    title(['Power: ' num2str(unique_powers(i)) 'mW'])
                end
                if i == 1 && k == 1
                    xlabel('vertical')
                    set(gca,'Yticklabel',[]) 
                end
                if k == 2
                    ylabel('horizontal')
                end

                if i == 1 && k == 2
                    zlabel('axial/horizontal')
                end
                if i == 2 && k == 1
                    set(gca,'Xticklabel',[]) 
                    set(gca,'Yticklabel',[]) 
                end
                if i == 2 && k == 2
                    set(gca,'Zticklabel',[]) 
                end
                if i == 3 && k == 1
                    set(gca,'Yticklabel',[]) 
                    set(gca,'Xticklabel',[]) 
                end
                if i == 3 && k == 2
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
            count = count - 1;
        end
%         subplot(3,length(unique_powers)+1,length(unique_powers)+1+(length(unique_powers)+1)*(k-1))
%         scatter3(result_ground_truth_set(j).spike_targ_pos(:,1),result_ground_truth_set(j).spike_targ_pos(:,2),result_ground_truth_set(j).spike_targ_pos(:,3),'filled')
%         xlim([min(result_ground_truth_set(j).spike_targ_pos(:,1)) max(result_ground_truth_set(j).spike_targ_pos(:,1))])
%         ylim([min(result_ground_truth_set(j).spike_targ_pos(:,2)) max(result_ground_truth_set(j).spike_targ_pos(:,2))])
%         zlim([min(result_ground_truth_set(j).spike_targ_pos(:,3)) max(result_ground_truth_set(j).spike_targ_pos(:,3))])
%         xlabel('vertical')
%         ylabel('horizontal')
%         zlabel('axial/horizontal')
%         title(sprintf('Pair %d, All Targets',j))
%         view(views(k,1),views(k,2))
        
    end
end

