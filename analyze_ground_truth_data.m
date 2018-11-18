%%

filenames = {'10_25_slice1_cell1_2.mat','10_25_15_33_data.mat'
             '10_25_slice1_cell2_4.mat','10_25_15_54_data.mat'
             '10_25_slice1_cell2_5.mat','10_25_15_54_data.mat'
             };      
                                                    
ch1_cell_type = [2,2,2]; %0-no cell, 1-spike cell, 2-psc cell
ch2_cell_type = [1,1,1];
               
map_trials = {3:5,3:5,3:5};
connection_check_trials = {1,1,1};
post_aspiration_trials = {6,6,6};


%%

for j = 2:size(filenames,1)
    
    j
    
    % load data
    load(filenames{j,2});
    load(filenames{j,1}); 
    experiment_setup = exp_data.experiment_setup;
    experiment_setup.trials.min_time = 30;
    experiment_setup.trials.max_time = 200;
    data_trials = map_trials{j};
    [result_full_nrp(j).traces_c1, result_full_nrp(j).traces_c2, this_seq, result_full_nrp(j).stim_traces, result_full_nrp(j).full_stim_key] = get_traces(data,data_trials);
    result_full_nrp(j).targ_pos = result_full_nrp(j).full_stim_key([this_seq.precomputed_target_index],:);
    result_full_nrp(j).targ_power = [this_seq.target_power];
    
    filename_base = ['/media/shababo/data/' experiment_setup.exp_id];
    % detect spikes
    disp('cell 1')
    if ch1_cell_type == 1
        

        result_full_nrp(j).c1_targ_pos = bsxfun(@minus,result_full_nrp(j).full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc);
        result_full_nrp(j).c1_pos = experiment_setup.patched_cell_loc;
        
        cell_spike_times_c1 = detect_peaks(-bsxfun(@minus,result_full_nrp(j).traces_c1,result_full_nrp(j).traces_c1(:,1)),20,30,0,Inf,-Inf,0,0,1);
        result_full_nrp(j).spike_times_c1 = zeros(size(cell_spike_times_c1));
        for i = 1:length(cell_spike_times_c1)
            if ~isempty(cell_spike_times_c1{i})
                result_full_nrp(j).spike_times_c1(i) = cell_spike_times_c1{i};
            else
                result_full_nrp(j).spike_times_c1(i) = NaN;
            end
        end
    elseif ch1_cell_type == 2
        
        %result_full_nrp(j).c1_targ_pos = bsxfun(@minus,result_full_nrp(j).full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc);
        result_full_nrp(j).c1_pos = experiment_setup.patched_cell_loc;
        
        fullsavepath = [filename_base '_traces.mat'];
        oasis_out_path = [filename_base '_traces_detect.mat'];
        traces = result_full_nrp(j).traces_c1;
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
                result_full_nrp(j).event_times_c1(jj) = ...
                    find(oasis_data(jj,...
                        experiment_setup.trials.min_time:experiment_setup.trials.max_time),1) + ...
                        experiment_setup.trials.min_time - 1;
            else
                result_full_nrp(j).event_times_c1(jj) = NaN;
            end
        end
    end 
        disp('cell 21')
    if ch2_cell_type == 1
        

        result_full_nrp(j).c2_targ_pos = bsxfun(@minus,result_full_nrp(j).full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc_2);
        result_full_nrp(j).c2_pos = experiment_setup.patched_cell_loc_2;
        
        cell_spike_times_c2 = detect_peaks(-bsxfun(@minus,result_full_nrp(j).traces_c2,result_full_nrp(j).traces_c2(:,1)),20,30,0,Inf,-Inf,0,0,1);
        result_full_nrp(j).spike_times_c2 = zeros(size(cell_spike_times_c2));
        for i = 1:length(cell_spike_times_c2)
            if ~isempty(cell_spike_times_c2{i})
                result_full_nrp(j).spike_times_c2(i) = cell_spike_times_c2{i};
            else
                result_full_nrp(j).spike_times_c2(i) = NaN;
            end
        end
    elseif ch2_cell_type == 2
        
        result_full_nrp(j).c2_targ_pos = bsxfun(@minus,result_full_nrp(j).full_stim_key([this_seq.precomputed_target_index],:),experiment_setup.patched_cell_loc_2);
        result_full_nrp(j).c2_pos = experiment_setup.patched_cell_loc_2;
        
        fullsavepath = [filename_base '_traces.mat'];
        oasis_out_path = [filename_base '_traces_detect.mat'];
        traces = result_full_nrp(j).traces_c2;
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
            result_full_nrp(j).event_times_c2(jj) = ...
                find(oasis_data(jj,...
                        experiment_setup.trials.min_time:experiment_setup.trials.max_time),1) + ...
                        experiment_setup.trials.min_time - 1;
        end
        
    end
    
    % detect pscs
    
    
    
end   
    
    
 %%   
    
    

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
paired_trials = ~isnan(result_full_nrp(1).spike_times_c2') & ~isnan(result_full_nrp(1).event_times_c1) & result_full_nrp(1).spike_times_c2' < result_full_nrp(1).event_times_c1;


figure;
subplot(121)
scatter(result_full_nrp(1).spike_times_c2(paired_trials),...
    result_full_nrp(1).event_times_c1(paired_trials));
hold on
% plot(0:1:200,0:1:200)
% subplot(132); histogram(result_full_nrp(1).event_times_c1(paired_trials))
% hold on
% histogram(result_full_nrp(1).spike_times_c2(paired_trials))

subplot(122);
histogram(result_full_nrp(1).event_times_c1(paired_trials)'...
    - result_full_nrp(1).spike_times_c2(paired_trials))
%%

figure;
unique_powers = unique(result_full_nrp(1).targ_power);
jitter_amt = .5;
for i = 1:length(unique_powers)
    subplot(1,length(unique_powers),i)
    trials = result_full_nrp(1).targ_power == unique_powers(i) & ~isnan(result_full_nrp(1).event_times_c1) & ...
        result_full_nrp(1).event_times_c1 > 50  & result_full_nrp(1).event_times_c1 < 160;
    scatter3(result_full_nrp(1).c2_targ_pos(trials,1)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        result_full_nrp(1).c2_targ_pos(trials,2)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        -result_full_nrp(1).c2_targ_pos(trials,3)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
    80,200 - result_full_nrp(1).event_times_c1(trials),'filled','markerfacealpha',.75)

end
hold on
scatter3(experiment_setup.local_nuc_locs(:,1),experiment_setup.local_nuc_locs(:,2),-experiment_setup.local_nuc_locs(:,3))
%%

figure;

unique_powers = unique(result_full_nrp(1).targ_power);
jitter_amt = .5;
for i = 1:length(unique_powers)
    subplot(1,length(unique_powers),i)
%     subplot(1,length(unique_powers),i)
    trials = result_full_nrp(1).targ_power' == unique_powers(i) & ~isnan(result_full_nrp(1).spike_times_c2);
%     scatter3(experiment_setup.local_nuc_locs(:,1),experiment_setup.local_nuc_locs(:,2)-6,-(experiment_setup.local_nuc_locs(:,3)-6),100,[1 0 0],'filled','markerfacealpha',0.15)
    hold on
    scatter3(result_full_nrp(1).targ_pos(:,1),...
        result_full_nrp(1).targ_pos(:,2),...
        -result_full_nrp(1).targ_pos(:,3),...
        5,[0 0 1],'filled','markerfacealpha',.25)
    hold on
    scatter3(result_full_nrp(1).targ_pos(trials,1)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        result_full_nrp(1).targ_pos(trials,2)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        -result_full_nrp(1).targ_pos(trials,3)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        50,200 - result_full_nrp(1).spike_times_c2(trials),'filled','markerfacealpha',.75)
    xlabel('vertical')
    ylabel('horiz')
    zlabel('axial')
    caxis([0 200])
    axis image
end

%%

% figure;

unique_powers = unique(result_full_nrp(1).targ_power);
jitter_amt = .5;
for i = 1:length(unique_powers)
    subplot(1,length(unique_powers),i)%+length(unique_powers))
%     subplot(1,length(unique_powers),i)
    trials = result_full_nrp(1).targ_power' == unique_powers(i) & ~isnan(result_full_nrp(1).spike_times_c2);
%     scatter3(experiment_setup.local_nuc_locs(:,1),experiment_setup.local_nuc_locs(:,2),-experiment_setup.local_nuc_locs(:,3),100,[1 0 0],'filled','markerfacealpha',0.15)
    hold on
    scatter3(result_full_nrp(1).c2_targ_pos(:,1),...
        result_full_nrp(1).c2_targ_pos(:,2),...
        -result_full_nrp(1).c2_targ_pos(:,3),...
        5,[0 0 1],'filled','markerfacealpha',.25)
    hold on
    scatter3(result_full_nrp(1).c2_targ_pos(trials,1)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        result_full_nrp(1).c2_targ_pos(trials,2)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        -result_full_nrp(1).c2_targ_pos(trials,3)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        50,250 - result_full_nrp(1).event_times_c1(trials),'filled','markerfacealpha',.75)
    caxis([0 200])
    axis image
end


%%
figure
unique_powers = unique(result_full_nrp(1).targ_power);
jitter_amt = .5;
for i = 1:length(unique_powers)
    subplot(1,length(unique_powers),i)%+2*length(unique_powers)
%     subplot(1,length(unique_powers),i)
    trials = result_full_nrp(1).targ_power == unique_powers(i) & result_full_nrp(1).event_times_c1 < 200  & result_full_nrp(1).event_times_c1 > 50;
    scatter3(experiment_setup.local_nuc_locs(:,1),experiment_setup.local_nuc_locs(:,2),-experiment_setup.local_nuc_locs(:,3),125,[1 0 0])
    hold on
%     scatter3(result_full_nrp(1).c2_targ_pos(:,1),...
%         result_full_nrp(1).c2_targ_pos(:,2),...
%         -result_full_nrp(1).c2_targ_pos(:,3),...
%         5,[0 0 1],'filled','markerfacealpha',.25)
%     hold on
    scatter3(result_full_nrp(1).targ_pos(trials,1)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        result_full_nrp(1).targ_pos(trials,2)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        -result_full_nrp(1).targ_pos(trials,3)+rand(sum(trials),1)*jitter_amt-jitter_amt/2,...
        50,200 - result_full_nrp(1).event_times_c1(trials),'filled','markerfacealpha',.75)
    caxis([0 150])
    axis image
end


