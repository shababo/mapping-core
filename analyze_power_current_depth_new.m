77777777%%
    

filenames = {'7_19_slice1_cell1.mat', '7_19_13_6_data.mat'
%              '7_19_slice1_cell2.mat', '7_19_13_24_data.mat'
             '7_19_slice1_cell4.mat', '7_19_13_47_data.mat'
             '7_19_slice1_cell5.mat', '7_19_13_51_data.mat'
             '7_19_slice2_cell3.mat', '7_19_15_0_data.mat'
             '7_19_slice2_cell4.mat', '7_19_15_17_data.mat'
             '7_19_slice3_cell1.mat', '7_19_15_28_data.mat'
             '7_19_slice3_cell2.mat', '7_19_15_36_data.mat'
             '7_19_slice3_cell4.mat', '7_19_15_53_data.mat'
             '7_19_slice3_cell5.mat', '7_19_15_57_data.mat'
             '7_19_slice3_cell6.mat', '7_19_16_5_data.mat'
             '7_19_slice3_cell7.mat', '7_19_16_11_data.mat'
             };     
         
spike_trials = {2,2,2,2,3,2,2,2,2,2,3};

current_trials = {[],[],[],[],[],[],[],[],[],[],[]};

colors = lines(size(filenames,1));


%%


% clear result_gain_depth


for j = 1:11
    

    
    load(filenames{j,2});
    load(filenames{j,1}); 
    experiment_setup = exp_data.experiment_setup;
    
    if ~isfield(experiment_setup, 'patched_cell_fluor')
        disp('fixing detection')
        experiment_setup.select_cell_index_full = exp_data.stack_viewer_output.selected_cell_index_full;
        do_detect = 0;
        [experiment_setup.nuclear_locs,experiment_setup.fluor_vals, ...
            experiment_setup.nuclear_locs_image_coord] = ...
                        detect_nuclei([experiment_setup.exp_id '_stack'],[],[125; 119],[],do_detect,[],0);
            experiment_setup.patched_cell_pos = exp_data.stack_viewer_output.selected_cell_pos;
            
            experiment_setup.patched_cell_loc = experiment_setup.nuclear_locs(experiment_setup.select_cell_index_full,:);
            experiment_setup.patched_cell_fluor = experiment_setup.fluor_vals(experiment_setup.select_cell_index_full);
            save([filenames{j,2}(1:end-4) '_bu.mat'],'exp_data')
            exp_data.experiment_setup = experiment_setup;
            save(filenames{j,2},'exp_data')
            
    end
    
    result_gain_depth(j).fluor_val = experiment_setup.patched_cell_fluor;
    result_gain_depth(j).cell_pos = experiment_setup.patched_cell_loc;
    
    [tissue_depth, plane] = fit_tissue_surface(experiment_setup.nuclear_locs,50);
    result_gain_depth(j).cell_tissue_depth = tissue_depth(experiment_setup.select_cell_index_full);

    if ~isempty(spike_trials{j})
        [result_gain_depth(j).spike_traces, ~, this_seq, result_gain_depth(j).spike_stim_traces, stim_key] = get_traces(data,spike_trials{j});
%         result_gain_depth(j).spike_targ_pos = [];

        result_gain_depth(j).spike_targ_pos = ...
            bsxfun(@minus,stim_key([this_seq.precomputed_target_index],:),result_gain_depth(j).cell_pos);
        result_gain_depth(j).spike_targ_power =[this_seq.target_power]';
        
        result_gain_depth(j).spike_times = detect_peaks(-bsxfun(@minus,result_gain_depth(j).spike_traces,result_gain_depth(j).spike_traces(:,1)),75,30,0,Inf,-Inf,0,0,1);
        cell_spike_times = result_gain_depth(j).spike_times;
        result_gain_depth(j).spike_times = zeros(size(result_gain_depth(j).spike_times));

        for i = 1:length(cell_spike_times)
            if ~isempty(cell_spike_times{i})
                result_gain_depth(j).spike_times(i) = cell_spike_times{i};
            else
                result_gain_depth(j).spike_times(i) = NaN;
            end
        end
        
        result_gain_depth(j).spike_unique_powers = unique(result_gain_depth(j).spike_targ_power);
        
        for i = 1:length(result_gain_depth(j).spike_unique_powers)

            these_trials = result_gain_depth(j).spike_targ_power == result_gain_depth(j).spike_unique_powers(i);
            result_gain_depth(j).spike_time_means(i) = nanmean(result_gain_depth(j).spike_times(these_trials));
            result_gain_depth(j).spike_time_jitter(i) = nanstd(result_gain_depth(j).spike_times(these_trials));
            result_gain_depth(j).prob_spike(i) = sum(~isnan(result_gain_depth(j).spike_times(these_trials)))/sum(these_trials);
        end
    
    end
    
    if ~isempty(current_trials{j})   
       
        [result_gain_depth(j).current_traces, ~,this_seq, result_gain_depth(j).current_stim_traces] = get_traces(data,current_trials{j});
        result_gain_depth(j).current_targ_pos = ...
            bsxfun(@minus,data.trial_metadata(current_trials{j}).stim_key([this_seq.precomputed_target_index],:),result_gain_depth(j).cell_pos);
        result_gain_depth(j).current_targ_power =[this_seq.target_power];
        result_gain_depth(j).max_curr = result_gain_depth(j).current_traces(:,1) - min(result_gain_depth(j).current_traces,[],2);
        
        result_gain_depth(j).curr_unique_powers = unique(result_gain_depth(j).current_targ_power);
        
        for i = 1:length(result_gain_depth(j).curr_unique_powers)

            these_trials = result_gain_depth(j).spike_targ_power == result_gain_depth(j).curr_unique_powers(i);
            result_gain_depth(j).max_curr_means(i) = nanmean(result_gain_depth(j).max_curr(these_trials));
            result_gain_depth(j).max_curr_std(i) = nanstd(result_gain_depth(j).max_curr(these_trials));
            
        end

    end
    

end

%%

figure;
norm_fluor = [result_gain_depth.fluor_val] - min([result_gain_depth.fluor_val])
norm_fluor = norm_fluor/max(norm_fluor)
for j = 1:size(filenames,1)
    
    plot(result_gain_depth(j).spike_unique_powers,result_gain_depth(j).spike_time_means/20,'color',[norm_fluor(j) 0 1-norm_fluor(j)]);
    hold on
%     scatter(result_gain_depth(j).spike_targ_power,result_gain_depth(j).spike_times/20,[],colors(j,:),'.')
%     scatter(result_gain_depth(j).spike_targ_power,result_gain_depth(j).spike_times/20,[],colors(j,:),'o')
    
end

%%



for j = 1:size(filenames,1)
    
    result_gain_depth(j).proxy_gain = -result_gain_depth(j).spike_time_means(3)/20;
    
end

%%
figure
scatter3([result_gain_depth.cell_tissue_depth],[result_gain_depth.fluor_val],[result_gain_depth.proxy_gain])
xlabel('tissue depth'); ylabel('fluor val'); zlabel('proxy gain')
% zlim([-120 10])

%%
figure
plotmatrix([[result_gain_depth.cell_tissue_depth];[result_gain_depth.fluor_val];[gain_mle(1:11)]]')
%%

%figure
hold on
scatter([result_gain_depth.cell_tissue_depth],[result_gain_depth.fluor_val],'rx')
xlabel('tissue depth'); ylabel('fluor val'); 

%% imaging calibration stuff

pockels_vals = {'026','050','076','100','126'};

for i = 3
    [pockels_test(i).nuclear_locs,pockels_test(i).fluor_vals,pockels_test(i).nuclear_locs_image_coord] = ...
        detect_nuclei(['/media/shababo/data/pockels_' pockels_vals{i} '_C0'],[],[],[],1,[],0);
end

%%
figure

subplot(121)
histogram(pockels_test(4).fluor_vals,0:20:1000);
hold on
histogram(pockels_test(5).fluor_vals,0:20:1000);

subplot(122)
histogram(pockels_test(4).nuclear_locs(:,3));
hold on
histogram(pockels_test(5).nuclear_locs(:,3));

%%

figure
scatter(pockels_test(4).nuclear_locs(:,3),pockels_test(4).fluor_vals)
hold on
scatter(pockels_test(5).nuclear_locs(:,3),pockels_test(5).fluor_vals)

%%
figure;
scatter3(pockels_test(4).nuclear_locs(:,1),pockels_test(4).nuclear_locs(:,2),-pockels_test(4).nuclear_locs(:,3),'bo')
hold on
scatter3(pockels_test(5).nuclear_locs(:,1),pockels_test(5).nuclear_locs(:,2),-pockels_test(5).nuclear_locs(:,3),'rx')

