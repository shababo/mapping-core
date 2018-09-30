%%
    

filenames = {'4_26_slice1_cell1.mat', '4_26_17_24_data.mat'
             '4_26_slice2_cell1.mat', '4_26_17_40_data.mat'
             '4_26_slice1_cell3.mat', '4_26_17_53_data.mat'};     
         
spike_trials = {2,[2 3],2};
current_trials = {[],[],[]};

colors = lines(size(filenames,1));


%%


clear result_gain_depth


for j = 1:size(filenames,1)

    
    load(filenames{j,2});
    load(filenames{j,1}); 
    experiment_setup = exp_data.experiment_setup;

    result_gain_depth(j).fluor_val = experiment_setup.fluor_vals(experiment_setup.select_cell_index_full);
    result_gain_depth(j).cell_pos = experiment_setup.patched_cell_loc;

    if ~isempty(spike_trials{j})
        [result_gain_depth(j).spike_traces, ~, this_seq, result_gain_depth(j).spike_stim_traces, stim_key] = get_traces(data,spike_trials{j});
        result_gain_depth(j).spike_targ_pos = [];

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

for j = 1:size(filenames,1)
    
    plot(result_gain_depth(j).spike_unique_powers,result_gain_depth(j).spike_time_means/20,'color',colors(j,:));
    hold on
    scatter(result_gain_depth(j).spike_targ_power,result_gain_depth(j).spike_times/20,[],colors(j,:),'.')
    scatter(result_gain_depth(j).spike_targ_power,result_gain_depth(j).spike_times/20,[],colors(j,:),'o')
    
end

