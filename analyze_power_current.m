% function result_tmp = analyze_power_current(data)

%%

filenames = {'1_23_slice1_cell1.mat'
             '1_23_slice1_cell2.mat'
             '1_23_slice2_cell2.mat'
             '1_25_slice1_cell1.mat'
             '1_25_slice1_cell3.mat'
             '1_25_slice2_cell1.mat'
             '1_25_slice2_cell2.mat'
             '1_26_slice1_cell1.mat'
             '1_26_slice1_cell2.mat'
             '1_26_slice1_cell3.mat'
             '1_27_slice1_cell1.mat'
             '1_27_slice1_cell2.mat'
             '1_27_slice1_cell3.mat'};        
        
spike_trial_inds = {1,1,1,1,1,1,1,[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
current_trial_inds = {4,4,4,4,4,4,4,5,5,5,5,5,4};

% pos_order = [1 2 4 3 5 6
%              3 2 5 1 4 6
%              6 1 3 4 5 2
%              4 6 3 5 2 1];

% colors = {'r','b','g','k','c','y','m'};


spike_figure = figure;
current_figure = figure;
both_figure = figure;
%% 

for j = 1:size(filenames,1)
    
    
    
    load(filenames{j}); 
    
    result_tmp_spike = analyze_current_diffraction_map(data,1,spike_trial_inds(j),'spikes');
    result_tmp_current = analyze_current_diffraction_map(data,1,current_trial_inds(j),'current');

    cell_spike_times = result_tmp_spike.spike_times{1};
    result_tmp_spike.spike_times{1} = zeros(size(result_tmp_spike.spike_times{1}));
    for i = 1:length(cell_spike_times)
        if ~isempty(cell_spike_times{i})
            result_tmp_spike.spike_times{1}(i) = cell_spike_times{i};
        else
            result_tmp_spike.spike_times{1}(i) = NaN;
        end
    end
    
    result_tmp_spike.power = {zeros(size(result_tmp_spike.stim_size{1}))};
    result_tmp_current.power = {zeros(size(result_tmp_current.stim_size{1}))};
    
    unique_powers_spike = unique(round(result_tmp_spike.stim_size{1},2));
    result_tmp_spike.these_powers = unique_powers_spike;
    result_tmp_spike.power_means = zeros(size(unique_powers_spike))';
    
    unique_powers_current = unique(round(result_tmp_current.stim_size{1},2));
    result_tmp_current.these_powers = unique_powers_current;
    result_tmp_current.power_means = zeros(size(unique_powers_current))';
    

    for i = 1:length(unique_powers_spike)
        these_trials = find(abs(result_tmp_spike.stim_size{1} - unique_powers_spike(i)) < .005);
        result_tmp_spike.power{1}(these_trials) = pockels_map_fine(abs(pockels_map_fine(:,1) - unique_powers_spike(i)) < .004,2);

%         upper_ind = find(pockels_map_fine(:,1) > unique_powers(i),1,'first');
%         if isempty(upper_ind)
%             upper_ind = size(pockels_map_fine,1);
%         end
%         lower_ind = upper_ind - 1;
%         step_size_volt = pockels_map_fine(upper_ind,1) - pockels_map_fine(lower_ind,1);
%         interp_dist = (unique_powers(i) - pockels_map_fine(lower_ind,1))/step_size_volt;
%         step_size_power = pockels_map_fine(upper_ind,2) - pockels_map_fine(lower_ind,2);
%         result_tmp.power{1}(these_trials) = pockels_map_fine(lower_ind,2) + interp_dist*step_size_power;
        

                result_tmp_spike.power_means(i) = nanmean(result_tmp_spike.spike_times{1}(these_trials));
                result_tmp_spike.power_ranges(i) = max(result_tmp_spike.spike_times{1}(these_trials)) - min(result_tmp_spike.spike_times{1}(these_trials));
                result_tmp_spike.power_jitter(i) = nanstd(result_tmp_spike.spike_times{1}(these_trials));
                result_tmp_spike.prob_spike(i) = sum(~isnan(result_tmp_spike.spike_times{1}(these_trials)))/length(these_trials);
                
        
    end
    
    for i = 1:length(unique_powers_current)
        these_trials = find(abs(result_tmp_current.stim_size{1} - unique_powers_current(i)) < .005);
        result_tmp_current.power{1}(these_trials) = pockels_map_fine(abs(pockels_map_fine(:,1) - unique_powers_current(i)) < .004,2);

%         upper_ind = find(pockels_map_fine(:,1) > unique_powers(i),1,'first');
%         if isempty(upper_ind)
%             upper_ind = size(pockels_map_fine,1);
%         end
%         lower_ind = upper_ind - 1;
%         step_size_volt = pockels_map_fine(upper_ind,1) - pockels_map_fine(lower_ind,1);
%         interp_dist = (unique_powers(i) - pockels_map_fine(lower_ind,1))/step_size_volt;
%         step_size_power = pockels_map_fine(upper_ind,2) - pockels_map_fine(lower_ind,2);
%         result_tmp.power{1}(these_trials) = pockels_map_fine(lower_ind,2) + interp_dist*step_size_power;
        

        result_tmp_current.power_means(i) = mean(result_tmp_current.max_curr{1}(these_trials));

    end
%         
    
    
%     subplot(211)

            figure(current_figure)
            gca
            hold on
            scatter(result_tmp_current.power{1}.^1,result_tmp_current.max_curr{1},15,'jitter','on','jitteramount',.3);
            result_current(j) = result_tmp_current;
   
            figure(spike_figure)
            gca
            hold on
            scatter(result_tmp_spike.power{1}.^1,result_tmp_spike.spike_times{1}/20,15,'jitter','on','jitteramount',.3);
            result_spikes(j) = result_tmp_spike;
    
    
    hold on
    
    
end

        figure(current_figure)
            gca
        title('Power vs. Peak Current')
        xlabel('Power (mW)')
        xlim([10 90])
        ylabel('Peak Current (pA)')
        ylim([0 2500])

        figure(spike_figure)
            gca
        title('Power vs. Spike Times')
        xlabel('Power (mW)')
        xlim([10 90])
        ylabel('Spike Times (msec)')
        ylim([0 15])


% legend({'cell 1','cell 2','cell 3'})

%%

% plot(26*ones(size(0:3000)),0:3000,'g--')
% plot(25.5*ones(size(0:3000)),0:3000,'b--')
plot(0:85,522.5*ones(size(0:85)),'b--')
plot(0:85,361.8*ones(size(0:85)),'g--')
plot(0:85,678.5*ones(size(0:85)),'r--')
% plot(14.9*ones(size(0:3000)),0:3000,'r--')

%%
figure
% figure(both_figure)
for i = 1:length(result_current)
    gca
    [shared_powers, current_i, spikes_i] = intersect(result_current(i).these_powers,result_spikes(i).these_powers);
    plot(result_current(i).power_means(current_i),result_spikes(i).prob_spike(spikes_i),'o')
    hold on
end

title('peak current vs. p(spike)')
xlabel('mean peak current (pA)')
ylabel('prob spike')
xlim([0 2500])
% ylim([0 15])