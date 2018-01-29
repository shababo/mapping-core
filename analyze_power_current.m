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
    
    result_tmp_spike = analyze_current_diffraction_map(data,1,sweep_trial,measure_type);
    switch measure_type
        case 'spikes'
            cell_spike_times = result_tmp.spike_times{1};
            result_tmp.spike_times{1} = zeros(size(result_tmp.spike_times{1}));
            for i = 1:length(cell_spike_times)
                if ~isempty(cell_spike_times{i})
                    result_tmp.spike_times{1}(i) = cell_spike_times{i};
                else
                    result_tmp.spike_times{1}(i) = NaN;
                end
            end
    end
    result_tmp.power = {zeros(size(result_tmp.stim_size{1}))};
    
    unique_powers = unique(round(result_tmp.stim_size{1},2));
    unique_powers(unique_powers > .49) = [];
%     unique_powers = sort(unique(round(trial_powers)));
    result_tmp.these_powers = unique_powers;
    result_tmp.power_means = zeros(size(unique_powers))';
    

    for i = 1:length(unique_powers)
        these_trials = find(abs(result_tmp.stim_size{1} - unique_powers(i)) < .005);
        result_tmp.power{1}(these_trials) = pockels_map_fine(abs(pockels_map_fine(:,1) - unique_powers(i)) < .004,2);

%         upper_ind = find(pockels_map_fine(:,1) > unique_powers(i),1,'first');
%         if isempty(upper_ind)
%             upper_ind = size(pockels_map_fine,1);
%         end
%         lower_ind = upper_ind - 1;
%         step_size_volt = pockels_map_fine(upper_ind,1) - pockels_map_fine(lower_ind,1);
%         interp_dist = (unique_powers(i) - pockels_map_fine(lower_ind,1))/step_size_volt;
%         step_size_power = pockels_map_fine(upper_ind,2) - pockels_map_fine(lower_ind,2);
%         result_tmp.power{1}(these_trials) = pockels_map_fine(lower_ind,2) + interp_dist*step_size_power;
        
        switch measure_type
            case 'current'
                result_tmp.power_means(i) = mean(result_tmp.max_curr{1}(these_trials));
            case 'spikes'
                result_tmp.power_means(i) = nanmean(result_tmp.spike_times{1}(these_trials));
                result_tmp.power_ranges(i) = max(result_tmp.spike_times{1}(these_trials)) - min(result_tmp.spike_times{1}(these_trials));
                result_tmp.power_jitter(i) = nanstd(result_tmp.spike_times{1}(these_trials));
                result_tmp.prob_spike(i) = sum(~isnan(result_tmp.spike_times{1}(these_trials)))/length(these_trials);
                
        end
    end
%         
    
    
%     subplot(211)
    switch measure_type
        case 'current'
            figure(current_figure)
            gca
            scatter(result_tmp.power{1}.^1,result_tmp.max_curr{1},15,'jitter','on','jitteramount',.3);
            result_current(j) = result_tmp;
        case 'spikes'
            figure(spike_figure)
            gca
            scatter(result_tmp.power{1}.^1,result_tmp.spike_times{1}/20,15,'jitter','on','jitteramount',.3);
            result_spikes(j) = result_tmp;
    end
    
    hold on
    
    
%     subplot(212)
%     scatter(result_tmp.stim_size{1},result_tmp.max_curr{1},5,colors{j},'filled');
% %     plot([0; these_powers],[0, result_tmp.power_means(result_tmp.power{1} ~= 0)]);
%     hold on

%     
%     [sorted_power,order_power] = sort(disk_power_meas_new_post(subset_pos));
%     [sorted_fluor2D,order_fluor2D] = sort(disk_fluor_meas2D(subset_pos));
%     [sorted_fluor3D,order_fluor3D] = sort(disk_fluor_meas3D(subset_pos));
% 
%     subplot(311)
% 	plot(sorted_power,this_cell_cur(subset_pos(order_power)),'-','color',colors{j})
%     title('Power Measure vs. Current')
%     ylabel('Peak Current (pA)')
%     xlabel('Power (mW)')
%     hold on
% %     subplot(312)
% % 	plot(sorted_fluor2D,this_cell_cur(subset_pos(order_fluor2D)),'-','color',colors{j})
% %     hold on
%     subplot(312)
% 	plot(sorted_fluor3D,this_cell_cur(subset_pos(order_fluor3D)),'-','color',colors{j})
%     title('Fluor Measure vs. Current')
%     ylabel('Peak Current (pA)')
%     xlabel('Fluor (a.u.)')
%     hold on
%     subplot(313)
% 	plot(sorted_power.^2,this_cell_cur(subset_pos(order_power)),'-','color',colors{j})
%     hold on
%     title('Power Measure Squared vs. Current')
%     ylabel('Peak Current (pA)')
%     xlabel('Power Measure Squared (mW^2)')
end
switch measure_type
    case 'current'
        title('Power vs. Peak Current')
        xlabel('Power (mW)')
        xlim([10 90])
        ylabel('Peak Current (pA)')
        ylim([0 2500])

    case 'spikes'
        title('Power vs. Spike Times')
        xlabel('Power (mW)')
        xlim([10 90])
        ylabel('Spike Times (msec)')
        ylim([0 15])
end

% legend({'cell 1','cell 2','cell 3'})

%%

% plot(26*ones(size(0:3000)),0:3000,'g--')
% plot(25.5*ones(size(0:3000)),0:3000,'b--')
plot(0:85,522.5*ones(size(0:85)),'b--')
plot(0:85,361.8*ones(size(0:85)),'g--')
plot(0:85,678.5*ones(size(0:85)),'r--')
% plot(14.9*ones(size(0:3000)),0:3000,'r--')

%%

figure(both_figure)
for i = 1:length(result_current)
    gca
    [shared_powers, current_i, spikes_i] = intersect(result_current(i).these_powers,result_spikes(i).these_powers);
    plot(result_current(i).power_means(current_i),result_spikes(i).power_means(spikes_i)/20,'o-')
    hold on
end

title('peak current vs. spike time range')
xlabel('mean peak current (pA)')
ylabel('spike time range (msec)')
xlim([0 2500])
ylim([0 15])