filenames = {'1_5_slice2_cell1.mat','1_5_18_8_data.mat'
            '1_5_slice2_cell2.mat','1_5_18_49_data.mat'
            '1_5_slice2_cell3.mat','1_5_19_8_data.mat'
            '1_5_slice3_cell1.mat','1_5_19_39_data.mat'
            '1_5_slice3_cell2.mat','1_5_20_0_data.mat'};
        
trial_inds = {3:8,3:8,3:8,4:8,3:6};

colors = {'r','b','g','k','c'};

%%

filenames = {'1_19_slice1_cell1.mat','1_19_17_46_data.mat'
            '1_19_slice1_cell2.mat','1_19_18_9_data.mat'
            '1_19_slice2_cell1.mat','1_19_18_32_data.mat'
            '1_19_slice2_cell2.mat','1_19_18_47_data.mat'
            '1_19_slice3_cell1.mat','1_19_19_7_data.mat'};
        
trial_inds = {2:7,3:8,3:8,3:8,3:8};

pos_order = [1 2 4 3 5 6
             3 2 5 1 4 6
             6 1 3 4 5 2
             4 6 3 5 2 1];

colors = {'r','b','g','k','c'};
%%

filenames = {'1_19_slice1_cell1.mat','1_19_17_46_data.mat'
            '1_19_slice1_cell2.mat','1_19_18_9_data.mat'
            '1_19_slice2_cell1.mat','1_19_18_32_data.mat'
            '1_19_slice2_cell2.mat','1_19_18_47_data.mat'
            '1_19_slice3_cell1.mat','1_19_19_7_data.mat'};
        
trial_inds = {2:7,3:8,3:8,3:8,3:8};

% pos_order = [1 2 4 3 5 6
%              3 2 5 1 4 6
%              6 1 3 4 5 2
%              4 6 3 5 2 1];

colors = {'r','b','g','k','c'};
%%

filenames = {'1_21_slice1_cell1.mat','1_21_12_32_data.mat'
            '1_21_slice1_cell2.mat','1_21_12_56_data.mat'
            '1_21_slice2_cell1.mat','1_21_13_24_data.mat'};
        
trial_inds = {3:8,3:8,3:8};

% pos_order = [1 2 4 3 5 6
%              3 2 5 1 4 6
%              6 1 3 4 5 2
%              4 6 3 5 2 1];

colors = {'r','b','g','k','c'};

%%

filenames = {'1_21_slice1_cell1.mat','1_21_12_32_data.mat'
            '1_21_slice1_cell2.mat','1_21_12_56_data.mat'
            '1_21_slice2_cell1.mat','1_21_13_24_data.mat'};
        
trial_inds = {3:8,3:8,3:8};

% pos_order = [1 2 4 3 5 6
%              3 2 5 1 4 6
%              6 1 3 4 5 2
%              4 6 3 5 2 1];

colors = {'r','b','g','k','c'};

%%

filenames = {'1_22_slice2_cell1.mat','1_22_14_34_data.mat'
            '1_22_slice2_cell2.mat','1_22_14_49_data.mat'
            'tmp-2018122155.mat','1_22_15_5_data.mat'
            '1_22_slice3_cell1.mat','1_22_15_41_data.mat'
            '1_22_slice3_cell2.mat','1_22_15_58_data.mat'};
        
trial_inds = {4:9,2:7,3:8,3:8,3:8};

% pos_order = [1 2 4 3 5 6
%              3 2 5 1 4 6
%              6 1 3 4 5 2
%              4 6 3 5 2 1];

colors = {'r','b','g','k','c'};



%%

figure


    
for j = 1:size(filenames,1)
    
    these_trials = trial_inds{j};
    
    load(filenames{j,2}); load(filenames{j,1}); 
    
    j
    pos_order(j,:) = experiment_setup.pos_order
%     
    
    result = analyze_current_diffraction_map(data,pos_order(j,:),these_trials);
    this_cell_cur = zeros(size(these_trials));
    
    subset_pos = sort(pos_order(j,1:6));
%     subset_pos = 1:3;
    for i = 1:length(subset_pos)
        pos_ind = subset_pos(i);
        these_cur = result.max_curr{pos_ind}(abs(result.stim_size{pos_ind} - .25) < .05);
        this_cell_cur(pos_ind) = mean(these_cur);
        subplot(311)
        scatter(disk_power_meas_new_post(pos_ind)*ones(size(these_cur)),these_cur,5,colors{j},'filled');
        hold on
%         subplot(312)
%         scatter(disk_fluor_meas2D(pos_ind)*ones(size(these_cur)),these_cur,5,colors{j},'filled');
%         hold on
        subplot(312)
        scatter(disk_fluor_meas3D(pos_ind)*ones(size(these_cur)),these_cur,5,colors{j},'filled');
        hold on
        subplot(313)
        scatter(disk_power_meas_new_post(pos_ind).^2*ones(size(these_cur)),these_cur,5,colors{j},'filled');
        hold on
%         scatter(result.stim_size{i},result.max_curr{i},10,colors{j}); 
%         hold on
    end
    
    [sorted_power,order_power] = sort(disk_power_meas_new_post(subset_pos));
    [sorted_fluor2D,order_fluor2D] = sort(disk_fluor_meas2D(subset_pos));
    [sorted_fluor3D,order_fluor3D] = sort(disk_fluor_meas3D(subset_pos));

    subplot(311)
	plot(sorted_power,this_cell_cur(subset_pos(order_power)),'-','color',colors{j})
    title('Power Measure vs. Current')
    ylabel('Peak Current (pA)')
    xlabel('Power (mW)')
    hold on
%     subplot(312)
% 	plot(sorted_fluor2D,this_cell_cur(subset_pos(order_fluor2D)),'-','color',colors{j})
%     hold on
    subplot(312)
	plot(sorted_fluor3D,this_cell_cur(subset_pos(order_fluor3D)),'-','color',colors{j})
    title('Fluor Measure vs. Current')
    ylabel('Peak Current (pA)')
    xlabel('Fluor (a.u.)')
    hold on
    subplot(313)
	plot(sorted_power.^2,this_cell_cur(subset_pos(order_power)),'-','color',colors{j})
    hold on
    title('Power Measure Squared vs. Current')
    ylabel('Peak Current (pA)')
    xlabel('Power Measure Squared (mW^2)')
end

figure

 subplot(311)
plot(sorted_power/max(sorted_power),this_cell_cur(subset_pos(order_power))/max(this_cell_cur(subset_pos)),'-','color',colors{j})
title('Power Measure vs. Current')
ylabel('Peak Current (pA)')
xlabel('Power (mW)')
hold on
%     subplot(312)
% 	plot(sorted_fluor2D,this_cell_cur(subset_pos(order_fluor2D)),'-','color',colors{j})
%     hold on
subplot(312)
plot(sorted_fluor3D/max(sorted_fluor3D),this_cell_cur(subset_pos(order_fluor3D))/max(this_cell_cur(subset_pos)),'-','color',colors{j})
title('Fluor Measure vs. Current')
ylabel('Peak Current (pA)')
xlabel('Fluor (a.u.)')
hold on
subplot(313)
plot((sorted_power.^2)/max((sorted_power.^2)),this_cell_cur(subset_pos(order_power)/max(this_cell_cur(subset_pos))),'-','color',colors{j})
hold on
title('Power Measure Squared vs. Current')
ylabel('Peak Current (pA)')
xlabel('Power Measure Squared (mW^2)')

