%%
     
stacknames = {'/media/shababo/data/01232018images/s2c1-pre - 1_C0'
              '/media/shababo/data/01232018images/s2c1-pre - 2_C0'
              '/media/shababo/data/01232018images/s2c1-pre - 5_C0'
              '/media/shababo/data/01252018images/s2c1-pre_C0'
              '/media/shababo/data/01252018images/s2c1-pre - 1_C0'
              '/media/shababo/data/01252018images/s2c1-pre - 3_C0'
              '/media/shababo/data/01252018images/s2c1-pre - 4_C0'
              '/media/shababo/data/0126272018images/s2c1-pre_C0'
              '/media/shababo/data/0126272018images/s2c1-pre - 1_C0'
              '/media/shababo/data/0126272018images/s2c1-pre - 2_C0'
              '/media/shababo/data/0126272018images/s2c1-pre - 10_C0'
              '/media/shababo/data/0126272018images/s2c1-pre - 11_C0'
              '/media/shababo/data/0126272018images/s2c1-pre - 12_C0'};

filenames = {'1_31_slice1_cell2.mat', '1_31_15_3_data.mat'
             '1_31_slice1_cell3.mat', '1_31_15_13_data.mat'
             '2_21_slice1_cell1.mat', '2_21_16_57_data.mat'
             '2_22_slice1_cell1.mat', '2_22_15_17_data.mat'
             '2_23_slice1_cell1.mat', '2_23_12_15_data.mat'
             '2_23_slice1_cell2.mat', '2_23_12_27_data.mat'
             '2_23_slice2_cell1.mat', '2_23_12_59_data.mat'
             '2_23_slice2_cell2.mat', '2_23_13_24_data.mat'
             '2_23_slice3_cell1.mat', '2_23_14_5_data.mat'
             '2_23_slice3_cell2.mat', '2_23_14_28_data.mat'
             '2_23_slice3_cell3.mat', '2_23_14_40_data.mat'};     
         
spike_trials = {1,1,1,1,1,1,1,1,1,1,1};
current_trials = {4,4,4,4,4,4,4,4,4,4,4};
z_center = [30 30 200 200 200 200 200 200 200 200 200];

colors = jet(length(filenames));
% colors(:) = 0;
% colors(1:size(filenames,1)) = 1; %red

%%

clear result_z
z_pos_vs_cur_and_spike_time = figure;
new_fig = 0;
if ~exist('curr_vs_time','var')
    new_fig = 1;
    curr_vs_time = figure;
end
for j = 1:size(filenames,1)
    
%     clear result_z(j)
    load(filenames{j,2});
    load(filenames{j,1}); 
    
    if ~isnan(spike_trials{j})
        [result_z(j).spike_traces, ~,~, result_z(j).spike_stim_traces] = get_traces(data,spike_trials{j});
        result_z(j).spike_z_pos = [data.trial_metadata(spike_trials{j}).sequence.piezo_z] - z_center(j);
        result_z(j).spike_times = detect_peaks(-bsxfun(@minus,result_z(j).spike_traces,result_z(j).spike_traces(:,1)),20,30,0,Inf,-Inf,0,0,1);
        cell_spike_times = result_z(j).spike_times;
        result_z(j).spike_times = zeros(size(result_z(j).spike_times));
        for i = 1:length(cell_spike_times)
            if ~isempty(cell_spike_times{i})
                result_z(j).spike_times(i) = cell_spike_times{i};
            else
                result_z(j).spike_times(i) = NaN;
            end
        end
    end
    
    if ~isnan(current_trials{j})
        [result_z(j).current_traces, ~,~, result_z(j).current_stim_traces] = get_traces(data,current_trials{j});
        result_z(j).current_z_pos = [data.trial_metadata(current_trials{j}).sequence.piezo_z] - z_center(j);
        result_z(j).max_curr = result_z(j).current_traces(:,1) - min(result_z(j).current_traces,[],2);
    end 
    
    if ~isnan(spike_trials{j}) && ~isnan(current_trials{j})
        tested_pos = unique(result_z(j).spike_z_pos);
        if ~isempty(setdiff(tested_pos,unique(result_z(j).current_z_pos)))
            disp([filenames{j,1} ' spike and current pos sets differ'])
            continue
        end


        for i = 1:length(tested_pos)
            these_trials = result_z(j).spike_z_pos == tested_pos(i);
            result_z(j).spike_time_means(i) = nanmean(result_z(j).spike_times(these_trials));
            result_z(j).spike_time_jitter(i) = nanstd(result_z(j).spike_times(these_trials));
            these_trials = result_z(j).current_z_pos == tested_pos(i);
            result_z(j).max_curr_means(i) = nanmean(result_z(j).max_curr(these_trials));
        end

        figure(curr_vs_time)
        subplot(121)
        [ordered_curr, curr_order] = sort(result_z(j).max_curr_means);
        plot(result_z(j).max_curr_means(curr_order),result_z(j).spike_time_means(curr_order)/20,'-k','LineWidth',2);
        hold on
        if new_fig
            xlabel('Peak Current Mean (pA)')
            ylabel('Spike Time Mean (msec)')
            title('Peak Current vs. Spike Time')
        end
        subplot(122)
        semilogy(result_z(j).max_curr_means(curr_order),result_z(j).spike_time_jitter(curr_order)/20,'-k','LineWidth',2);
        hold on
        if new_fig
            xlabel('Peak Current Mean (pA)')
            ylabel('Spike Time Std. Dev. (msec)')
            title('Peak Current vs. Spike Time Jitter')
        end
    end
    
    
    
    figure(z_pos_vs_cur_and_spike_time)
    if ~isnan(spike_trials{j})
        subplot(121)
        scatter(result_z(j).spike_z_pos,result_z(j).spike_times/20,[],colors(j,:));
        hold on
        plot(tested_pos,result_z(j).spike_time_means/20,'color',colors(j,:))
        xlabel('z distance (um)')
        ylabel('spike time (msec)')
        title('Z Distance vs. Spike Times')
    end
    if ~isnan(current_trials{j})
        subplot(122)
        scatter(result_z(j).current_z_pos,result_z(j).max_curr,[],colors(j,:));
        hold on
        plot(tested_pos,result_z(j).max_curr_means,'color',colors(j,:))
        ylim([0 2500])
        xlabel('z distance (um)')
        ylabel('peak current (pA)')
        title('Z Distance vs. Peak Current')
    end
    
end

