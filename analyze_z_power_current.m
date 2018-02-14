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
             '2_9_slice2_cell3.mat', '2_9_17_54_data.mat'};     
         
spike_trials = {1,1};
current_trials = {4,4};
         
% z_pos = [-8 0 8];
%%



for j = 1
    
    clear result(j)
    load(filenames{j,2});
    load(filenames{j,1}); 
    
    [result(j).spike_traces, ~,~, result(j).spike_stim_traces] = get_traces(data,spike_trials{j});
    result(j).spike_z_pos = [data.trial_metadata(spike_trials{j}).sequence.piezo_z];
    [result(j).current_traces, ~,~, result(j).current_stim_traces] = get_traces(data,current_trials{j});
    result(j).current_z_pos = [data.trial_metadata(current_trials{j}).sequence.piezo_z];
    result(j).max_curr = result(j).current_traces(:,1) - min(result(j).current_traces,[],2);
    result(j).spike_times = detect_peaks(-bsxfun(@minus,result(j).spike_traces,result(j).spike_traces(:,1)),20,30,0,Inf,-Inf,0,0,1);
    cell_spike_times = result(j).spike_times;
    result(j).spike_times = zeros(size(result(j).spike_times));
    for i = 1:length(cell_spike_times)
        if ~isempty(cell_spike_times{i})
            result(j).spike_times(i) = cell_spike_times{i};
        else
            result(j).spike_times(i) = NaN;
        end
    end
    
    figure;
    subplot(121)
    scatter(result(j).spike_z_pos,result(j).spike_times);
    subplot(122)
    scatter(result(j).current_z_pos,result(j).max_curr);
    
    
end

% colors = {'r','b','g','k','c','y','m'};
colors = jet(length(filenames));

spike_figure = figure;
current_figure = figure;
both_figure = figure;