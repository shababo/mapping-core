% function result_tmp = analyze_power_current(data)

%%

filenames = {'1_23_slice1_cell1.mat', '/media/shababo/data/01232018images/s2c1-pre - 1_C0'
             '1_23_slice1_cell2.mat', '/media/shababo/data/01232018images/s2c1-pre - 2_C0'
             '1_23_slice2_cell1.mat', '/media/shababo/data/01232018images/s2c1-pre - 4_C0' % spkes 1 and currents 4
             '1_23_slice2_cell2.mat', '/media/shababo/data/01232018images/s2c1-pre - 5_C0'
             '1_25_slice1_cell1.mat', ''
             '1_25_slice1_cell2.mat', '/media/shababo/data/01252018images/s2c1-pre_C0' % spikes 1 only
             '1_25_slice1_cell3.mat', '/media/shababo/data/01252018images/s2c1-pre - 1_C0'
             '1_25_slice2_cell1.mat', '/media/shababo/data/01252018images/s2c1-pre - 3_C0'
             '1_25_slice2_cell2.mat', '/media/shababo/data/01252018images/s2c1-pre - 4_C0' % MAYBE THIS CELL'S IMAGE??
             '1_25_slice2_cell3x.mat', '/media/shababo/data/01252018images/s2c1-pre - 6_C0' % spikes 1 only
             '1_25_slice2_cell4.mat', '/media/shababo/data/01252018images/s2c1-pre - 7_C0'  % spikes 1 only
             '1_25_slice3_cell1.mat', '/media/shababo/data/01252018images/s2c1-pre - 8_C0'  % spikes 1, kinda shitty currents 3
             '1_25_slice3_cell1next.mat', '/media/shababo/data/01252018images/s2c1-pre - 9_C0'  % spikes 1
             '1_26_slice1_cell1.mat', '/media/shababo/data/0126272018images/s2c1-pre_C0'
             '1_26_slice1_cell2.mat', '/media/shababo/data/0126272018images/s2c1-pre - 1_C0'
             '1_26_slice1_cell3.mat', '/media/shababo/data/0126272018images/s2c1-pre - 2_C0'
             '1_26_slice2_cell2.mat', '/media/shababo/data/0126272018images/s2c1-pre - 4_C0'  % spikes [1 2] only
             '1_26_slice2_cell3.mat', '/media/shababo/data/0126272018images/s2c1-pre - 6_C0'  % spikes [1 2] only
             '1_26_slice3_cell1next.mat', '/media/shababo/data/0126272018images/s2c1-pre - 7_C0'  % spikes [1 2] only
             '1_26_slice3_cell2next.mat', '/media/shababo/data/0126272018images/s2c1-pre - 8_C0'  % spikes [1 2] only
             '1_26_slice3_cell3.mat', '/media/shababo/data/0126272018images/s2c1-pre - 9_C0'  % spikes [1 2] only
             '1_27_slice1_cell1.mat', '/media/shababo/data/0126272018images/s2c1-pre - 10_C0'
             '1_27_slice1_cell2.mat', '/media/shababo/data/0126272018images/s2c1-pre - 11_C0'
             '1_27_slice1_cell3.mat', '/media/shababo/data/0126272018images/s2c1-pre - 12_C0' % bad vm/ih
             '1_27_slice2_cell1.mat', '/media/shababo/data/0126272018images/s2c1-pre - 13_C0'  % spikes [1 2] only 
             '1_27_slice2_cell2.mat', '/media/shababo/data/0126272018images/s2c1-pre - 14_C0'  % spikes [1 2], current 4 - bad holding
             '1_27_slice3_cell1.mat', '/media/shababo/data/0126272018images/s2c1-pre - 16_C0'  % spikes [1 2] only
             '1_27_slice3_cell2.mat', '/media/shababo/data/0126272018images/s2c1-pre - 17_C0'  % spikes [1 2] and current 5 - no evoked spikes
             };


cell_pos = [118 135 10];
        
spike_trial_inds = {1,1,1,1,1,1,1,1,1,1,1,1,1,...
    [1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
current_trial_inds = {4,4,4,4,4,[],4,4,4,[],[],3,[],...
    5,5,5,[],[],[],[],[],5,5,4,[],4,[],5};


% pos_order = [1 2 4 3 5 6
%              3 2 5 1 4 6
%              6 1 3 4 5 2
%              4 6 3 5 2 1]z

% colors = {'r','b','g','k','c','y','m'};

colors = jet(size(filenames,1));

% spike_figure = figure;
% current_figure = figure;
sumary_fig = figure;

do_detect = 1;

%%
clear result_current result_spikes
%% 

for j = 1:size(filenames,1)
    
    
    
    load(filenames{j,1}); 
    
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
    result_tmp_spike.spike_time_means = zeros(size(unique_powers_spike))';
    
    unique_powers_current = unique(round(result_tmp_current.stim_size{1},2));
    result_tmp_current.these_powers = unique_powers_current;
    result_tmp_current.peak_current_means = zeros(size(unique_powers_current))';
    

    for i = 1:length(unique_powers_spike)
        these_trials = find(abs(result_tmp_spike.stim_size{1} - unique_powers_spike(i)) < .005);
        result_tmp_spike.power{1}(these_trials) = pockels_map_fine(abs(pockels_map_fine(:,1) - unique_powers_spike(i)) < .004,2);
        unique_powers_spike(i) = pockels_map_fine(abs(pockels_map_fine(:,1) - unique_powers_spike(i)) < .004,2);
%         upper_ind = find(pockels_map_fine(:,1) > unique_powers(i),1,'first');
%         if isempty(upper_ind)
%             upper_ind = size(pockels_map_fine,1);
%         end
%         lower_ind = upper_ind - 1;
%         step_size_volt = pockels_map_fine(upper_ind,1) - pockels_map_fine(lower_ind,1);
%         interp_dist = (unique_powers(i) - pockels_map_fine(lower_ind,1))/step_size_volt;
%         step_size_power = pockels_map_fine(upper_ind,2) - pockels_map_fine(lower_ind,2);
%         result_tmp.power{1}(these_trials) = pockels_map_fine(lower_ind,2) + interp_dist*step_size_power;
        

                result_tmp_spike.spike_time_means(i) = nanmean(result_tmp_spike.spike_times{1}(these_trials));
                result_tmp_spike.spike_time_ranges(i) = max(result_tmp_spike.spike_times{1}(these_trials)) - min(result_tmp_spike.spike_times{1}(these_trials));
                result_tmp_spike.spike_time_jitter(i) = nanstd(result_tmp_spike.spike_times{1}(these_trials));
                result_tmp_spike.prob_spike(i) = sum(~isnan(result_tmp_spike.spike_times{1}(these_trials)))/length(these_trials);
                
        
    end
    
    for i = 1:length(unique_powers_current)
        these_trials = abs(result_tmp_current.stim_size{1} - unique_powers_current(i)) < .005;
        unique_powers_current(i) = pockels_map_fine(abs(pockels_map_fine(:,1) - unique_powers_current(i)) < .004,2);
        result_tmp_current.power{1}(these_trials) = unique_powers_current(i);
        
%         upper_ind = find(pockels_map_fine(:,1) > unique_powers(i),1,'first');
%         if isempty(upper_ind)
%             upper_ind = size(pockels_map_fine,1);
%         end
%         lower_ind = upper_ind - 1;
%         step_size_volt = pockels_map_fine(upper_ind,1) - pockels_map_fine(lower_ind,1);
%         interp_dist = (unique_powers(i) - pockels_map_fine(lower_ind,1))/step_size_volt;
%         step_size_power = pockels_map_fine(upper_ind,2) - pockels_map_fine(lower_ind,2);
%         result_tmp.power{1}(these_trials) = pockels_map_fine(lower_ind,2) + interp_dist*step_size_power;
        these_trials_powerthresh = these_trials & result_tmp_current.max_curr{1} < 2500;
        result_tmp_current.peak_current_means(i) = mean(result_tmp_current.max_curr{1}(these_trials_powerthresh));

    end
%         
%     if ~isempty(filenames{j,2})
%         [nuclear_locs,fluor_vals,nuclear_locs_image_coord] = detect_nuclei(filenames{j,2},[],[],[],do_detect);
%         offsets = nuclear_locs_image_coord - cell_pos';
% 
%         [targ_error, index] = min(sqrt(sum(offsets.^2,1)));
%         result_tmp_current.fluor_val = fluor_vals(index);
%         result_tmp_spike.fluor_val = fluor_vals(index);
%         result_tmp_current.cell_pos = nuclear_locs_image_coord(:,index);
%         result_tmp_spike.cell_pos = nuclear_locs_image_coord(:,index);
%     else
%         result_tmp_current.fluor_val = NaN;
%         result_tmp_spike.fluor_val = NaN;
%         result_tmp_current.cell_pos = NaN;
%         result_tmp_spike.cell_pos = NaN;
%     end
%     subplot(211)

%     figure(current_figure)
%     gca
    subplot(222)
    hold on
    these_trials = result_tmp_current.max_curr{1} < 2500 & result_tmp_current.power{1} < 80;
    scatter(result_tmp_current.power{1}(these_trials),result_tmp_current.max_curr{1}(these_trials),15,colors(j,:),'jitter','on','jitteramount',.3);
    these_trials = unique_powers_current < 80;
    plot(unique_powers_current(these_trials),result_tmp_current.peak_current_means(these_trials),'-','color',colors(j,:),'Linewidth',1)
    result_current(j) = result_tmp_current;

%     figure(spike_figure)
%     gca
    subplot(221)
    hold on
    these_trials = result_tmp_spike.power{1} < 80;
    scatter(result_tmp_spike.power{1}.^1,result_tmp_spike.spike_times{1}/20,15,colors(j,:),'jitter','on','jitteramount',.3);
    these_trials = unique_powers_spike < 80;
    plot(unique_powers_spike,result_tmp_spike.spike_time_means/20,'-','color',colors(j,:),'Linewidth',1)
    result_spikes(j) = result_tmp_spike;

    
    
    
    
end

%         figure(current_figure)
%             gca
        subplot(222)
        title('Power vs. Peak Current')
        xlabel('Power (mW)')
        xlim([0 90])
        ylabel('Peak Current (pA)')
        ylim([0 2500])

%         figure(spike_figure)
%             gca
        subplot(221)
        title('Power vs. Spike Times')
        xlabel('Power (mW)')
        xlim([0 90])
        ylabel('Spike Times (msec)')
        ylim([0 15])


% legend({'cell 1','cell 2','cell 3'})

%%

for i = 1:size(filenames,1)
    load(filenames{j,1});
    [traces, ~,~, stim_traces] = get_traces(data,trials{i});
end
%%

% plot(26*ones(size(0:3000)),0:3000,'g--')
% plot(25.5*ones(size(0:3000)),0:3000,'b--')
plot(0:85,522.5*ones(size(0:85)),'b--')
plot(0:85,361.8*ones(size(0:85)),'g--')
plot(0:85,678.5*ones(size(0:85)),'r--')
% plot(14.9*ones(size(0:3000)),0:3000,'r--')

%%
% new_fig = 0;
% if ~exist('curr_vs_time','var')
%     new_fig = 1;
%     curr_vs_time = figure;
% else
%     figure(curr_vs_time)
% end
for i = 1:length(result_current)
    subplot(223)
    [shared_powers, current_i, spikes_i] = intersect(result_current(i).these_powers,result_spikes(i).these_powers);
    plot(result_current(i).peak_current_means(current_i),result_spikes(i).spike_time_means(spikes_i)/20,'o-','color',colors(i,:),'Linewidth',1)
    hold on
    title('peak current vs. spike time')
    xlabel('mean peak current (pA)')
    ylabel('spike time (msec)')
    xlim([0 2500])
    ylim([0 15])
    subplot(224)
    semilogy(result_current(i).peak_current_means(current_i),result_spikes(i).spike_time_jitter(spikes_i)/20,'o-','color',colors(i,:),'Linewidth',1)
    hold on
    title('peak current vs. spike time')
    xlabel('mean peak current (pA)')
    ylabel('spike time std. dev. (msec)')
    xlim([0 2500])
end


% ylim([0 15])
%%

% range = 8:13;
fluor_vals = [result_spikes.fluor_val];
transformed_flour_base = round(log([result_spikes.fluor_val]/max([result_spikes.fluor_val]))*100);
transformed_flour = transformed_flour_base - min(transformed_flour_base) + 1;
[fluor_order,sort_order] = sort(transformed_flour);
max_fluor = max(fluor_order);
colors = jet(28);
% transformed_flour = max_fluor - transformed_flour + 1;
% close all


% current_figure = figure;
% spike_figure = figure;
% transformed_fluor = [result_spikes.fluor_val];
% for i = [1:16 19:21 23:28]%setdiff(1:28,[1:16 19:21 23:28])
%     
% %     if result_current(i).inj_ratio == 10
%         figure(current_figure)
%         gca
%         hold on
%         scatter(result_current(i).power{1}.^1 * transformed_flour(i),result_current(i).max_curr{1},15,colors(min(round(transformed_flour(i)),max_fluor),:),'jitter','on','jitteramount',.3);
%         
% 
%         figure(spike_figure)
%         gca
%         hold on
%         scatter(result_spikes(i).power{1}.^1 * transformed_flour(i),result_spikes(i).spike_times{1}/20,15,colors(min(round(transformed_flour(i)),max_fluor),:),'jitter','on','jitteramount',.3);
%         ylabel('spike time (msec)')
%         xlabel('input power (mW)'); 
%         title('power vs. spike times, color is a function of log(fluor) - more gold is higher fluor')
% %     end
% end

figure
for ii = 1:length(orig_cells)%[1:16 19:21 23:27]%setdiff(1:28,[1:16 19:21 23:28])
    i = orig_cells(ii)
        subplot(121)
        hold on
        scatter(result_spikes(i).power{1}.^1,result_spikes(i).spike_times{1}/20,15,colors(i,:));%,colors(min(round(transformed_flour(i)),max_fluor),:));
        hold on
        all_powers = unique(round(result_spikes(i).power{1},1));
        plot(all_powers,result_spikes(i).spike_time_means/20,'color',colors(i,:))
        ylabel('spike time (msec)')
        xlabel('input est (a.u.)'); 
        title('power vs. spike times')
        
        subplot(122)
        hold on
        scatter(result_spikes(i).power{1}.^1 * gain_mle(ii),result_spikes(i).spike_times{1}/20,15,colors(i,:));%colors(min(round(transformed_flour(i)),max_fluor),:));
        hold on
        all_powers = unique(round(result_spikes(i).power{1},1));
        plot(all_powers * gain_mle(ii),result_spikes(i).spike_time_means/20,'color',colors(i,:))
        ylabel('spike time (msec)')
        xlabel('input intensity (power x gain, a.u.)'); 
        title('power x gain vs. spike times')
        
%         subplot(133)
%         hold on
%         scatter(result_spikes(i).power{1}.^1 * transformed_flour(i),result_spikes(i).spike_times{1}/20,15,colors(min(round(transformed_flour(i)),max_fluor),:),'jitter','on','jitteramount',.3);
%         ylabel('spike time (msec)')
%         xlabel('input est (a.u.)'); 
%         title('power*log(fluor)*100 vs. spike times')
end
%%

% range = 8:13;

[fluor_order,sort_order] = sort([result_spikes.fluor_val])
colors = copper(20);
% close all
current_figure = figure;
spike_figure = figure;
for i = 1:7
    min(round(result_current(i).vm_rest+80),20)
%     if result_current(i).inj_ratio == 10
        figure(current_figure)
        gca
        hold on
        scatter(result_current(i).power{1}.^1,result_current(i).max_curr{1},15,colors(min(round(result_current(i).vm_rest+80),20),:),'jitter','on','jitteramount',.3);


        figure(spike_figure)
        gca
        hold on
        scatter(result_spikes(i).power{1}.^1,result_spikes(i).spike_times{1}/20,15,colors(min(round(result_current(i).vm_rest+80),20),:),'jitter','on','jitteramount',.3);
%     end
end

%%

figure; 
% hold on
cell_select = [1:12]; 
plot([result_current_bu(cell_select).vm_rest],[result_current(cell_select).fluor_val],'or')
% xlim([-80 -50])

%% plot one cell trials

cell_id = 2;

[shared_powers, current_i, spikes_i] = intersect(result_current(cell_id).power{1},result_spikes(cell_id).power{1});

low_power = min(shared_powers); high_power = max(shared_powers);
med_power = shared_powers(floor(length(shared_powers)/2));
medlow_power = shared_powers(floor(length(shared_powers)/4));

colors = lines(4);


figure

time = (0:length(result_current(cell_id).traces{1}(1,:))-1)/20;

high_trial = find(result_current(cell_id).power{1} == high_power);
med_trial = find(result_current(cell_id).power{1} == med_power);
medlow_trial = find(result_current(cell_id).power{1} == medlow_power);
low_trial = find(result_current(cell_id).power{1} == low_power);

subplot(121)
hold on
plot(time,result_current(cell_id).traces{1}(low_trial,:)','linewidth',2,'color',colors(1,:))
plot(time,result_current(cell_id).traces{1}(med_trial,:)','linewidth',2,'color',colors(2,:))
plot(time,result_current(cell_id).traces{1}(medlow_trial,:)','linewidth',2,'color',colors(3,:))
plot(time,result_current(cell_id).traces{1}(high_trial,:)','linewidth',2,'color',colors(4,:))

title('Optically Evoked Currents')
xlabel('time (msec)')
ylabel('peak current (pA)')


high_trial = find(result_spikes(cell_id).power{1} == high_power);
med_trial = find(result_spikes(cell_id).power{1} == med_power);
medlow_trial = find(result_spikes(cell_id).power{1} == medlow_power);
low_trial = find(result_spikes(cell_id).power{1} == low_power);

subplot(122)
hold on
plot(time,result_spikes(cell_id).traces{1}(low_trial,:)','linewidth',2,'color',colors(1,:))
plot(time,result_spikes(cell_id).traces{1}(med_trial,:)','linewidth',2,'color',colors(2,:))
plot(time,result_spikes(cell_id).traces{1}(medlow_trial,:)','linewidth',2,'color',colors(3,:))
plot(time,result_spikes(cell_id).traces{1}(high_trial,:)','linewidth',2,'color',colors(4,:))

title('Optically Evoked Spikes')
xlabel('time (msec)')
ylabel('cell-attached recording (pA)')

legend({num2str(round(low_power)),num2str(round(medlow_power)),num2str(round(med_power)),num2str(round(high_power))})


%% plot one cell trials

cell_id = 4;
round_level = 3;
trial_count = 1;
cells_to_run = [1 2 3 4 5 6];
shared_powers = round(result_xy_bu(cell_id).these_x_power,round_level);

low_power = min(shared_powers); high_power = max(shared_powers);
med_power = shared_powers(floor(length(shared_powers)/2));
medlow_power = shared_powers(floor(length(shared_powers)/4));

colors = parula(length(cells_to_run)+2);
% colors = parula(ceil(length(shared_powers)/2)+1);
% colors = colors(randperm(length(shared_powers)),:);
% colors = [colors; colors(floor(length(shared_powers)/2):-1:1,:)];




figure

time = (0:length(result_xy_bu(cell_id).current_traces(1,:))-1)/20;

for ii = 1:length(cells_to_run)
    
    i = cells_to_run(ii);
    % high_trial = find(result_xy_bu(cell_id).curr_targ_power == high_power);
    % med_trial = find(result_xy_bu(cell_id).curr_targ_power == med_power);
    % medlow_trial = find(result_xy_bu(cell_id).curr_targ_power == medlow_power);
    % low_trial = find(result_xy_bu(cell_id).curr_targ_power == low_power);

    these_trials = find(round(result_xy_bu(cell_id).curr_targ_power,round_level) == shared_powers(i),trial_count);


    subplot(121)
    hold on
    % plot(time,result_xy_bu(cell_id).current_traces(low_trial,:)','linewidth',2,'color',colors(1,:))
    % plot(time,result_xy_bu(cell_id).current_traces(med_trial,:)','linewidth',2,'color',colors(2,:))
    % plot(time,result_xy_bu(cell_id).current_traces(medlow_trial,:)','linewidth',2,'color',colors(3,:))
    % plot(time,result_xy_bu(cell_id).current_traces(high_trial,:)','linewidth',2,'color',colors(4,:))

    plot(time,result_xy_bu(cell_id).current_traces(these_trials,:)','linewidth',2,'color',colors(ii,:))

    title('Optically Evoked Currents')
    xlabel('time (msec)')
    ylabel('peak current (pA)')


%     high_trial = find(result_xy_bu(cell_id).spike_targ_power == high_power);
%     med_trial = find(result_xy_bu(cell_id).spike_targ_power == med_power);
%     medlow_trial = find(result_xy_bu(cell_id).spike_targ_power == medlow_power);
%     low_trial = find(result_xy_bu(cell_id).spike_targ_power == low_power);

    these_trials = find(round(result_xy_bu(cell_id).spike_targ_power,round_level) == shared_powers(i),trial_count);

    subplot(122)
    hold on
%     plot(time,result_xy_bu(cell_id).spike_traces(low_trial,:)','linewidth',2,'color',colors(1,:))
%     plot(time,result_xy_bu(cell_id).spike_traces(med_trial,:)','linewidth',2,'color',colors(2,:))
%     plot(time,result_xy_bu(cell_id).spike_traces(medlow_trial,:)','linewidth',2,'color',colors(3,:))
%     plot(time,result_xy_bu(cell_id).spike_traces(high_trial,:)','linewidth',2,'color',colors(4,:))

    plot(time,result_xy_bu(cell_id).spike_traces(these_trials,:)','linewidth',2,'color',colors(ii,:))

    title('Optically Evoked Spikes')
    xlabel('time (msec)')
    ylabel('cell-attached recording (pA)')
end

% legend({num2str(round(low_power)),num2str(round(medlow_power)),num2str(round(med_power)),num2str(round(high_power))})




















