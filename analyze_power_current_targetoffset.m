%%

filenames_shape_offset = {'5_2_slice1_cell1next.mat','5_2_16_4_data.mat'};     
         

z_pos = [-10 0 10];         
spike_trials = {[]};
current_trials = {4:6};
% intrinsic_trials = {2,2,2,2,3,3,3,3,3,3,3};
% z_center = z_pos(2)*ones(size(filenames_shape_offset,1),1);

colors = jet(4);
% colors(:) = 0;
% colors(size(filenames_shape_offset,1)+1:size(filenames_shape_offset,1)*2) = 1; %green
% z_pos = [-8 0 8];


%%

clear result_shape_offset

for j = 1:size(filenames_shape_offset,1)
    
    load(filenames_shape_offset{j,2});
    load(filenames_shape_offset{j,1}); 
    
%     result_shape_offset(j).quadrant = experiment_setup.quadrant;
    
    if ~isempty(spike_trials{j})
        [result_shape_offset(j).spike_traces, ~, this_seq, result_shape_offset(j).spike_stim_traces, result_shape_offset(j).spike_full_stim_key] = get_traces(data,spike_trials{j});
        result_shape_offset(j).spike_targ_pos = result_shape_offset(j).spike_full_stim_key([this_seq.precomputed_target_index],:);
    
        result_shape_offset(j).spike_times = detect_peaks(-bsxfun(@minus,result_shape_offset(j).spike_traces,result_shape_offset(j).spike_traces(:,1)),20,30,0,Inf,-Inf,0,0,1);
        cell_spike_times = result_shape_offset(j).spike_times;
        result_shape_offset(j).spike_times = zeros(size(result_shape_offset(j).spike_times));

        for i = 1:length(cell_spike_times)
            if ~isempty(cell_spike_times{i})
                result_shape_offset(j).spike_times(i) = cell_spike_times{i};
            else
                result_shape_offset(j).spike_times(i) = NaN;
            end
        end
    end
    
    [result_shape_offset(j).current_traces, ~,this_seq, result_shape_offset(j).current_stim_traces, result_shape_offset(j).curr_full_stim_key] = get_traces(data,current_trials{j});
    result_shape_offset(j).current_targ_pos = ...
        [result_shape_offset(j).curr_full_stim_key([this_seq.precomputed_target_index],:)];
    result_shape_offset(j).max_curr = result_shape_offset(j).current_traces(:,1) - min(result_shape_offset(j).current_traces,[],2);

    
%     
%     tested_pos = unique(result_shape_offset(j).spike_targ_pos(:,1));
%     if ~isempty(setdiff(tested_pos,unique(result_shape_offset(j).current_targ_pos(:,1))))
%         disp([filenames_shape_offset{j,1} ' spike and current pos sets differ'])
%         continue
%     elseif ~isempty(setdiff(tested_pos,unique(result_shape_offset(j).spike_targ_pos(:,2))))
%         disp([filenames_shape_offset{j,1} ' y and x pos sets differ'])
%         continue
%     elseif ~isempty(setdiff(tested_pos,unique(result_shape_offset(j).current_targ_pos(:,2))))
%         disp([filenames_shape_offset{j,1} ' y and x pos sets differ across spike curr'])
%         continue
%     end
%     
%     these_x_power = x_measurements(:,result_shape_offset(j).quadrant);
%     these_y_power = y_measurements(:,result_shape_offset(j).quadrant);
%     result_shape_offset(j).spike_targ_power = zeros(size(result_shape_offset(j).spike_times));
%     result_shape_offset(j).curr_targ_power = zeros(size(result_shape_offset(j).max_curr));
%     for i = 1:length(tested_pos)
%         these_trials = result_shape_offset(j).spike_targ_pos(:,1) == 0 & result_shape_offset(j).spike_targ_pos(:,2) == tested_pos(i);
%         result_shape_offset(j).x_spike_time_means(i) = nanmean(result_shape_offset(j).spike_times(these_trials));
%         result_shape_offset(j).x_spike_time_jitter(i) = nanstd(result_shape_offset(j).spike_times(these_trials));
%         result_shape_offset(j).spike_targ_power(these_trials) = these_x_power(i);
%         these_trials = result_shape_offset(j).current_targ_pos(:,1) == 0 & result_shape_offset(j).current_targ_pos(:,2) == tested_pos(i);
%         result_shape_offset(j).x_max_curr_means(i) = nanmean(result_shape_offset(j).max_curr(these_trials));
%         result_shape_offset(j).curr_targ_power(these_trials) = these_x_power(i);
%         
%         these_trials = result_shape_offset(j).spike_targ_pos(:,2) == 0 & result_shape_offset(j).spike_targ_pos(:,1) == tested_pos(i);
%         result_shape_offset(j).y_spike_time_means(i) = nanmean(result_shape_offset(j).spike_times(these_trials));
%         result_shape_offset(j).y_spike_time_jitter(i) = nanstd(result_shape_offset(j).spike_times(these_trials));
%         result_shape_offset(j).spike_targ_power(these_trials) = these_y_power(i);
%         these_trials = result_shape_offset(j).current_targ_pos(:,2) == 0 & result_shape_offset(j).current_targ_pos(:,1) == tested_pos(i);
%         result_shape_offset(j).y_max_curr_means(i) = nanmean(result_shape_offset(j).max_curr(these_trials));
%         result_shape_offset(j).curr_targ_power(these_trials) = these_y_power(i);
%     end
%     
%     figure(pow_vs_cur_and_spike_time)
%     subplot(221)
% %     these_trials = result_shape_offset(j).spike_targ_pos(:,1) == 0;
%     scatter(result_shape_offset(j).spike_targ_power,result_shape_offset(j).spike_times/20,[],colors(experiment_setup.quadrant,:));
%     hold on
% %     plot(tested_pos,result_shape_offset(j).x_spike_time_means/20,'color',colors(experiment_setup.quadrant,:))
% %     xlim([-20 20])
% %     xlabel('Horizontal Distance (um)')
% %     ylabel('Spike Time (msec)')
% %     title('Horizontal Distance vs. Spike Time')
%     subplot(122)
% %     these_trials = result_shape_offset(j).current_targ_pos(:,1) == 0;
%     scatter(result_shape_offset(j).curr_targ_power,result_shape_offset(j).max_curr,[],colors(experiment_setup.quadrant,:));
% %     plot(tested_pos,result_shape_offset(j).x_max_curr_means,'color',colors(experiment_setup.quadrant,:))
%     hold on
% %     xlim([-20 20])
% %     xlabel('Horizontal Distance (um)')
% %     ylabel('Peak Current (pA)')
% %     title('Horizontal Distance vs. Peak Current')
%     
%     figure(pos_vs_cur_and_spike_time)
%     subplot(221)
%     these_trials = result_shape_offset(j).spike_targ_pos(:,1) == 0;
%     scatter(result_shape_offset(j).spike_targ_pos(these_trials,2),result_shape_offset(j).spike_times(these_trials)/20,[],colors(experiment_setup.quadrant,:));
%     hold on
%     plot(tested_pos,result_shape_offset(j).x_spike_time_means/20,'color',colors(experiment_setup.quadrant,:))
%     xlim([-20 20])
%     xlabel('Horizontal Distance (um)')
%     ylabel('Spike Time (msec)')
%     title('Horizontal Distance vs. Spike Time')
%     subplot(222)
%     these_trials = result_shape_offset(j).current_targ_pos(:,1) == 0;
%     scatter(result_shape_offset(j).current_targ_pos(these_trials,2),result_shape_offset(j).max_curr(these_trials),[],colors(experiment_setup.quadrant,:));
%     plot(tested_pos,result_shape_offset(j).x_max_curr_means,'color',colors(experiment_setup.quadrant,:))
%     hold on
%     xlim([-20 20])
%     xlabel('Horizontal Distance (um)')
%     ylabel('Peak Current (pA)')
%     title('Horizontal Distance vs. Peak Current')
%     subplot(223)
%     these_trials = result_shape_offset(j).spike_targ_pos(:,2) == 0;
%     scatter(result_shape_offset(j).spike_targ_pos(these_trials,1),result_shape_offset(j).spike_times(these_trials)/20,[],colors(experiment_setup.quadrant,:));
%     hold on
%     plot(tested_pos,result_shape_offset(j).y_spike_time_means/20,'color',colors(experiment_setup.quadrant,:))
%     xlim([-20 20])
%     xlabel('Vertical Distance (um)')
%     ylabel('Spike Time (msec)')
%     title('Vertical Distance vs. Spike Time')
%     subplot(224)
%     these_trials = result_shape_offset(j).current_targ_pos(:,2) == 0;
%     scatter(result_shape_offset(j).current_targ_pos(these_trials,1),result_shape_offset(j).max_curr(these_trials),[],colors(experiment_setup.quadrant,:));
%     hold on
%     plot(tested_pos,result_shape_offset(j).y_max_curr_means,'color',colors(experiment_setup.quadrant,:))
%     xlim([-20 20])
%     xlabel('Vertical Distance (um)')
%     ylabel('Peak Current (pA)')
%     title('Vertical Distance vs. Peak Current')
%     
%     figure(curr_vs_time)
%     subplot(121)
%     [ordered_curr, curr_order_x] = sort(result_shape_offset(j).x_max_curr_means);
%     plot(result_shape_offset(j).x_max_curr_means(curr_order_x),result_shape_offset(j).x_spike_time_means(curr_order_x)/20,'-r','Linewidth',2);
%     hold on
%     [ordered_curr, curr_order_y] = sort(result_shape_offset(j).y_max_curr_means);
%     plot(result_shape_offset(j).y_max_curr_means(curr_order_y),result_shape_offset(j).y_spike_time_means(curr_order_y)/20,'-r','Linewidth',2);
%     ylim([0 15])
%     if new_fig
%         xlabel('Peak Current Mean (pA)')
%         ylabel('Spike Time Mean (msec)')
%         title('Peak Current vs. Spike Time')
%     end
%     subplot(122)
%     semilogy(result_shape_offset(j).x_max_curr_means(curr_order_x),result_shape_offset(j).x_spike_time_jitter(curr_order_x)/20,'-r','Linewidth',2);
%     hold on
%     semilogy(result_shape_offset(j).y_max_curr_means(curr_order_y),result_shape_offset(j).y_spike_time_jitter(curr_order_y)/20,'-r','Linewidth',2);
%     if new_fig
%         xlabel('Peak Current Mean (pA)')
%         ylabel('Spike Time Std. Dev. (msec)')
%         title('Peak Current vs. Spike Time Jitter')
%     end

    
    
    
    
end

%%
do_detect = 0;

for j = 1:size(filenames_shape_offset,1)
    if ~isempty(filenames_shape_offset{j,3})
        
        disp('...')
        load(filenames_shape_offset{j,2});
        [nuclear_locs,fluor_vals,nuclear_locs_image_coord] = detect_nuclei(filenames_shape_offset{j,3},[],[],[],do_detect,[],0);
        offsets = nuclear_locs - [experiment_setup.center_pos_um(1:2) 30];

        [targ_error, index] = min(sqrt(sum(offsets.^2,2)));
        result_shape_offset(j).fluor_val = fluor_vals(index);
        result_shape_offset(j).cell_pos = nuclear_locs(index,:);
        result_shape_offset(j).exp_cell_pos = [experiment_setup.center_pos_um(1:2) 30];
        result_shape_offset(j).err_cell_pos = result_shape_offset(j).exp_cell_pos - result_shape_offset(j).cell_pos;

        
        figure
        plot_nuclear_detect_3D([filenames_shape_offset{j,3} '.tif'],nuclear_locs_image_coord);
        hold on
        scatter(result_shape_offset(j).cell_pos(2)/1.89 + 118.2,result_shape_offset(j).cell_pos(1)/1.89 + 124.5,'go')
        scatter(result_shape_offset(j).exp_cell_pos(2)/1.89 + 118.2,result_shape_offset(j).exp_cell_pos(1)/1.89 + 124.5,'bx')
        max_curr_tmp = result_shape_offset(j).max_curr;
        max_curr_tmp(max_curr_tmp > 2000) = 0;
        [max_curr_spatial, max_curr_ind] = max(max_curr_tmp);
        max_curr_loc = result_shape_offset(j).current_targ_pos(max_curr_ind,1:2) + experiment_setup.center_pos_um(1:2);
        scatter(max_curr_loc(2)/1.89 + 118.2,max_curr_loc(1)/1.89 + 124.5,'mo')
        title(['Cell, ' num2str(j) ', Error (um): ' mat2str(result_shape_offset(j).err_cell_pos)])
        caxis([0 400])
             
    else
        
        result_shape_offset(j).fluor_val = NaN;
        result_shape_offset(j).cell_pos = [experiment_setup.center_pos_um(1:2) experiment_setup.piezo_center];

    end
end
%%

z_pos = [30 40 50];%[160 200 240];
figure
count = 1;

% [ha,~] = tight_subplot(length(z_pos),size(filenames_shape_offset,1),.01,.2,.2);
for i = 1:length(z_pos)
    

    for j = 1:length(result_shape_offset)
        count
        subplot(length(z_pos),length(result_shape_offset),count)
%         axes(ha(count))
        scatter3(result_shape_offset(j).current_targ_pos(round(result_shape_offset(j).current_targ_pos(:,3)) == z_pos(i) & result_shape_offset(j).max_curr < 2000,1), ...
                 result_shape_offset(j).current_targ_pos(round(result_shape_offset(j).current_targ_pos(:,3)) == z_pos(i) & result_shape_offset(j).max_curr < 2000,2), ...
                 result_shape_offset(j).current_targ_pos(round(result_shape_offset(j).current_targ_pos(:,3)) == z_pos(i) & result_shape_offset(j).max_curr < 2000,3),[], ...
                 result_shape_offset(j).max_curr(round(result_shape_offset(j).current_targ_pos(:,3)) == z_pos(i) & result_shape_offset(j).max_curr < 2000)...%/...
                 )%max(result_shape_offset(j).max_curr(result_shape_offset(j).max_curr < 2000)))
%              zlim([0 2000])
             caxis([0 1400])
             view([90 90])
             
         count = count + 1;

    end
end

%%

z_pos = [160 200 240];

count = 1;
[xq, yq] = meshgrid(-35:5:35,-35:5:35);
spatial_maps = zeros(length(xq),length(yq),length(z_pos),length(result_shape_offset));
spatial_maps_per_cell = struct();

all_curr = [];
all_pos = [];

    

for j = 1:length(result_shape_offset)
        
    for i = 1:length(z_pos)
        these_xy = unique(result_shape_offset(j).current_targ_pos(result_shape_offset(j).current_targ_pos(:,3) == z_pos(i) & result_shape_offset(j).max_curr < 2000,1:2),'rows');
        curr_means = zeros(size(these_xy,1),1);
        x = [];
        y = [];
        z = [];
        curr = [];
        
        for k = 1:size(these_xy,1)
            curr_mean = nanmean(result_shape_offset(j).max_curr(result_shape_offset(j).current_targ_pos(:,3) == z_pos(i) & result_shape_offset(j).max_curr < 2000 & ...
                result_shape_offset(j).current_targ_pos(:,1) == these_xy(k,1) & result_shape_offset(j).current_targ_pos(:,2) == these_xy(k,2)));
            if ~isnan(curr_mean)
                x = [x these_xy(k,1)];
                y = [y these_xy(k,2)];
                z = [z z_pos(i)];
                curr = [curr curr_mean];
            end
        end
        if z_pos(i) == 200
            [max_curr, max_pos] = max(curr);
            max_pos = [x(max_pos) y(max_pos) z(max_pos)];
        end
       spatial_maps(:,:,i,j) = griddata(x,y,curr,xq,yq);
    end
    spatial_maps_per_cell(j).curr = curr;
    spatial_maps_per_cell(j).pos = [x; y; z]';% - max_pos;
    all_curr = [all_curr spatial_maps_per_cell(j).curr];
    all_pos = [all_pos; spatial_maps_per_cell(j).pos];
end

%%

[all_unique_pos, all_unique_pos_inds, all_pos_inds] = unique(all_pos,'rows');
mean_curr_all_cells = zeros(size(all_unique_pos,1),1);
var_curr_all_cells = zeros(size(all_unique_pos,1),1);

for i = 1:length(mean_curr_all_cells)
    
    mean_curr_all_cells(i) = nanmean(all_curr(all_pos_inds == i));
    var_curr_all_cells(i) = nanvar(all_curr(all_pos_inds == i));
    
end

mean_curr_all_cells = mean_curr_all_cells/max(mean_curr_all_cells);

mean_curr_map = zeros(length(xq),length(yq),length(z_pos));
var_curr_map = zeros(length(xq),length(yq),length(z_pos));

z_pos = [160 200 240];
for i = 1:length(z_pos)
    these_z_pos = all_unique_pos(:,3) == z_pos(i);
    these_xy = all_unique_pos(these_z_pos,1:2);
    mean_curr_map(:,:,i) = griddata(these_xy(:,1),these_xy(:,2),mean_curr_all_cells(these_z_pos),xq,yq);
    var_curr_map(:,:,i) = griddata(these_xy(:,1),these_xy(:,2),var_curr_all_cells(these_z_pos),xq,yq);
end

figure;
for i = 1:size(mean_cell,3)
    subplot(3,2,(i-1)*2+1)
    imagesc(mean_curr_map(:,:,i)')
    caxis([0 max(mean_curr_map(:))])
    title(['Mean Cell Shape (Z = ' num2str(z_pos(i) - 200) ' um)'])
    
    subplot(3,2,(i-1)*2+2)
    imagesc(var_curr_map(:,:,i)')
    caxis([0 max(var_curr_map(:))])
    title(['Variance Cell Shape (Z = ' num2str(z_pos(i) - 200) ' um)'])
end
%%

figure
count = 1;
for i = [1 2 3]%1:size(spatial_maps,3)
    for j = [8 10 7 9 6]%1:length(result_shape_offset)
         
        subplot(3,5,count)
        b = bar3(spatial_maps(:,:,i,j)')
        
        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
        end
        
        caxis([0 max(max(max(spatial_maps(:,:,:,j))))])
        zlim([0 1800])
        
%         if count > 5
            view(0,90)
%         end
        axis off 
        count = count + 1;
    end
end

%%

% for i = 1:size(spatial_maps,4)
%     
%     spatial_maps_norm(:,:,:,i) = spatial_maps(:,:,:,i)/max(max(max(spatial_maps(:,:,:,i))));
% end



mean_cell = mean(spatial_maps,4);
std_cell = var(spatial_maps,[],4);
figure;
for i = 1:size(mean_cell,3)
    subplot(3,2,(i-1)*2+1)
    imagesc(mean_cell(:,:,i)')
    caxis([0 max(mean_cell(:))])
    title(['Mean Cell Shape (Z = ' num2str(z_pos(i) - 200) ' um)'])
    
    subplot(3,2,(i-1)*2+2)
    imagesc(std_cell(:,:,i)')
    caxis([0 max(std_cell(:))])
    title(['Variance Cell Shape (Z = ' num2str(z_pos(i) - 200) ' um)'])
end
    
