filenames = {'2_9_slice2_cell2.mat', '2_9_17_35_data.mat'
             '2_9_slice2_cell3.mat', '2_9_17_54_data.mat'};     
         
trials = {3:8,3:8,3:8};
         
z_pos = [-8 0 8];
%%

images_fig = figure

for j = 1:2
    
    clear result(j)
    load(filenames{j,2});
    load(filenames{j,1}); 
    
    [result(j).traces, ~,~, result(j).stim_traces] = get_traces(data,trials{j});
    result(j).max_curr = result(j).traces(:,1) - min(result(j).traces,[],2);
    
    result(j).curr_map = zeros(11,11,3);
    xy_id_order = repmat(experiment_setup.xy_id_order,1,2);
    z_id_order = repmat(experiment_setup.z_id_order,1,2);
    for i = 1:size(result(j).max_curr)
        [x1 x2] = ind2sub([11 11],xy_id_order(i));
        result(j).curr_map(x1,x2,z_id_order(i)) = result(j).max_curr(i);
    end
    
    figure(images_fig)
    for i = 1:3
        subplot(size(filenames,1),3,(j-1)*3 + i)
        [X,Y] = meshgrid(unique(experiment_setup.all_xy_targets(:,1)),...
            unique(experiment_setup.all_xy_targets(:,2)));
        surf(X,Y,result(j).curr_map(:,:,i))
        caxis([0 max(result(j).max_curr)])
        xlabel('x dist (um)'); ylabel('y dist (um)'), zlabel('current (pA)');
        title(['Cell ' num2str(j) ': z = ' num2str(z_pos(i))])
    end
    
%     all_cells_result(j) = result(j);
    
end

%%
z_pos = [-8 0 8];
figure;
colors = [0 1 0; 0 0 1];
for j = 1:2
    d
    
    max_curr = max(result(j).max_curr);
    for i = 1:size(result(j).max_curr)
        scatter3(experiment_setup.all_xy_targets(xy_id_order(i),2),...
            experiment_setup.all_xy_targets(xy_id_order(i),1),...
        	z_pos(experiment_setup.z_id_order(i)),...
            round(20*result(j).max_curr(i)/max_curr)+1,...
            colors(j,:));
            hold on;
    end
    
end
    
    
%%
    
    
    
    
    