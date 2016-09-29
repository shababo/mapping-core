figure;

offset = 100;

t = (0:size(traces,2)-1)/20000;

plot(t,mean(traces(33:38,:)),'r')
hold on;
plot(t,mean(traces(39:43,:))-offset,'g')
hold on;
plot(t,mean(traces(44:49,:)) - 2*offset,'c')
hold on;
plot(t,mean(traces(92:100,:)) - 3*offset,'b')
% hold on
% plot(t,mean(traces(142:146,:)) - 4*offset,'k')
hold off

legend({'950','900','850','800'})
title('Chronos Spectrum (one cell)')
xlabel('Time (sec)')
ylabel('Current (pA)')

%%

hack(3,1,'2','4_5','min','','pulseamp','100')
hack(3,1,'2','4_5','min','','pulseamp','50')
hack(3,1,'2','4_5','min','','pulseamp','25')
hack(4,1,'2','4_5','min','','pulseamp','25')
hack(4,1,'2','4_5','min','','pulseamp','50')
hack(4,1,'2','4_5','min','','pulseamp','100')

%%

figure
subplot(131)
plot_trace_stack_grid(traces_by_location_4_5_s3c1_r3_50mw,3,1,0);

num_traces = 10;
currents_diff_locations = zeros(num_traces,size(traces_by_location_4_5_s3c1_r3_50mw{1,1},2));
location_ids = cell(num_traces,1);

currents_diff_locations(1,:) = traces_by_location_4_5_s3c1_r3_50mw{6,1}(1,:);
location_ids{1} = '6,1';
currents_diff_locations(2,:) = traces_by_location_4_5_s3c1_r3_50mw{1,6}(1,:);
location_ids{2} = '1,6';
currents_diff_locations(3,:) = traces_by_location_4_5_s3c1_r3_50mw{5,6}(1,:);
location_ids{3} = '5,6';
currents_diff_locations(4,:) = traces_by_location_4_5_s3c1_r3_50mw{8,4}(1,:);
location_ids{4} = '8,4';
currents_diff_locations(5,:) = traces_by_location_4_5_s3c1_r3_50mw{6,6}(1,:);
location_ids{5} = '6,6';
currents_diff_locations(6,:) = traces_by_location_4_5_s3c1_r3_50mw{11,5}(1,:);
location_ids{6} = '11,5';
currents_diff_locations(7,:) = traces_by_location_4_5_s3c1_r3_50mw{10,6}(1,:);
location_ids{7} = '10,6';
currents_diff_locations(8,:) = traces_by_location_4_5_s3c1_r3_50mw{4,7}(1,:);
location_ids{8} = '4,7';
currents_diff_locations(9,:) = traces_by_location_4_5_s3c1_r3_50mw{9,11}(1,:);
location_ids{9} = '9,11';
currents_diff_locations(10,:) = traces_by_location_4_5_s3c1_r3_50mw{6,3}(1,:);
location_ids{10} = '6,3';

[~,order] = sort(min(currents_diff_locations,[],2),'descend');
subplot(132); plot_trace_stack(currents_diff_locations(order,:),150,zeros(size(currents_diff_locations,1),3),'-');
legend(location_ids(order))

norm_currents_diff_locations = bsxfun(@minus,currents_diff_locations,mean(currents_diff_locations(:,1:99),2));
norm_currents_diff_locations = bsxfun(@rdivide,norm_currents_diff_locations,min(norm_currents_diff_locations,[],2)+5);
colors = hsv(10);

subplot(133); 
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(-norm_currents_diff_locations(order,:)');
legend(location_ids(order))
axis off

%%
set(gcf, 'Color', 'w');
    
export_fig(gcf, 'figures/st-chr2-spatial-shape-comparison.pdf','-pdf')

%%


