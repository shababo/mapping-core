load('5_13_s2c2_r4_tracegrid-2000.mat')
load('5_13_s2c2_r4_tracegrid.mat')
results_5_13_s2c2_r4_tracegrid_2000 = results;
traces_5_13_s2c2_r4_tracegrid = traces;


time_posteriors = zeros(length(results_5_13_s2c2_r4_tracegrid_2000),2000);


for i = 1:length(results_5_13_s2c2_r4_tracegrid_2000)
    time_posteriors(i,:) = histcounts(results_5_13_s2c2_r4_tracegrid_2000(i).trials.times,0:2000);
end

time_posteriors(:,1:50) = 0;
posteriors_grid_5_13_s2c2_r4_2000 = unstack_traces(time_posteriors,params.rebuild_map);

figure; compare_trace_stack_grid_overlap({posteriors_grid_5_13_s2c2_r4_2000,traces_5_13_s2c2_r4_tracegrid},3,1,[],0,{'L4','L5'},1)

%%

som_posterios = {posteriors_grid_5_13_s2c1_r4_2000,posteriors_grid_5_13_s2c2_r4_2000};
pv_posterios = {posteriors_grid_5_12_s2c1_r4_2000,posteriors_grid_5_12_s2c2_r4_2000};

%%

detection_grid_5_12_s2c2_r4_2000 = cell(size(posteriors_grid_5_12_s2c2_r4_2000));



all_points = [];
for i = 1:size(posteriors_grid_5_12_s2c2_r4_2000,1)
    for j = 1:size(posteriors_grid_5_12_s2c2_r4_2000,2)
        for k = 1:size(posteriors_grid_5_12_s2c2_r4_2000{i,j},1)
            all_points = [all_points posteriors_grid_5_12_s2c2_r4_2000{i,j}(k,:)];
        end
    end
end

std_all_points = std(all_points);
delete all_points

for i = 1:size(posteriors_grid_5_12_s2c2_r4_2000,1)
    for j = 1:size(posteriors_grid_5_12_s2c2_r4_2000,2)
        
        detection_grid_5_12_s2c2_r4_2000{i,j} = zeros(size(posteriors_grid_5_12_s2c2_r4_2000{i,j}));
        
        for k = 1:size(posteriors_grid_5_12_s2c2_r4_2000{i,j},1)
            [~, event_times] = findpeaks(posteriors_grid_5_12_s2c2_r4_2000{i,j}(k,:),'MinPeakHeight',1.0*std_all_points,'MinPeakDistance',60);
            detection_grid_5_12_s2c2_r4_2000{i,j}(k,event_times) = 50;
        end
    end
end

figure; compare_trace_stack_grid_overlap({detection_grid_5_12_s2c2_r4_2000,traces_5_12_s2c2_r4_tracegrid},3,1,[],0,{'L4','L5'},1)
















%%

event_timeseries2 = get_event_times_init(results_5_12_s2c2_r4_tracegrid_2000,2000,1,[]);
event_timeseries1 = get_event_times_init(results_5_12_s2c1_r4_tracegrid_2000,2000,1,[]);
event_timeseries1_smooth = smoothts(event_timeseries1,'g',100,20);
event_timeseries2_smooth = smoothts(event_timeseries2,'g',100,20);
events_ts_grid1_smooth = unstack_traces(event_timeseries1_smooth*500,params.rebuild_map);
events_ts_grid2_smooth = unstack_traces(event_timeseries2_smooth*500,params.rebuild_map);
% figure; compare_trace_stack_grid([traces_by_location_5_12_s2c1_2_r4(:); {events_ts_grid1_smooth}; {events_ts_grid2_smooth}],3,1,[],0,{'L4','L5'},1)
figure; compare_trace_stack_grid_overlap({events_ts_grid1_smooth,traces_5_12_s2c1_r4_tracegrid},3,1,[],0,{'L4','L5'},1)

%%
all_points = [];
for i = 1:length(results)
    all_points = [all_points results(i).filtered_trace];
end

std_all_points = std(all_points);
delete all_points


results_new = results;
min_window = 20;
threshold = 2.25;


for i = 1:length(results_new)
    [~, results_new(i).event_times_init] = findpeaks(results_new(i).filtered_trace,'MinPeakHeight',threshold*std_all_points,'MinPeakDistance',min_window);
%     results_new(i).event_sizes_init = zeros(size(results_new(i).event_times_init));
%     for j = 1:length(results_new(i).event_times_init)
%         baseline = min(trace(max(1,results_new(i).event_times_init(j)-40):results_new(i).event_times_init(j)));
%         results_new(i).event_sizes_init(j) = max(trace(results_new(i).event_times_init(j):min(results_new(i).event_times_init(j)+200,length(trace)))) - 10 - baseline;
%     end

    results_new(i).event_times_init(results_new(i).event_times_init < 5) = [];
%     results_new(i).event_sizes_init(results_new(i).event_sizes_init < 5) = [];
end

