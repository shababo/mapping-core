load('data/5_13_s2c1_r4_tracegrid-1000.mat')
load('data/5_13_s2c1_r4_tracegrid.mat')
results_5_12_s2c1_r4_tracegrid_1000 = results;
traces_5_12_s2c1_r4_tracegrid_1000 = traces;


time_posteriors = zeros(length(results_5_12_s2c1_r4_tracegrid_1000),2000);


for i = 1:length(results_5_12_s2c1_r4_tracegrid_1000)
    time_posteriors(i,:) = histcounts(results_5_12_s2c1_r4_tracegrid_1000(i).trials.times,0:2000);
end

time_posteriors(:,1:50) = 0;
posteriors_grid = unstack_traces(time_posteriors/5,params.rebuild_map);

figure; compare_trace_stack_grid_overlap({posteriors_grid,traces_5_12_s2c1_r4_tracegrid_1000},3,1,[],0,{'L4','L5'},1)
