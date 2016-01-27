function trace_array = build_trace_grid(traces, metadata)

[traces, traces_metadata] = get_sweeps_dir('data','12_3_slice6_cell2',0,1,0,Inf,'run_count',8);
get_sweeps_dir(dirname,match_string,recursive,ch_ind,plot,max_sweep,varargin)
[traces, traces_metadata] = get_sweeps_dir('data','12_3_slice6_cell2',0,1,0,Inf,'run_count',8);
p.Results
[traces, traces_metadata] = get_sweeps_dir('data','12_3_slice6_cell2',0,1,0,Inf,'run_count',8);
traces = traces{1};
traces_metadata = traces_metadata{1};
traces_by_location = cell(11,11);
start_ind = 20000*.280; end_ind = 20000*.400;
for i = 1:size(traces,1)
ind1 = traces_metadata(i).relative_position(1)/10 + 6;
ind2 = traces_metadata(i).relative_position(2)/10 + 6;
traces_by_location{ind1,ind2} = [traces_by_location{ind1,ind2} traces(i,start_ind:end_ind)'];
end
figure; plot_trace_stack_grid(traces_by_location,5,10,0)
for i = 1:size(traces,1)
ind1 = traces_metadata(i).relative_position(1)/10 + 6;
ind2 = traces_metadata(i).relative_position(2)/10 + 6;
traces_by_location{ind1,ind2} = [traces_by_location{ind1,ind2} traces(i,start_ind:end_ind)'];
end
for i = 1:size(traces,1)
ind1 = traces_metadata(i).relative_position(1)/10 + 6;
ind2 = traces_metadata(i).relative_position(2)/10 + 6;
traces_by_location{ind1,ind2} = [traces_by_location{ind1,ind2}; traces(i,start_ind:end_ind)];
end
clear traces_by_location
for i = 1:size(traces,1)
ind1 = traces_metadata(i).relative_position(1)/10 + 6;
ind2 = traces_metadata(i).relative_position(2)/10 + 6;
traces_by_location{ind1,ind2} = [traces_by_location{ind1,ind2}; traces(i,start_ind:end_ind)];
end
traces_by_location = cell(11,11);
for i = 1:size(traces,1)
ind1 = traces_metadata(i).relative_position(1)/10 + 6;
ind2 = traces_metadata(i).relative_position(2)/10 + 6;
traces_by_location{ind1,ind2} = [traces_by_location{ind1,ind2}; traces(i,start_ind:end_ind)];
end
figure; plot_trace_stack_grid(traces_by_location,5,10,0)
evoked_traces = traces_by_location{9,4};
figure; plot(mean(evoked_traces))
evoked_traces(end+1,:) = mean(evoked_traces);
figure; plot_trace_stack(evoked_traces,25,zeros(size(evoked_traces,1),3),'-')
figure; plot_trace_stack(evoked_traces,75,zeros(size(evoked_traces,1),3),'-')