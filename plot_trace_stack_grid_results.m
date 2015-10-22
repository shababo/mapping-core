function plot_trace_stack_grid_results(traces_array,results)

% figure
[num_rows, num_cols] = size(traces_array);

trace_count = 1;
plot_count = 1;
for i = 1:num_rows
    for j = 1:num_cols
        num_traces = size(traces_array{i,j},2);
          subplot(num_rows,num_cols,plot_count)
%  figure
%         these_traces = traces_array{i,j}';
%         plot_trace_stack(these_traces,zeros(size(these_traces)),[],zeros(length(these_traces),3),[],size(these_traces,2)-1,125)
        plot_detection_results(results(trace_count:trace_count+num_traces-1))
        trace_count = trace_count + num_traces;
        plot_count = plot_count + 1;
    end
    
end