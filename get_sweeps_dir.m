function [traces, traces_metadata] = get_sweeps_dir(dirname,match_string,recursive,ch_ind,plot,max_sweep,varargin)

if iscellstr(match_string)
    match_string = cellstr2regexp(match_string);
end
if recursive
    dirinfo = rdir(dirname);
else
    dirinfo = dir(dirname);
end

traces = {};
traces_metadata = {};
num_traces = 0;
count = 1;
for i = 1:length(dirinfo)
    
    if ~dirinfo(i).isdir && ~isempty(regexpi(dirinfo(i).name,match_string))
        [traces{count}, traces_metadata{count}] = ...
            get_sweeps([dirname '/' dirinfo(i).name],ch_ind,[],0,max_sweep,varargin{:});
        
        num_traces = num_traces + size(traces{count},1);
        count = count + 1;
    end
end

if plot

    traces_plot = zeros(num_traces,size(traces{1},2));
    colors = zeros(num_traces,3);
    first = 1;
    last = 0;
    count = 1;
    for i = 1:length(traces)
        if ~isempty(traces{i})
            last = last + size(traces{i},1);
            traces_plot(first:last,:) = traces{i};
            colors(first:last,:) = repmat([0 0 mod(count,2)],last-first+1,1);
            first = first + size(traces{i},1);
            count = count + 1;
        end
    end

    figure
    plot_trace_stack(traces_plot,100,colors,'-');
end

    