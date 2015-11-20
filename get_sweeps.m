function [traces, traces_metadata] = get_sweeps(filename,ch_ind,trace_inds,plot_stack,varargin)

load(filename,'sweeps','data')

if length(varargin) > 0
    max_sweep = varargin{1};
else
    max_sweep = Inf;
end

p = inputParser;
addParameter(p,'amp_units','ignore',@(x) ischar(x) || cellstr(x));
addParameter(p,'pulseamp',NaN,@isvector);
addParameter(p,'pulseduration',NaN,@isvector);
addParameter(p,'pulsefrequency',NaN,@isvector);
addParameter(p,'pulse_starttime',NaN,@isvector);
addParameter(p,'stim_type','ignore',@(x) ischar(x) || cellstr(x));
addParameter(p,'hologram_id','ignore',@(x) ischar(x) || cellstr(x));
addParameter(p,'note','ignore',@(x) ischar(x) || cellstr(x));
addParameter(p,'lut_used',NaN,@isnumeric);
addParameter(p,'run_count',NaN,@isnumeric);
addParameter(p,'clamp_type','ignore',@(x) ischar(x) || cellstr(x));

parse(p,varargin{2:end})

assignin('base','p',p)

if ~exist('data','var')
    data.sweeps = sweeps;
end

if ~isempty(trace_inds)
    trace_inds = trace_inds(trace_inds <= max_sweep);
    trace_array = data.sweeps(trace_inds);
    traces_metadata = data.trial_metadata(trace_inds);
else
    if length(data.sweeps) > max_sweep
       trace_array = data.sweeps(1:max_sweep);
        traces_metadata = data.trial_metadata(1:max_sweep);
    else
        trace_array = data.sweeps;
        traces_metadata = data.trial_metadata;
    end
end

start_ind = 1;
end_ind = length(trace_array{1});

traces = zeros(length(trace_array),length(trace_array{1}(start_ind:end_ind,ch_ind)'));

if ~exist('sweeps','var')
    match = zeros(length(trace_array),1);
end

for i = 1:length(trace_array)
    
    if ~exist('sweeps','var')
        match(i) = match_trial(p.Results, data.trial_metadata(i));
    end
    traces(i,:) = trace_array{i}(start_ind:end_ind,ch_ind)'; 
    
end


if ~exist('sweeps','var')
    match = logical(match);
    traces = traces(match,:);
    traces_metadata = traces_metadata(match);
end

% traces = bsxfun(@minus,traces,mean(traces,1));

if plot_stack
    figure; plot_trace_stack(traces,100,zeros(length(traces),3),'-')
end

