function [traces, traces_metadata] = get_sweeps(filename,ch_ind,trace_inds,plot_stack,varargin)

load(filename,'sweeps','data')


if ~isempty(varargin)
    max_sweep = varargin{1};
else
    max_sweep = Inf;
end

p = inputParser;
addParameter(p,'amp_units','ignore',@(x) ischar(x) || iscellstr(x));
addParameter(p,'pulseamp',NaN,@isvector);
addParameter(p,'pulseduration',NaN,@isvector);
addParameter(p,'pulsefrequency',NaN,@isvector);
addParameter(p,'pulse_starttime',NaN,@isvector);
addParameter(p,'stim_type','ignore',@(x) ischar(x) || iscellstr(x));
addParameter(p,'hologram_id','ignore',@(x) ischar(x) || iscellstr(x));
addParameter(p,'note','ignore',@(x) ischar(x) || iscellstr(x));
addParameter(p,'lut_used',NaN,@isnumeric);
addParameter(p,'run_count',NaN,@isvector);
addParameter(p,'clamp_type','ignore',@(x) ischar(x) || iscellstr(x));
addParameter(p,'relative_to_start_position',[NaN NaN NaN],@isnumeric);


% varargin{:}
p.parse(varargin{2:end})

assignin('base','p',p)

if ~exist('data','var')
    data.sweeps = sweeps;
end

if ~isempty(trace_inds)
    trace_inds = trace_inds(trace_inds <= max_sweep);
    trace_array = data.sweeps(trace_inds);
    if isfield(data,'trial_metadata')
        traces_metadata = data.trial_metadata(trace_inds);
    else
        traces_metadata = [];
    end
else
    if length(data.sweeps) > max_sweep
       trace_array = data.sweeps(1:max_sweep);
       if isfield(data,'trial_metadata')
           traces_metadata = data.trial_metadata(1:max_sweep);
       else
           traces_metadata = [];
       end
    else
        trace_array = data.sweeps;
        if isfield(data,'trial_metadata')
            traces_metadata = data.trial_metadata;
        else
            traces_metadata = [];
        end
    end
end

% start_ind = 1;
% end_ind = Inf;


% traces = zeros(length(trace_array),1);


if ~exist('sweeps','var') && isfield(data,'trial_metadata')
    match = zeros(length(trace_array),1);
end

count = 1;

for i = 1:length(trace_array)
    
    if ~exist('sweeps','var')  && isfield(data,'trial_metadata')
        match(i) = match_trial(p.Results, data.trial_metadata(i)) && ...
            ((~isfield(data.trial_metadata,'test_trial') || ~data.trial_metadata(i).test_trial) && (~isfield(data.trial_metadata,'is_good_trial') || data.trial_metadata(i).is_good_trial));
%         if ~data.trial_metadata(i).is_good_trial
%             disp(match(i))
%         end
        if match(i) 
            
            if (exist('traces','var') && length(trace_array{i}(:,ch_ind)) == size(traces,2)) || ~exist('traces','var')
                traces(count,:) = trace_array{i}(:,ch_ind)';
                count = count + 1;
            end
        end
    
    else
        traces(i,:) = trace_array{i}(:,ch_ind)';
    end

     
    
end



if ~exist('sweeps','var')  && isfield(data,'trial_metadata')
    match = logical(match);
%     traces = traces(match,:);
    if isfield(data,'trial_metadata')
        traces_metadata = traces_metadata(match);
    end
end

% traces = bsxfun(@minus,traces,mean(traces,1));

if plot_stack
    figure; plot_trace_stack(traces,100,zeros(length(traces),3),'-')
end


