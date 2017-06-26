function trial_metadata = revise_trial_metadata(filename,do_bu,varargin)

% load data
load(filename)

% bu data
if do_bu
    [pathstr,name,ext] = fileparts(filename);
    if isempty(pathstr)
        bu_name = [name '_bu' ext];
    else
        bu_name = [pathstr '/' name '_bu' ext];
    end
    save(bu_name,'data','defaults')
end

% edit metadata

trial_metadata = data.trial_metadata;

num_params_to_edit = length(varargin)/2;
num_trials = length(trial_metadata);

for i = 1:num_params_to_edit
    
    param_name = varargin{2*(i-1)+1};
    param_val = varargin{2*i};
    
    for j = 1:num_trials
        
        if length(param_val) == num_trials
            trial_metadata(j).(param_name) = param_val{j};
        else
            trial_metadata(j).(param_name) = param_val;
        end
    end
end
data.trial_metadata = trial_metadata;
save(filename,'data','defaults')   
