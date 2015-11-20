function match_inds = match_trials(params, trial_metadata)

match_inds = [];
for i = 1:length(trial_metadata)
    
    if match_trial(params, trial_metadata(i))
        match_inds = [match_inds i];
    end
    
end