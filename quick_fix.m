function quick_fix(filename,bad_trials)

load(filename)
bu_filename = [filename(1:end-4) '-bu.mat'];
save(bu_filename,'data','defaults')

% if isempty(trials)
    
    trials = 1:length(data.trial_metadata);
    
% end
    

for j = 1:length(trials)
    
    
    is_good_trial = ~any(bad_trials == j);
    data.trial_metadata(j).is_good_trial = is_good_trial;
        

end


save(filename,'data','defaults')