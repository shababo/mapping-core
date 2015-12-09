function quick_fix(filename,trials)

load(filename)
bu_filename = [filename(1:end-4) '-bu.mat'];
save(bu_filename,'data','defaults')

if isempty(trials)
    
    trials = 1:length(data.trial_metadata);
    
end
    
for j = 1:length(trials)
    
    ind = trials(j);
    temp = data.trial_metadata(ind).relative_position;
    data.trial_metadata(ind).relative_position = [temp(1:2) 0];
        

end


save(filename,'data','defaults')