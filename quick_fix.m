function quick_fix(filename,trials)

load(filename)
bu_filename = [filename(1:end-4) '-bu.mat'];
save(bu_filename,'data','defaults')

if isempty(trials)
    
    trials = 1:length(data.trial_metadata);
    
end
    
cell_position = data.trial_metadata(1).obj_position + [3 0 0];

for j = 1:length(trials)
    
    ind = trials(j);
    data.trial_metadata(ind).cell_position = cell_position;
        

end


save(filename,'data','defaults')