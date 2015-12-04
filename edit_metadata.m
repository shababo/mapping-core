function edit_metadata(filename,fields,values,trials)

load(filename)
bu_filename = [filename(1:end-4) '-bu.mat'];
save(bu_filename,'data','defaults')

for i = 1:length(fields)
    
    this_field = fields{i};
    these_trials = trials{i};
    these_values = values{i};
    
    for j = 1:length(these_trials)
        
        data.trial_metadata(these_trials(j)).(this_field) = these_values{j};
        
    end
end

save(filename,'data','defaults')