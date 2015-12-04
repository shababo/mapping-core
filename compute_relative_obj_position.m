function compute_relative_obj_position(filename,trials)

load(filename)
bu_filename = [filename(1:end-4) '-bu.mat'];
save(bu_filename,'data','defaults')

if isempty(trials)
    
    trials = 1:length(data.trial_metadata);
    
end
    
for j = 1:length(trials)

    if isfield(data.trial_metadata,'obj_position') && isfield(data.trial_metadata,'cell_position')
        data.trial_metadata(trials(j)).relative_position = data.trial_metadata(trials(j)).obj_position - data.trial_metadata(trials(j)).cell_position;
    end
        

end


save(filename,'data','defaults')