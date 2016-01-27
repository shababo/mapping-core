function match = match_trial(params,trial_metadata)

    assignin('base','params',params)

    param_names = fieldnames(params);

    match = 1;

    for i = 1:length(param_names)

        these_params = params.(param_names{i});
        if iscellstr(these_params)
            these_params = cellstr2regexp(these_params);
        elseif iscell(these_params)
            disp('in')
            match_tmp = 0;
            for j = 1:length(these_params)
                match_tmp = match_tmp || match_vec(these_params{j},trial_metadata,param_names{i});
            end
            match = match && match_tmp;
            continue
        end

        if ischar(these_params) && ~strcmp(these_params,'ignore')
            match = match && match_str(these_params,trial_metadata,param_names{i});
        elseif ~ischar(these_params) && isvector(these_params)
            match = match && match_vec(these_params,trial_metadata,param_names{i});
        end

    end
    
end

function match = match_str(param_set,trial_metadata,param_name)

    match = 0;
    
    if strcmp(param_set,'ignore')
        match = 1;
        return
    end

    if isfield(trial_metadata,param_name)

        inds = logical(regexpi(trial_metadata.(param_name),param_set));
        match = ~isempty(inds);
        
    end


end


function match = match_vec(param_set,trial_metadata,param_name)

    match = 0;
    
    if isnan(param_set)
        match = 1;
        return
    end

    if isfield(trial_metadata,param_name)
        
        match = isequal(param_set,trial_metadata.(param_name));

    end
    
end
