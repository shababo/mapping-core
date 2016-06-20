function match = match_trial(params,trial_metadata)

%     assignin('base','params',params)

    param_names = fieldnames(params);

    match = 1;

    for i = 1:length(param_names)
        these_params = params.(param_names{i});

        if iscellstr(these_params)
            these_params = cellstr2regexp(these_params);
        elseif iscell(these_params)
            match_tmp = 0;
            for j = 1:length(these_params)
                match_tmp = match_tmp || match_vec(these_params{j},trial_metadata,param_names{i});
            end
            match = match && match_tmp;
            continue
        end

        if ischar(these_params) && ~strcmp(these_params,'ignore')
            match = match && match_str(these_params,trial_metadata,param_names{i});
        elseif ~ischar(these_params) && isscalar(trial_metadata.(param_names{i}))
            match = match && match_vec(these_params,trial_metadata,param_names{i});
        elseif ~ischar(these_params) && isvector(trial_metadata.(param_names{i}))
            match = match && match_matrix(these_params,trial_metadata,param_names{i});
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
    
    if isfield(trial_metadata,param_name)
    
        for i = 1:length(param_set)

            this_param = param_set(i);

            if isnan(this_param)
                
                match = match || 1;

            else

                match = match || isequal(param_set(i),trial_metadata.(param_name));

            end
        end
    else
        match = 1;
    end
    
end

function match = match_matrix(param_set,trial_metadata,param_name)

    if isfield(trial_metadata,param_name)
        
        match = 0;
        metadata_param = trial_metadata.(param_name);
    
        for i = 1:size(param_set,1)

            this_vec = param_set(i,:);
            local_match = 1;

            for j = 1:length(this_vec)

                if isnan(this_vec(j))
                    
                    local_match = local_match && 1;

                else

                    local_match = local_match && isequal(this_vec(j),metadata_param(j));

                end
            end
            
            match = match || local_match;
        end
    else
        match = 1;
    end
    
end
