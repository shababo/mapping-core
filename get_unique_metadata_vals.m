function unique_values = get_unique_metadata_vals(metadata,param_name)

if ischar(metadata(1).(param_name))
    disp('here')
    unique_values = {};
    
    for i = 1:length(metadata)
    
        match = 0;
        
        for j = 1:length(unique_values)
            
            if isequal(metadata(i).(param_name),unique_values{j})
                
                match = 1;
                
            end
        end
        
        if ~match
            unique_values{end + 1} = metadata(i).(param_name);
        end
        
    end
        
else
    disp('here')
    unique_values = [];
    
    for i = 1:length(metadata)
    
        match = 0;
        
        for j = 1:size(unique_values,1)
            
            if isequal(metadata(i).(param_name),unique_values(j,:))
                
                match = 1;
                
            end
        end
        
        if ~match
            unique_values(end+1,:) = metadata(i).(param_name);
        end
        
    end
        
end

