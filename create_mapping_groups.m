function [groups, names] = create_mapping_groups(start_inds, run_length);

groups = cell(run_length,1);
names = cell(run_length,1);

for i = 1:run_length
    
    groups{i} = [start_inds + i - 1];
    names{i} = num2str(i);
    
end
    