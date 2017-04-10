function charge_map = get_charge_map(traces)

charge_map = zeros(size(traces));

for i = 1:size(traces,1)
    for j = 1:size(traces,1)
        
        if ~isempty(traces{i,j})
            charge_map(i,j) = mean(sum(bsxfun(@minus,traces{i,j},median(traces{i,j},2)),2));
        end
            
        
    end
end