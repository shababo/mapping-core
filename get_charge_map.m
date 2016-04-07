function charge_map = get_charge_map(traces)

charge_map = zeros(size(traces));

for i = 1:size(traces,1)
    for j = 1:size(traces,1)
        
        charge_map(i,j) = mean(sum(bsxfun(@minus,traces{i,j},traces{i,j}(:,1)),2));
        
    end
end