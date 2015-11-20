function regexp_str = cellstr2regexp(in_cellstr)

tmp = ['(^' in_cellstr{1}];
for j = 2:length(in_cellstr)
    tmp = [tmp '.*\.mat$)|(' in_cellstr{j}];
end
tmp = [tmp '.*\.mat$)'];
regexp_str = tmp;