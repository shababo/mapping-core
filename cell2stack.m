function traces = cell2stack(traces_in)

traces = zeros(length(traces_in),length(traces_in{1}));

for i = 1:length(traces_in)
    traces(i,:) = traces_in{i};
end