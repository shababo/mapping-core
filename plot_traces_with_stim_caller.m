function plot_traces_with_stim_caller(filepath,trace_inds,start_ind,end_ind,spacing,barsize)


load(filepath,'sweeps')

data = sweeps(trace_inds);
response = zeros(length(data),length(data{1}(start_ind:end_ind,1)'));

for i = 1:length(data), response(i,:) = data{i}(start_ind:end_ind,1)'; end
stim = data{1}(start_ind:end_ind,2)' > 25;
% response = bsxfun(@minus,response,mean(response,1));

figure; plot_traces_with_stim(response,stim,spacing,barsize)

