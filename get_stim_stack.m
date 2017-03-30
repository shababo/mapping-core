function [stack_ch1,stack_ch2] = get_stim_stack(data,trials,num_stims)

traces_ch1 = cell(1,length(trials));
traces_ch2 = cell(1,length(trials));
for i = 1:length(trials)
    
    trial_ind = trials(i); 
    traces = [data.sweeps{trial_ind}(:,1)'; data.sweeps{trial_ind}(:,2)'];
    stim = data.sweeps{trial_ind}(:,4)' > .05; sum(diff(stim) == 1)
    stim_starts = find(diff(stim) == 1);
    if length(stim_starts) ~= num_stims
        itis = diff(stim_starts);
        early_fails = find(itis > median(itis),length(itis),'first');
        stim_starts(early_fails) = [];
        itis = diff(stim_starts);
        late_fails = find(itis > median(itis),length(itis),'last');
        stim_starts(late_fails) = [];
    end
%     
    stacks = build_stim_stack_multi(traces,stim_starts,.1*20000);
    traces_ch1{i} = stacks{1};
    traces_ch2{i} = stacks{2};
end

stack_ch1 = stack_stacks(traces_ch1);
stack_ch2 = stack_stacks(traces_ch2);
