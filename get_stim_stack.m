function [traces_ch1,traces_ch2] = get_stim_stack(data,trials,num_stims)

traces_ch1 = cell(1,length(trials));
traces_ch2 = cell(1,length(trials));
for i = 1:length(trials)
    
    trial_ind = trials(i); 
    traces = [data.sweeps{trial_ind}(:,1)'; data.sweeps{trial_ind}(:,2)'];
    stim = data.sweeps{trial_ind}(:,4)' > .05; sum(diff(stim) == 1)
    stim_starts = find(diff(stim) == 1);
    if stim_starts > num_stims
        while length(stim_starts) > num_stims
            disp('in while')
            length(stim_starts)
            num_stims
            itis = diff(stim_starts);
            assignin('base','itis',itis)
    %         early_fails = find(abs(itis - median(itis)) > 40,20,'first');
    %         if early_fails(1) == 1
    %             early_fails(1:find(diff(early_fails) ~= 1,1,'first')) = [];
    % %             early_fails(early_fails > length(itis)/2) = [];
    %             stim_starts(early_fails) = [];
    %         end
    %         itis = diff(stim_starts);
    %         late_fails = find(abs(itis - median(itis)) > 40,20,'last');
    %         if late_fails(end) == length(itis)
    %             late_fails(1:find(diff(late_fails) ~= 1,1,'first')) = [];
    %             late_fails(late_fails < length(itis)/2) = [];
    %             stim_starts(late_fails+1) = [];
    %         end

             bad_itis = find(abs(itis - median(itis)) > 10);
             if ~isempty(bad_itis)
                 if bad_itis(1) == 1 && bad_itis(2) ~= 2
                     stim_starts(1) = [];
                 else
                     stim_starts(bad_itis(1)+1) = [];
                 end

             end
%              return
        end
    end
    stacks = build_stim_stack_multi(traces,stim_starts,.1*20000);
    traces_ch1{i} = stacks{1};
    traces_ch2{i} = stacks{2};
end

traces_ch1 = stack_stacks(traces_ch1);
traces_ch2 = stack_stacks(traces_ch2);
