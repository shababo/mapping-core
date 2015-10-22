function plot_trace_stack_multi(traces,stims,label,colors,events)



offset = 0;
stim_start = 1;
time_after_stim = 1; %1000
offset_step = 50;

for trial = 1:size(traces,1)
    
    offset = offset - 2 * min(traces(trial,:));
    
end

stim_top = 2*max(traces(1,:));
stim_bottom = -offset;

num_stims = length(find(diff(stims(1,:))))/2;
change_points = find(diff(stims(1,:)));

if ~isempty(change_points)

    trial_length = change_points(end) - change_points(1) + stim_start + time_after_stim;

    for i = 1:num_stims

        stim_length = change_points(2*i) - change_points(2*i - 1);
        this_start = change_points(2*i-1)- change_points(1) + stim_start;
        rectangle('Position', [this_start stim_bottom stim_length stim_top-stim_bottom],'FaceColor','b','EdgeColor','b')
        hold on
    end
else
    trial_length = 2999;

end

offset = 0;
for trial = 1:size(traces,1)
    
    if isempty(change_points)
        this_trial_start = stim_start;
    else
        this_trial_start = find(stims(trial,:),1,'first') - stim_start;
    end
    
    for trace_i = 1:size(traces,3)
        plot((1:trial_length+1)/10000,traces(trial,this_trial_start:this_trial_start+trial_length,trace_i) - offset,'LineWidth',2,'Color',colors(ceil(trial/2),:,trace_i))
        hold on
    end
    
    if ~isempty(events)
        scatter((events{trial} - stim_start)/10000,(traces(trial,this_trial_start) - offset + offset_step/3)*ones(size(events{trial})),[],colors(ceil(trial/2),:),'filled')
        hold on
    end
    
    offset = offset + offset_step;
    
    
end


ylim([-200 0])

% title(label)
axis tight
axis off


hold off

