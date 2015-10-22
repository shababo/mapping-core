function plot_raster(spikes,stims,label)

offset = 0;
stim_start = 500;
time_after_stim = 1000;


for trial = 1:size(spikes,1)
    
    offset = offset - 2 * min(spikes(trial,:));
    
end

stim_top = 2*max(spikes(1,:));
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
    trial_length = 13000;

end

offset = 0;
for trial = 1:size(spikes,1)
    
    if isempty(change_points)
        this_trial_start = 6000;
    else
        this_trial_start = find(stims(trial,:),1,'first') - stim_start;
    end
    
    plot(spikes(trial,this_trial_start:this_trial_start+trial_length) - offset,'k','LineWidth',2)
    offset = offset - 2 * min(spikes(trial,:));
    hold on
    
end

title(label)
axis tight
axis off