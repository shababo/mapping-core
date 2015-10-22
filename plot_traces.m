function plot_traces(traces,default_length,offset_step)

offset = 0;

for trial = 1:size(traces,1)
    
    offset = offset - 2 * min(traces(trial,:));
    
end


trial_length = default_length;


offset = 0;

    
for trial = 1:size(traces,1)
    
    
    
    trace_to_plot = traces(trial,this_trial_start:this_trial_start+trial_length);
    
    plot((1:trial_length+1)/20000,trace_to_plot - offset - median(trace_to_plot),'LineWidth',2)
    hold on
   

    
    offset = offset + offset_step;
    
    
end


ylim([-100 50])

% title(label)
axis tight
axis off


hold off

