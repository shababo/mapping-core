function [traces_offset, offsets]= get_trace_stack(traces,default_length,offset_step,down_sample_rate)

offset = 0;
stim_start = 1;
time_after_stim = 1; %1000

for trial = 1:size(traces,1)
    
    offset = offset - 2 * min(traces(trial,:));
    
end



    trial_length = default_length;



offset = 0;

    
traces_offset = [];
offsets = [];
for trial = 1:size(traces,1)
    
    this_trial_start = stim_start;
    
    trace_to_plot = traces(trial,this_trial_start:this_trial_start+trial_length);
    
    traces_offset = [traces_offset; downsample(trace_to_plot - offset - trace_to_plot(1),down_sample_rate)];
    offsets = [offsets  -offset];
   
    offset = offset + offset_step;
    
    
end
