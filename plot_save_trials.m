function [mean_trace, traces] =  plot_save_trials(filename,trials,tag,plot_mean,plot_traces)


traces = get_sweeps(filename, trials, plot_traces);
mean_trace = mean(traces);

full_tag = [tag '-' mat2str(trials)];

if plot_mean
    figure; 
    plot(mean_trace-mean_trace(3475)); title(full_tag)
    savefig(tag)
end