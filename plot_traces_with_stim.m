function plot_traces_with_stim(traces,stims, offset_step, barsize)

offset = 0;

for trial = 1:size(traces,1)
    
    offset = offset - 2 * min(traces(trial,:));
    
end

% stim_top = max(traces(1,:)) - offset - median(traces(1,:)) + 20


num_stims = length(find(diff(stims)))/2;
change_points = find(diff(stims));

if ~isempty(change_points)


    for i = 1:num_stims

        plot(change_points((i-1)*2+1:(i-1)*2+2),[1 1]*75,'b','Linewidth',4)
        hold on
    end

end



% scatter(find(stims),max(max(traces))*ones(1,length(find(stims))),'.b')
hold on

offset = 0;
for trial = 1:size(traces,1)
    
    
    trace_to_plot = traces(trial,:);
    

    plot(trace_to_plot - offset - median(trace_to_plot),'LineWidth',2,'Color',[0 0 0])
    hold on
   
    
    offset = offset + offset_step;
    
    
end

plot(floor(length(trace_to_plot)/5)*[1 1], -[offset offset + barsize],'Linewidth',2,'Color',[0 0 0])
hold on
plot(floor(length(trace_to_plot)/5)*[1 1] + [0 1000], -[offset + barsize offset + barsize],'Linewidth',2,'Color',[0 0 0])

% ylim([-200 0])

% title(label)
axis tight
axis off


hold off

