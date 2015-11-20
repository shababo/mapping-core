%get_sweeps_dir(dirname,match_string,recursive,ch_ind,plot,max_sweep,varargin)

%% get all 2P chrimson
[traces_chrimson_2P, traces_metadata_chrimson_2P] = get_sweeps_dir('data',{'11_10','11_11'},0,1,1,Inf,'hologram_id','^ROI[1-2]$',...
                                               'stim_type','2P'); 


%% plot all traces
                                           
axis on; xlim([.3 .375]); title('Chrimson, All 2P')


%% get subsets for main conditions: [50, 100, 150]mW X [.05 .01 .02]sec


% trial params to iterate over
pulseamps = [25 50 100 150];
pulsedurations = [.005 .01 .02];
% trial params to keep constant (all others ignored)
params.pulsenumber = 1;

traces_by_condition = cell(length(traces_chrimson_2P),length(pulseamps),length(pulsedurations));
trace_metadata_by_condition = cell(length(traces_chrimson_2P),length(pulseamps),length(pulsedurations));
spike_times_by_condition = cell(length(traces_chrimson_2P),length(pulseamps),length(pulsedurations));



test_traces = zeros(10,20000);

for i = 1:length(traces_chrimson_2P)
    
    these_traces = traces_chrimson_2P{i};
    these_metadata = traces_metadata_chrimson_2P{i};
    
    if ~isempty(these_traces)
        for j = 1:length(pulseamps)
            for k = 1:length(pulsedurations)

                params.pulseamp = pulseamps(j);
                params.pulseduration = pulsedurations(k);

                match_inds = match_trials(params,these_metadata);
                traces_by_condition{i,j,k} = these_traces(match_inds,:);
                trace_metadata_by_condition{i,j,k} = these_metadata(match_inds);
                spike_times_by_condition{i,j,k} = cell(size(traces_by_condition{i,j,k},1),1);
                for l = 1:size(traces_by_condition{i,j,k},1)
                    spike_times_by_condition{i,j,k}{l} = get_spike_times_cell_attached(traces_by_condition{i,j,k}(l,:),.2*20000);
                end
            end
        end
        
        test_traces(i,:) = traces_by_condition{i,j,k}(1,:);
    end
end

%% test spike detection

spike_times = cell(10,1);
offset = 300;
figure
plot_trace_stack(test_traces,offset,zeros(10,3),'-');
hold on
for i = 1:10
    
    spike_times{i} = get_spike_times_cell_attached(test_traces(i,:),.2*20000);
    scatter(spike_times{i},-(i-1)*offset*ones(length(spike_times{i}),1),'r*','linewidth',2)
    hold on
    
end

%% get first spike latencies and stds

num_first_spikes_per_cond = zeros(length(traces_chrimson_2P),length(pulseamps),length(pulsedurations));
cond_by_cell_spike_time_latencies = zeros(length(pulseamps),length(pulsedurations),length(traces_chrimson_2P));
cond_by_cell_spike_time_stds = zeros(size(cell_spike_time_latencies));

max_spike_time = 

for i = 1:length(traces_chrimson_2P)
    
    for j = 1:length(pulseamps)
        for k = 1:length(pulsedurations)
            
            first_spike_times = [];
            
            for l = 1:length(spike_times_by_condition)
                
            end
        end
    end
end

%% plot all conditions for each cell

min_traces = 5;
for i = 1:length(traces_chrimson_2P)
    
    figure('position',[100 100 1500 1000])
    set(gcf,'Units','Inches');
    pos = get(gcf,'position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    i
    count = 1;
    for j = 1:length(pulsedurations)
        for k = 1:length(pulseamps)
            subplot(length(pulsedurations),length(pulseamps),count)
            if ~isempty(traces_by_condition{i,k,j})
                if size(traces_by_condition{i,k,j},1) > max_traces
                    traces_to_plot = traces_by_condition{i,k,j}(1:max_traces,:);
                else
                    traces_to_plot = traces_by_condition{i,k,j};
                end
                plot_trace_stack(traces_to_plot,50,zeros(size(traces_by_condition{i,k,j},1),3),'-')
                axis on;
                set(gca,'yticklabel','')
                xlim([.25 .45])
            else
                axis off
            end
            
            
            title([num2str(pulseamps(k)) 'mW, ' num2str(pulsedurations(j)) ' sec'])
            count = count + 1;
        end
    end
    
    saveas(gcf,['data/figures/chrimson-cell-' num2str(i) '.png'])
    close gcf
end

%% plot long runs of trials
min_traces = 100;
max_traces = 50;
for i = 1:length(traces_chrimson_2P)

    %     figure('position',[100 100 1500 1000])
%     set(gcf,'Units','Inches');
%     pos = get(gcf,'position');
%     set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    i
%     count = 1;
    for j = 1:length(pulsedurations)
        for k = 1:length(pulseamps)
%             subplot(length(pulsedurations),length(pulseamps),count)
            if ~isempty(traces_by_condition{i,k,j})
                if size(traces_by_condition{i,k,j},1) > min_traces
                    if size(traces_by_condition{i,k,j},1) > max_traces
                        traces_to_plot = traces_by_condition{i,k,j}(1:max_traces,:);
                    else
                        traces_to_plot = traces_by_condition{i,k,j};
                    end
                    figure
                    plot_trace_stack(traces_to_plot,50,zeros(size(traces_by_condition{i,k,j},1),3),'-')
                    axis on;
                    set(gca,'yticklabel','')
                    xlim([.25 .45])
                    title([num2str(pulseamps(k)) 'mW, ' num2str(pulsedurations(j)) ' sec'])
                    saveas(gcf,...
                        ['data/figures/chrimson-cell-' num2str(i) '-long-run-' num2str(pulseamps(k)) '-mW-' num2str(pulsedurations(j)) '-sec.png'])
                    close(gcf)
                end
            
            end
            
%             count = count + 1;
        end
    end
    
    
    close gcf
end


                
                
                
                


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

