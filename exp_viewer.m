function exp_viewer( datapath , trials, type, varargin)
%EXP_VIEWER Summary of this function goes here
%   Detailed explanation goes here

% load the experiment file
load(datapath,'sweeps','data');

if exist('sweeps')
    % how many trials
    N_trials = length(sweeps);

    % data.spikes = zeros(N_trials,length(sweeps{5}));
    % data.stims = zeros(size(data.spikes));
    % data.traces = cell();

    data.sweeps = sweeps;

    % 
    % for n = 1:N_trials
    % 
    %     data.traces(n,:) = sweeps{n}(1:20000,1)';
    % %     data.spikes(n,:) = find_cell_attached_spikes(sweeps{n}(1:20000,1)',4000);
    %     data.stims(n,:) = sweeps{n}(1:20000,2)' > 30;  
    % 
    % end


    data.n_trials = N_trials;
    data.dt = 1/20000;
    % data.time = (0:length(sweeps{5})-1)*data.dt;
else
    data.n_trials = length(data.sweeps);
    data.dt = 1/20000;
end
    


if strcmp(type,'trial_by_trial')
    

    trial_by_trial_gui(data);
    
elseif strcmp(type,'detection')
    

    trial_by_trial_detection_gui(varargin{1},varargin{2});

 
elseif strcmp(type,'groups')
    
    plot_grid = varargin{1};
    plot_type = varargin{2}; % raster or traces
    titles = varargin{3};
    
    figure;
    for g = 1:length(trials)
        subplot(plot_grid(1),plot_grid(2),g)
        if strcmp(plot_type,'raster')
            plot_raster(data.spikes(trials{g},:),data.stims(trials{g},:),titles{g})
        elseif strcmp(plot_type,'traces')
            plot_trace_stack(data.traces(trials{g},:),data.stims(trials{g},:),titles{g})
        end
    end
    
elseif strcmp(type,'spatial-special')
    
    trace_ind = varargin{1};
    trace = data.sweeps{trace_ind}(:,1);
    
    thresh = varargin{2};
    event_times = find_square_events(data.sweeps{trace_ind}(:,2),thresh);
    
    event_ids = varargin{3};
    location_plot_inds = varargin{4};
    
    ids = unique(event_ids);
    event_length = 2000;
    
    all_traces = cell(length(ids),1);
    
    figure
    for g = 1:length(ids)
        
        this_id = ids(g)
        these_events = find(event_ids == this_id);
        
        these_traces = zeros(length(these_events),event_length+100);
        
        for i = 1:length(these_events)
            start_i = event_times(these_events(i));
            try
                these_traces(i,:) = trace(start_i - 100:start_i + event_length - 1);
            catch
            end
        end
            size(these_traces)
        loc = find(location_plot_inds == g);
        if ~isempty(loc)
            subplot(5,5,loc)
        else
%             subplot(5,5,find(location_plot_inds == 50))
        end
        plot_trace_stack(these_traces,zeros(size(these_traces)),[],zeros(length(these_events),3),[],event_length-1,50)
%         plot_trace_stack(traces,stims,label,colors,events,default_length,offset_step)
        all_traces{g} = these_traces;
        
    end
    
    [~, name] = fileparts(datapath);
    save([name  '_trace_' num2str(trace_ind) '_traces_' strrep(num2str(fix(clock)), ' ', '') '.mat'],'all_traces');
    
    
    
%     subplot(num_rows,num_cols,find(subplot_inds == num_rows*num_cols))
%     plot([600 1100],[-50 -50],'k','LineWidth',3)
%     hold on
%     plot([600 600],[-50 0],'k','LineWidth',3)
%     text(700,-75,'50 ms')
%     text(200,-25,'50 pA')
%     xlim([1 2000])
%     ylim([-200 0])
%     axis off

elseif strcmp(type,'spatial-special-charge')
    
    trace_ind = varargin{1};
    trace = data.sweeps{trace_ind}(:,1);
    
    thresh = varargin{2};
    event_times = find_square_events(sweeps{trace_ind}(:,2),thresh);
    
    event_ids = varargin{3};
    location_plot_inds = varargin{4};
    
    ids = unique(event_ids);
    event_length = 3000;
    
    all_traces = cell(length(ids),1);
    
    charge = zeros(7,7);
    
    figure
    for g = 1:length(ids)
        
        this_id = ids(g);
        these_events = find(event_ids == this_id);
        
        these_traces = zeros(length(these_events),event_length);
        
        for i = 1:length(these_events)
            start_i = event_times(these_events(i));
            try
                these_traces(i,:) = trace(start_i:start_i + event_length - 1);
            catch
            end
        end
            
        loc = find(location_plot_inds == g);
%         if ~isempty(loc)
%             subplot(7,7,loc)
%         else
%             subplot(7,7,find(location_plot_inds == 50))
%         end
        
        charge(loc) = -sum(sum(bsxfun(@minus,these_traces,median(these_traces,2))));
        
%         plot_trace_stack(these_traces,zeros(size(these_traces)),[],zeros(length(these_events),3),[],event_length-1)
        
%         all_traces{g} = these_traces;
        
    end
    
    imagesc(charge)
    
%     [~, name] = fileparts(datapath);
%     save([name  '_trace_' num2str(trace_ind) '_traces_' strrep(num2str(fix(clock)), ' ', '') '.mat'],'all_traces');
    
    
    
%     subplot(num_rows,num_cols,find(subplot_inds == num_rows*num_cols))
%     plot([600 1100],[-50 -50],'k','LineWidth',3)
%     hold on
%     plot([600 600],[-50 0],'k','LineWidth',3)
%     text(700,-75,'50 ms')
%     text(200,-25,'50 pA')
%     xlim([1 2000])
%     ylim([-200 0])
%     axis off
    
elseif strcmp(type,'spatial-special-events')
    
    trace_ind = varargin{1};
    trace = data.sweeps{trace_ind}(:,1)';
    
    thresh = varargin{2};
    event_times = find_square_events(data.sweeps{trace_ind}(:,2),thresh);
    
    event_ids = varargin{3};
    location_plot_inds = varargin{4};
    
    ids = unique(event_ids);
    event_length = 3000;
    
    tau = [20 100];
    
    detection_results = cell(length(ids),1);
    feature_mats = cell(length(ids),1);
    features = {'amplitude','tau1','tau2','delay'};
    
    pool = parpool(4)
    
    tic
    
    for g = 1:length(ids)
        
        this_id = ids(g);
        these_events = find(event_ids == this_id);
        
        target_results = struct();
        
        
        target_feature_mat = zeros(0,length(features));
        tmp_feature_mats = cell(length(these_events),1);
        
        

        parfor i = 1:length(these_events)
            
            start_i = event_times(these_events(i));
            trace_i = trace(start_i:start_i + event_length - 1);
            
            tGuess = find_pscs_new(trace_i, data.dt, .002, .25, 0, 0);
            disp(['Convolution Events: ' num2str(length(tGuess))]);
            
            trace_i = max(trace_i) - trace_i;
            [times, trials, mcmc, params] = sampleParams_0629tuned(trace_i,tau,tGuess,data.dt);

            errP = zeros(1,length(trials.curves));
            for j = 1:length(trials.curves)
                errP(j) = sum((trials.curves{j}-trace_i).^2);
            end

            [min_err min_err_ind] = min(errP);
            
            taus = zeros(0,2);
            disp(['Events Detected in this trace: ' num2str(length(trials.tau{min_err_ind}))])
            for j = 1:length(trials.tau{min_err_ind})
                taus = [taus; trials.tau{min_err_ind}{j}];
            end
                    
            new_features = [trials.amp{min_err_ind}' taus trials.times{min_err_ind}'];
            tmp_feature_mats{i} = new_features;
            target_feature_mat = [target_feature_mat; new_features];
            
            
            target_results(i).times = times;
            target_results(i).trials = trials;
            target_results(i).mcmc = mcmc;
            target_results(i).params = params;
            target_results(i).min_err = min_err;
            target_results(i).min_err_ind = min_err_ind;
            
        end
        
        
        
        feature_mats{g} = tmp_feature_mats;
        
        detection_results{g} = target_results;
                
        loc = find(location_plot_inds == g);
        if ~isempty(loc)
            subplot(7,7,loc)
        else
            subplot(7,7,find(location_plot_inds == 50))
        end
        
        plotmatrix(target_feature_mat);
        
    end
    toc
    delete(pool);
    
    [~, name] = fileparts(datapath);
    save([name '_trace_' num2str(trace_ind) '_detction_results_' strrep(num2str(fix(clock)), ' ', '') '.mat'],'detection_results','feature_mats');
        
    
    
    
%     subplot(num_rows,num_cols,find(subplot_inds == num_rows*num_cols))
%     plot([600 1100],[-50 -50],'k','LineWidth',3)
%     hold on
%     plot([600 600],[-50 0],'k','LineWidth',3)
%     text(700,-75,'50 ms')
%     text(200,-25,'50 pA')
%     xlim([1 2000])
%     ylim([-200 0])
%     axis off
    
elseif strcmp(type,'spatial')
    
    subplot_inds = varargin{1};
    plot_type = varargin{2};
    titles = varargin{3};
    
    
    num_rows = size(subplot_inds,1);
    num_cols = size(subplot_inds,2);
    
    subplot_inds = subplot_inds(:);
    
    data.traces = data.sweeps;
    
    figure
    for g = 1:length(trials)

        subplot(num_rows,num_cols,find(subplot_inds == g))
        if strcmp(plot_type,'raster')
            plot_raster(data.spikes(trials{g},:),data.stims(trials{g},:),titles{g})
        elseif strcmp(plot_type,'traces')
            plot_trace_stack(data.traces(trials{g},:),data.stims(trials{g},:),titles{g})
        end
    end
    
    subplot(num_rows,num_cols,find(subplot_inds == num_rows*num_cols))
    plot([600 1100],[-50 -50],'k','LineWidth',3)
    hold on
    plot([600 600],[-50 0],'k','LineWidth',3)
    text(700,-75,'50 ms')
    text(200,-25,'50 pA')
    xlim([1 2000])
    ylim([-200 0])
    axis off
    
end

    
end