function [trace_array, metadata_array] = build_trace_array(filename,run,trial_grid)

load(filename)

trace_array = cell(size(trial_grid));
metadata_array = cell(size(trial_grid));

for i = 1:size(trace_array,1)
    for j = 1:size(trace_array,2)
        
        
    
        [traces, traces_metadata] = get_sweeps_dir('.',filename,0,1,0,Inf,'run_count',run);
        size(traces)
        traces = traces{1};
    
        params1.run_count = run;
        match_inds = match_trials(params1, traces_metadata{1});
        traces = traces(match_inds,:);
        temp = traces_metadata{1};
        traces_metadata = temp(match_inds);
        
        params2.relative_position = trial_grid{i,j};
        match_inds = unique([match_trials(params2, traces_metadata)]);

        if isempty(match_inds)
            'none found'
        end
        
        trace_array{i,j} = traces(match_inds,:);
        metadata_array{i,j} = traces_metadata(match_inds);
        
    end
end