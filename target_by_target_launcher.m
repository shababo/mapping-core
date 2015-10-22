function target_by_target_launcher(traces_path,features_path)

load(traces_path)

data.traces = all_traces;

load(features_path)

data.feature_mats = feature_mats;

data.n_targets = length(feature_mats);
            %propose a uniform add

data.results = detection_results;

target_by_target_gui(data)