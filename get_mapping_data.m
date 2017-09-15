function [maps, power_curve_num, varargout] = get_mapping_data(data,trials,varargin)

mpp_pow = [];
if ~isempty(varargin) && ~isempty(varargin{1})
    mpp_pow = varargin{1};
end

if length(varargin) > 1 && ~isempty(varargin{2})
    Fs = varargin{2};
else
    Fs = 20000;
end

% load(filename)

this_seq = cell(length(trials),1);
this_stim_key = cell(length(trials),1);
power_curve_num = cell(length(trials),1);
stim_starts = cell(length(trials),1);

full_stim_key = [];

for i = 1:length(trials)
    cur_trial = trials(i);
    this_seq{i} = data.trial_metadata(cur_trial).sequence;
    stims_per_trial(i) = length(this_seq{i});
    this_stim_key{i} = data.trial_metadata(cur_trial).stim_key;
    power_curve_num{i} = unique([this_seq{i}.target_power]);
    stim_starts{i} = [data.trial_metadata(cur_trial).sequence.start];
    for j = 1:length(this_seq{i})
        if i == 1 && j == 1
            full_seq(1) = this_seq{i}(j);
        else
            full_seq(end+1) = this_seq{i}(j);
        end
        full_seq(end).precomputed_target_index = ...
            full_seq(end).precomputed_target_index + size(full_stim_key,1);
    end
    full_stim_key = [full_stim_key; this_stim_key{i}];
end
power_curve_num = unique([power_curve_num{:}]);
maps = cell(length(power_curve_num),1);
this_seq = [this_seq{:}];
max_trial = length(this_seq);
% max_trial = 1200;
[traces_ch1,traces_ch2] = ...
    get_stim_stack(data,trials,...
        stims_per_trial,stim_starts,Fs);
stim_inds = [full_seq.precomputed_target_index];
% on_cell_trials = isnan(full_stim_key(stim_inds,1,2));
on_cell_trials = ones(size(stim_inds))';
% power_curve_num = 150;
traces = [];
stim_pow = [];
target_locs = [];
stim_inds = [];
deorder = [];
num_trials = 0;
spacing = 1;
% power_curve_num = power_curve_num(end-1:end);
for i = 1:length(power_curve_num)
    
    
%     this_seq = this_seq(1:max_trial);
    traces_pow{1} = traces_ch1(on_cell_trials' & [full_seq.target_power] == power_curve_num(i),:);
%     traces = [traces; traces_pow{1}];
%     deorder = [deorder find(on_cell_trials' & [full_seq.target_power] == power_curve_num(i))]; 
    traces_pow{2} = traces_ch2(on_cell_trials' & [full_seq.target_power] == power_curve_num(i),:);
    this_seq_power = full_seq(on_cell_trials' & [full_seq.target_power] == power_curve_num(i));
%     mpp_pow = mpp(on_cell_trials' & [full_seq.target_power] == power_curve_num(i));
%     mpp_pow = mpp(num_trials+(1:length(this_seq_power)));
    mpp_pow = [];
    num_trials = num_trials + length(this_seq_power);
    [maps{i}, mpp_maps{i}] = see_grid_multi(traces_pow,mpp_pow,this_seq_power,full_stim_key,spacing,0);
%     title(['Power = ' num2str(power_curve_num(i)) ' mW'])
%     xlim(xlims); ylim(ylims);
%     get_mean_events_image(mpp_maps{i}, 2000, 1, 1);
%     title(['Event Counts, Power = ' num2str(power_curve_num(i)) ' mW'])
%     caxis([0 2])
end

if ~isempty(mpp_pow)
    varargout{1} = mpp_maps;
end