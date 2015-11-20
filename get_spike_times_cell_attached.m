function [spiketimes ] = get_spike_times_cell_attached( trace , start_ind)
% Extracts spikes times for a given sweep.
% Threshold is given in absolute pA or mV. The defaults 
%   Detailed explanation goes here




% high pass filter the data
thissweep=highpass_filter(trace);

% get absoulte value of threshold since sweep will be inverted
threshold = 5*std(thissweep(1:.045*20000));

% invert sweep 
thissweep = -thissweep;

% core line for extracting spike times
% set minimum time distance between found spike times to avoid multiple
% points per any given spike (default is 0.003 s)
if (~isempty(find(thissweep>threshold)));
[y,spiketimes]=findpeaks(thissweep, 'minpeakheight',threshold,'minpeakdistance',20000*0.003);
spiketimes = spiketimes';
spiketimes = spiketimes(spiketimes > start_ind);
spiketimes = spiketimes/20000;


else
    
    spiketimes = [];



end

