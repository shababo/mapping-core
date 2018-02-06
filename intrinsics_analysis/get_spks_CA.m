function [spikeno spiketimes]=get_spks_CA(trace,raw_spike_thresh)

Fc=5000; 
dt=5*10^-5;
%factordvdt=0.03 ;
%win=200; win2=100;

 mean_thresh=[];
 sd_thresh=[];
 thresvalue=[];
 thres_coords_orginal=[];
 dvdt=[];
 spiketrace=[];
 
% Filtering the trace
% %Fc =[5000]; 
% [b1,a1]=butter(2,2*Fc*dt, 'low'); % high pass butterworth filter
% trace2 = filter(b1,a1,trace);  % Filtered trace to minimize noise
trace2 = highpass_filter(trace);
trace2(1:6)=trace2(7);
trace2=trace2-mean(trace2);
raw_spike_thresh = raw_spike_thresh*std(trace2);

gim=-trace2;

clear('set_crossgi')
set_crossgi=find(gim(1:end) > raw_spike_thresh)  ;  % setting the threshold
clear('index_shift_neggi');clear('index_shift_pos');


if isempty(set_crossgi) ~= 1     % This to make sure there is a spike otherwise the code below gives problems. There is an empty else statement below.
    
    clear('set_cross_plusgi');clear('set_cross_minus')
    index_shift_posgi(1)=min(set_crossgi);
    index_shift_neggi(length(set_crossgi))=max(set_crossgi);
    
    for i=1:length(set_crossgi)-1
        if set_crossgi(i+1) > set_crossgi(i)+1 ;
            index_shift_posgi(i+1)=i;
            index_shift_neggi(i)=i;
        end
    end
    
    %Identifying up and down slopes:    
    set_cross_plusgi=  set_crossgi(find(index_shift_posgi));   % find(x) returns nonzero arguments.
    set_cross_minusgi=  set_crossgi(find(index_shift_neggi));   % find(x) returns nonzero arguments.
    set_cross_minusgi(length(set_cross_plusgi))= set_crossgi(end);
    nspikes= length(set_cross_plusgi); % Number of pulses, i.e. number of windows.
    
    for i=1:nspikes
        spikemax(i)=min(find(gim(set_cross_plusgi(i):set_cross_minusgi(i)) == max(gim(set_cross_plusgi(i):set_cross_minusgi(i))))) +set_cross_plusgi(i)-1;
    end
    
%     hold on
%     plot(gim)
%     plot(spikemax,gim(spikemax),'o')
    
else
    spikemax=[];
    display('no spikes in trace')
end
spikeno=length(spikemax);
spiketimes=spikemax;
spikeheight=gim(spikemax);