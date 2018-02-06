function [thresvalue, thres_coords_orginal, spikeheight, spikePeaktimes, halfMax, fullWidthhalfMax, halfMaxtimeRise_original,...
    halfMaxtimeFall_original, returnTothresh_original, returnToposVmderiv_original, widespiketrace,...
    widedvdt, dvdt, spiketrace]=analyze_spikes(trace,win, win2, factordvdt,raw_spike_thresh)

%% Adapted from Rune Berg 2015 by Alex Naka 2016

%%  Spike threshold detection alogrithm - Method I
%   This function calculates the spike threshold of a membrane potential
%   recording based on the method by Sekerli et al IEEE Trans Biom Eng
%   51(9): 1665-1672, 2004. In this paper they suggest two methods of
%   estimating threshold. Here we use the first method, which is based on
%   revealing the threshold in a phase plot of V versus the derivative
%   dVdt. The largest positive slope is an indication of the threshold. 


%%  INPUT:
%   trace: Vm recording containing the spikes for which the threshold to
%   determine
%   
%   win: the half-window size in data points around the spike to show, default is
%   200 points.
%   win2: the window from the peak backwards in time in which the threshold
%   is found. Default value is 100 datapoints
%   factordvdt: is the factor times the max derivative in Vm, as a minimum to
%   include in the analysis. This is to limit the divergent data points.
%   Default is 0.03
%   default settings:
%   dt: sampling interval in seconds
%   Fc: cut-off frequency for low pass filtering, default value is 5000 Hz
%%  OUTPUT
%   mean_thresh: mean of threshold values
%   sd_thresh: standard deviation of threshold values
%   N: number of spikes
%   thresvalue: the threshold values
%   thres_coords_orginal: the x-cords of the filtered trace, which is
%   approximately the same as the original.

%   Sample :  [mean_thresh sd_thresh N thresvalue thres_coords_orginal]=Spike_threshold_PS(trace,5*10^-4,100, 100, 5000,0.03)

%   Rune W. Berg 2015
%   University of Copenhagen
%   rune@berg-lab.net, runeb@sund.ku.dk
%   www.berg-lab.net

%Putting in some defaults

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
%Fc =[5000]; 
[b1,a1]=butter(2,2*Fc*dt, 'low'); % low pass butterworth filter
trace2 = filter(b1,a1,trace);  % Filtered trace to minimize noise
trace2(1:6)=trace2(7);



gim=trace;

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
N=length(spikemax);
spikecords=spikemax;
spikeheight=gim(spikemax);


base=1:2*win+1 ; g=zeros(N,2*win+1);
newspikecords=[];
for i = 1:N
    if spikecords(i)-win>0 &spikecords(i)+win<length(trace2)
        newspikecords=[newspikecords spikecords(i)];
    end
end

spikecords=newspikecords;
N=length(spikecords);

spikePeaktimes=spikecords;

for i=1:N
    
    widespiketrace(i,:)=trace2(spikecords(i)-win:spikecords(i)+win);     % Selecting trace around spike
    widedvdt(i,:)=[diff(widespiketrace(i,:),1) 0];                          % First derivative
    setNZ=find(widedvdt(i,:)>factordvdt*max(widedvdt(i,:)));       % This value limits the divergent points
    d2vdt2(i,:)=[diff(widespiketrace(i,:),2) 0 0];                      % second derivative
    g(i,setNZ)=d2vdt2(i,setNZ)./widedvdt(i,setNZ);                      % Metric for which the maximum indicate threshold
    setNZwin2=find(setNZ<win & (setNZ>win-win2))+win-win2 ;         % points in win2 before peak of spike
    thres_coord(i)=min(find(g(i,win-win2:win) == max(g(i,win-win2:win))))+win-win2-1 ;
    thresvalue(i)=widespiketrace(i,thres_coord(i));
    thres_coords_orginal(i)=spikecords(i)-win+thres_coord(i)-1;
    
    halfMax(i)=thresvalue(i)+(spikeheight(i)-thresvalue(i))/2;   
    [trash halfMaxtimeRise(i)]=min(abs(trace2(thres_coords_orginal(i):spikecords(i))-halfMax(i)));    
    [trash halfMaxtimeFall(i)]=min(abs(trace2(spikecords(i):spikecords(i)+40)-halfMax(i)));        
    halfMaxtimeRise_original(i)=thres_coords_orginal(i)+halfMaxtimeRise(i)-1;
    halfMaxtimeFall_original(i)=spikecords(i)+halfMaxtimeFall(i)-1;
    fullWidthhalfMax(i)=halfMaxtimeFall_original(i)-halfMaxtimeRise_original(i);    
    
    returnTothresh(i)=nan;
    if min(trace2(halfMaxtimeFall_original(i):spikecords(i)+win)-thresvalue(i))
        [trash returnTothresh(i)]=min(abs(trace2(halfMaxtimeFall_original(i):spikecords(i)+win)-thresvalue(i)));
        returnTothresh_original(i)=halfMaxtimeFall_original(i)+returnTothresh(i)-1;
    end
        returnToposVmderiv(i)=nan;
        vmdiff=diff(trace2(halfMaxtimeFall_original(i):spikecords(i)+win));
        try
        returnToposVmderiv(i)=min(find(vmdiff>0));
        catch
        returnToposVmderiv(i)=spikecords(i)+win;    
        end
        returnToposVmderiv_original(i)=halfMaxtimeFall_original(i)+returnToposVmderiv(i)-1;
     
        spiketrace{i}=trace2(thres_coords_orginal(i):returnTothresh_original(i));
        dvdt{i}=diff(trace2(thres_coords_orginal(i)-1:returnTothresh_original(i)));
       
end




%%
% figure(203)
% subplot(131);plot(spiketrace','color',[.7 .7 .7]);hold on
% for i=1:N
%     plot(thres_coord(i),spiketrace(i,thres_coord(i)),'.r')
% end
% plot(mean(spiketrace,1)','b');hold off; grid on
% xlabel('datapoints')
% subplot(132);plot(spiketrace',dvdt','color',[.7 .7 .7]);hold on
% plot(mean(spiketrace,1)',mean(dvdt,1)','r');hold off; grid on
% xlabel('V'); ylabel('dVdt')
% subplot(133);hist(thresvalue); xlabel('thresholds');
% 
% figure(205)
% plot(trace);hold on
% plot(thres_coords_orginal,trace2(thres_coords_orginal),'.r')
% plot([1 length(trace)], [mean_thresh mean_thresh],'k')
% plot([1 length(trace)], [mean_thresh-sd_thresh mean_thresh-sd_thresh],'--k')
% plot([1 length(trace)], [mean_thresh+sd_thresh mean_thresh+sd_thresh],'--k')
% hold off
end
