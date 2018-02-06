function [datastruct metadatastruct]=analyzeSpikeproperties(sweeps,datastruct,metadatastruct,varargin)

if isempty(varargin{1})
view=0;
else
view=1;
end

cellno=metadatastruct.cellno;
iStart=metadatastruct.iStart;
iAmps=metadatastruct.iAmps;
spikethreshvalue = metadatastruct.spikethreshvalue;
factordvdtsetting= metadatastruct.factordvdtsetting;

try
    metadatastruct.spikedetectionflag=0;
    for i = 1:length(metadatastruct.INsweeprange)
    if view
    figure
    hold on
    end
        for j = 1:length(iStart{i})
            try
                thisStart=(iStart{i}(j))*10-4010;
                thisEnd=thisStart+32000;
               
                thisTrace=sweeps{metadatastruct.INsweeprange(i)}(thisStart:thisEnd,cellno);
                currInjtraces{i,j}=thisTrace;
                if iAmps{i}(j)>0
                    [thisThreshvalues, thisThreshcoords, datastruct.spikeheight{i,j}, datastruct.spikePeaktimes{i,j},datastruct.halfMax{i,j}, FWHM{i,j}, datastruct.halfMaxtimeRise_original{i,j},...
                        datastruct.halfMaxtimeFall_original{i,j}, datastruct.returnTothresh_original{i,j}, datastruct.returnToposVmderiv_original{i,j}, datastruct.widespiketrace{i,j},...
                        datastruct.widedvdt{i,j}, dvdt{i,j}, spiketrace{i,j}]=analyze_spikes(thisTrace,60,50,factordvdtsetting,spikethreshvalue);
                    thisThreshtimes=(thisThreshcoords-thisStart)/20;
                    spikeOnsettimes{i,j}=thisThreshtimes;
                    threshValues{i,j} =thisThreshvalues;                    
                end
                if view
                hold on
                timebase=((1:length(thisTrace))+thisStart);
                plot(timebase,thisTrace)
                    for k = 1:length(thisThreshcoords)
                        try
                        endpoint=min(datastruct.returnToposVmderiv_original{i,j}(k),thisThreshcoords(k)+100);
                        thisTimes=(thisThreshcoords(k):endpoint);
                        end
                        try
                        plot(thisTimes+thisStart,thisTrace(thisTimes));
                        catch
                        print derp
                        end
                    end
               
                plot(thisThreshcoords+thisStart,thisTrace(thisThreshcoords),'or')
                plot(datastruct.halfMaxtimeRise_original{i,j}+thisStart,...
                thisTrace(datastruct.halfMaxtimeRise_original{i,j}),'xk')
%                 plot(datastruct.halfMaxtimeFall_original{i,j}+thisStart,...
%                 thisTrace(datastruct.halfMaxtimeFall_original{i,j}),'xk')
                plot(datastruct.spikePeaktimes{i,j}+thisStart,thisTrace(datastruct.spikePeaktimes{i,j}),'og')
                end

            end
        end
    end
    allThreshvalues=cell2mat(threshValues(:)');
    mThresh=mean(allThreshvalues);
    allFWHM=cell2mat(FWHM(:)');
    mFWHM=mean(allFWHM);
    datastruct.allFWHM=allFWHM;
    datastruct.FWHM=FWHM;
    datastruct.mFWHM=mFWHM;
    datastruct.allThreshvalues=allThreshvalues;
    datastruct.mThresh=mThresh;
    datastruct.spikeOnsettimes=spikeOnsettimes;
    datastruct.threshValues=threshValues;
    datastruct.spiketrace=spiketrace;
    datastruct.dvdt=dvdt;
    datastruct.currInjtraces=currInjtraces;
catch
    metadatastruct.spikedetectionflag=1;
end

try
    fFWHMs=[];
    for j = 1:size(FWHM,1)
        for k = 1:size(FWHM,2)
            if ~isempty(FWHM{j,k})
                fFWHMs=[fFWHMs FWHM{j,k}(1)];
            end
        end
    end
    
    mfFWHM=nanmean(fFWHMs);
    medfFWHM=nanmedian(fFWHMs);
    datastruct.fFWHMs=fFWHMs;
    datastruct.mfFWHM=mfFWHM;
    datastruct.medfFWHM=medfFWHM;
end


end