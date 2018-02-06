function [datastruct metadatastruct]=getCCinputResistance(sweeps,datastruct,metadatastruct,varargin)


cellno=metadatastruct.cellno;
iStart=metadatastruct.iStart;
iAmps=metadatastruct.iAmps;

    try
        metadatastruct.CCinputresistanceflag=0;
        allCCIR=[];
        posCCIR=[];
        CCIR=cell(1);
        meanCCIR=[];
        for i = 1:length(metadatastruct.INsweeprange)
            try
                for j = 1:length(iStart{i})
                    if iAmps{i}(j)<0;
                        %if negative inj, get vm from +700 ms to+950
                        injStart=(iStart{i}(j)+700)*10;
                        injEnd=(iStart{i}(j)+950)*10;
                        injvm=sweeps{metadatastruct.INsweeprange(i)}(injStart:injEnd,cellno);
                        baseStart=(iStart{i}(j)-50)*10;
                        baseEnd=(iStart{i}(j)-1)*10;
                        baselinevm=sweeps{metadatastruct.INsweeprange(i)}(baseStart:baseEnd,cellno);
                        CCIR{i}(j)=(median(injvm)-median(baselinevm))/iAmps{i}(j);
                    elseif iAmps{i}(j)>0 && isempty(spikeOnsettimes{i})
                        injStart=(iStart{i}(j)+700)*10;
                        injEnd=(iStart{i}(j)+950)*10;
                        injvm=sweeps{metadatastruct.INsweeprange(i)}(injStart:injEnd,cellno);
                        baseStart=(iStart{i}(j)-50)*10;
                        baseEnd=(iStart{i}(j)-1)*10;
                        baselinevm=sweeps{metadatastruct.INsweeprange(i)}(baseStart:baseEnd,cellno);
                        posCCIR{i}(j)=(median(injvm)-median(baselinevm))/iAmps{i}(j);
                    end
                end

            end
        end
        allCCIR=cell2mat(CCIR);
        datastruct.allCCIR=allCCIR;
        datastruct.CCIR=CCIR;
        datastruct.posCCIR=posCCIR;
        meanCCIR=nanmean(allCCIR);
        datastruct.meanCCIR=meanCCIR;
    catch
        metadatastruct.CCinputresistanceflag=1;
    end
