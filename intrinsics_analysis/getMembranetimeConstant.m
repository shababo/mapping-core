function [datastruct metadatastruct]= getMembranetimeConstant(zsweeps,datastruct,metadatastruct)

try
allCCtau=[];
cellno=metadatastruct.cellno;
iStart=metadatastruct.iStart;
metadatastruct.CCtauflag=0;
for i = 1:length(metadatastruct.INsweeprange)
    for j = 1:length(iStart{i})
        endStart=(iStart{i}(j)+2000)*10;
        endEnd=endStart+200;
        tautrace=zsweeps{cellno,metadatastruct.INsweeprange(i)}(endStart:endEnd);
        timebase=((0:endEnd-endStart)/20000)';
        [tauObj]=fit(timebase,tautrace,'exp2');
        %                             figure
        %                             plot(tauObj); hold on; plot(timebase,tautrace)
        output=coeffvalues(tauObj);
        etau{i}(j)=-1/min(output);
        startStart=(iStart{i}(j))*10-5;
        startEnd=startStart+200;
        tautrace=zsweeps{cellno,metadatastruct.INsweeprange(i)}(startStart:startEnd);
        timebase=((0:startEnd-startStart)/20000)';
        [tauObj]=fit(timebase,tautrace,'exp2');
        %                                             figure
        %                                             plot(tauObj); hold on; plot(timebase,tautrace)
        output=coeffvalues(tauObj);
        stau{i}(j)=-1/min(output);

    end
    tau{i}=[stau{i};etau{i}];
    allCCtau=[allCCtau; tau{i}(:)];

end


allCCtau=allCCtau(isfinite(allCCtau));
meanCCtau=mean(allCCtau);
medCCtau=median(allCCtau);
datastruct.allCCtau=allCCtau;
datastruct.meanCCtau=meanCCtau;
datastruct.stau=stau;
datastruct.etau=etau;
datastruct.tauCC=tau;
datastruct.medCCtau=medCCtau;
catch
metadatastruct.CCtauflag=1;
end
                