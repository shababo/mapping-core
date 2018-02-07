function [datastruct metadatastruct]=getVCintrinsicProps(zsweeps,cellstruct,datastruct,metadatastruct,varargin)


if isempty(varargin{1})
    view=0;
else
    view=1;
end



cellno=metadatastruct.cellno;

try
    metadatastruct.VCintrinsicsflag=0;
    IRsweeprange=find(cellstruct.series_r<30);

    try
        IRsweeprange=IRsweeprange(cellstruct.rectype<2);
    end
    
    
    allCap=[];
    allIR=[];
    VCtau=[];
    Fs=20000;
    
    for i = 1:length(IRsweeprange)
        thisTrace=zsweeps{cellno,IRsweeprange(i)};
        subTrace = thisTrace(0.10015*Fs:0.12*Fs) - median(thisTrace(0.12*Fs:0.15*Fs));
        capacitanceByCharge = (trapz(subTrace)/Fs)/4; % in nanoFarards
        Ir = -4/(mean(thisTrace((Fs*.080):(Fs*.095)))-mean(thisTrace((Fs*.01):(Fs*.02))));
        tauM = capacitanceByCharge*Ir;
        allCap(i)=capacitanceByCharge;
        allIR(i)=Ir;
        VCtau(i)=tauM;
    end
    
    goodCap=find(allCap>0);
    allIR=allIR(goodCap);
    VCtau=VCtau(goodCap);
    allCap=allCap(goodCap);
    meanIR=mean(allIR);
    meanCap=mean(allCap);
    meanVCtau=mean(VCtau);
    if view
        figure
        subplot(131)
        hist(allCap)
        xlabel('capacitance')
        subplot(132)
        hist(allIR*1000)
        xlabel('input resistance')
        subplot(133)
        plotSpread({VCtau},'showMM',1)
        ylabel('tau')
    end
    medVCtau=median(VCtau);
    datastruct.allCap=allCap;
    datastruct.allIR=allIR;
    datastruct.VCtau=VCtau;
    datastruct.meanVCtau=meanVCtau;
    datastruct.meanIR=meanIR;
    datastruct.meanCap=meanCap;
    datastruct.medVCtau=medVCtau;
    metadatastruct.IRsweeprange=IRsweeprange;
    
catch
    metadatastruct.VCintrinsicsflag=1;
end

end
