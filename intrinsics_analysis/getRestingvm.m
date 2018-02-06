function [datastruct metadatastruct]= getRestingvm(sweeps,metadatastruct,datastruct,varargin)   

if isempty(varargin{1})
    view=0;
else
    view=1;
end
try
    metadatastruct.Vmflag=0;
    cellno=metadatastruct.cellno;
    for i = 1:length(metadatastruct.INsweeprange)
        thiszeroInd=1:10000;
        Vm{i}=sweeps{metadatastruct.INsweeprange(i)}(thiszeroInd,cellno)';
    end
    allVm=cell2mat(Vm);
    restingVm=median(allVm);

    if view
        figure
        hist(allVm,-90:0.1:-20)
    end

    datastruct.restingVm=restingVm;
    datastruct.allVm=allVm;
catch
    metadatastruct.Vmflag=1;
end