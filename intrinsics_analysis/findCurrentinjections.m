function [metadatastruct] =findCurrentinjections(ExpStruct,metadatastruct)
%metadatastruct must have fields INsweeprange and cellno

iStart=cell(1);
cellno=metadatastruct.cellno;
for i = 1:length(metadatastruct.INsweeprange)
    iAmps{i}=unique(ExpStruct.stims{metadatastruct.INsweeprange(i)}{cellno+1})*400;
    iAmps{i}=iAmps{i}(iAmps{i}~=0);
    for j = 1:length(iAmps{i})
        tol = eps;
        tempInd=find(abs(ExpStruct.stims{metadatastruct.INsweeprange(i)}{cellno+1}-iAmps{i}(j)/400) < tol*3);
        iStart{i}(j)=tempInd(1);
    end
end

metadatastruct.iAmps=iAmps;
metadatastruct.iStart=iStart;