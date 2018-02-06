%metadatastruct must have fields INsweeprange and cellno
% INsweeprange should be current clamp trials with 1s current injections

%sweeps to analyze
MD.INsweeprange=[5 6 7];
MD.cellno=1;

% find current injections
MD=findCurrentinjections(ExpStruct,MD)

%zsweeps = baseline subtracted sweeps
% get membrane time constant
[D MD]=getMembranetimeConstant(zsweeps,D,MD)

% get sample of resting vm
[D MD]= getRestingvm(sweeps,MD,D,1)               

                
% get spike properties
MD.spikethreshvalue=5;
MD.factordvdtsetting=0.01;
[D,MD]=analyzeSpikeproperties(sweeps,D,MD,1);
                                
% get input resistance from CC traces
[D,MD]=getCCinputResistance(sweeps,D,MD)
                
%cellstruct = cell1 or cell2
% calculate intrinsic props for VC trials
[D MD]=getVCintrinsicProps(zsweeps,cellstruct,D,MD,1)