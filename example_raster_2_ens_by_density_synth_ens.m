%% Example code for ensemble detection from raster
%1.- Creating synthetic raster
N=100;
T=10000;
fr=0.2;
nens = 6;
ntimesperens = [0.1 0.1 0.1 0.1 0.1 0.1]; % probability of each ensemble
ncellsperens = [10 15 20 25 30 35]; % cells per ensemble

[ensmat_in,enscells_in,raster,frates] = MakeEnsembles_fix_rate(N,fr,T,nens,ncellsperens,ntimesperens);
%% 2.- Detecting ensembles
% Algorith parameters
pars.dc = 0.02; % cut-off for distances
pars.npcs = 6;
pars.minspk = 3; % minimum 3 spikes per pattern
pars.nsur = 100; % surrogates for core-cells, should be 1000 or more.
pars.prct = 99.9;% percentile on the surrogate core-cell distribution
pars.cent_thr = 99.9;
pars.inner_corr = 5;
pars.minsize = 3;
[ensmat_out,det_core_cells] =  raster2ens_by_density(raster,pars);
det_nens = size(ensmat_out,1);
%% 3.- Comparing inputs vs output ensembles
det_ens_corr = 1-pdist2(ensmat_in,ensmat_out,'correlation');
[~,ens_id] = max(det_ens_corr,[],2);
% If any, here we look for the output ensembles that were not plugged into the raster
if det_nens>nens
    ee = 1:det_nens;
    extraens = find(~ismember(ee,ens_id));
else
    extraens =[];
end
%% 4.- Plotting input vs output ensemble sequences
tt=1:T;
figure
subplot(211)
for e=1:nens
    plot(tt(ensmat_in(e,:)>0),ensmat_in(e,ensmat_in(e,:)>0)*e,'ro');hold on
    plot(tt(ensmat_out(ens_id(e),:)>0),ensmat_out(ens_id(e),ensmat_out(ens_id(e),:)>0)*(e+0.2),'bo');hold on
end
xlabel('Time (bins)')
ylabel('Ensemble ID')
set(gca,'ytick',1:nens)
ylim([0.5 nens+0.5])

subplot(212)
nn=1:N;
for e=1:nens
    plot(nn(enscells_in(:,e)),enscells_in(enscells_in(:,e),e)*e,'ro');hold on
    plot(nn(det_core_cells(:,ens_id(e))),det_core_cells(det_core_cells(:,ens_id(e)),ens_id(e))*(e+0.2),'bo');hold on    
end
xlabel('Neuron')
ylabel('Ensemble ID')
set(gca,'ytick',1:nens)
ylim([0.5 nens+0.5])
legend('Ground Truth','Detected')

