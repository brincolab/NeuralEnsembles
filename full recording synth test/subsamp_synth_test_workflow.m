%% Testing the ensembles detection method using subsampling
% 1.- Generating synthetic data
T = 50000;% recording time
N = 300; % number of neuron
nens = 10; % number of ensembles
fr = 0.2; % maximal neuron firing rate
ncellsperens = repmat(N*0.05,[1 nens]); % cells per ensemble
ntimesperens = repmat((0.7)/nens,[1 nens]); % probability of each ensemble
[ensmat_in,enscells_in,raster,frates] = MakeEnsembles_fix_rate(N,fr,T,nens,ncellsperens,ntimesperens);

%% Detecting using subsampling
pars.sr = 1; % Hz
pars.bin = 0.1; % seconds
pars.dc = 0.02; % cut-off for distances
pars.ishalo = 1; % not using halo refinement of the clusters
pars.npcs = 6;
pars.maxmem = 10^9; % maximal Gb of RAM 
pars.minspk = 3; % minimum 3 spikes per pattern
pars.nsur = 10; % surrogates for core-cells, should be 1000 or more.
pars.prct = 99;% percentile on the surrogate core-cell distribution
pars.ccthr = 0.3; % minimal correlation between template and pattern
pars.sampfac = 1000; % sampling factor; controls how many subsamples ara drawn.
pars.thrmet = 'fit'; % method to automatically detect centroids

tic
[~,ensId,centInd,delta,rho,Ts,sampCent,templates,predbounds] = ...
    subSample_ensembles_pca(raster,pars.minspk,pars.sampfac,pars.dc,...
    pars.ishalo,pars.thrmet,pars.maxmem,pars.npcs,'test');

toc
%% Performance analysis: AUC vs number of samples
nsamps = length(sampCent);
%ns = 10:5:length(sampCent);
ns = 10:10:1000;
nreps = 24;
auc_vs_samps = zeros(length(ns),nreps);
mean_err = auc_vs_samps;
ens_num = auc_vs_samps;
for ss=1:length(ns)
    numsamps = ns(ss);
    tic
    parfor r=1:nreps
        sampid = randi(nsamps,[numsamps 1]);
        repCents = sampCent(sampid);
        repCents = cat(2,repCents{:});
        [~,pcs] = pca(repCents','NumComponents',pars.npcs);
        centcor = pdist2(pcs,pcs); %        
        [Nens,ensId,centInd,~,delta,rho,predbounds]=find_ens_by_density(centcor,pars.dc,pars.thrmet,1,0);
        
        templates = zeros(N,Nens);
        for i=1:Nens
            templates(:,i) = nanmean(repCents(:,centInd==i),2); % final centroids
        end
        
        % Template Matching
        [ens_seq,maxcor] = ens_raster_from_templates(templates,raster,pars.ccthr);
        
        % ensemble raster        
        ensmat_out = bsxfun(@eq,ens_seq,(1:Nens))';
        
        % re-assigning cluster to fit the input clusters
        C = 1-pdist2(ensmat_in,ensmat_out,'correlation');
        [~,ens_id] = max(C,[],2);
        % If any, here we look for the output ensembles that were not plugged into the raster
        if Nens>nens
            ee = 1:Nens;
            extraens = find(~ismember(ee,ens_id));
        else
            extraens =[];
        end
        ensmat_out = ensmat_out(ens_id,:);
        %ensmat_out = ensmat_out([ens_id;extraens(:)],:);
        
        % computing performance
        tpr_times = zeros(nens,1);
        fpr_times = tpr_times;
        for n=1:nens
            try
                [tpr_times(n),fpr_times(n)] = tpr_fpr_pks(ensmat_in(n,:),ensmat_out(n,:)); % computing tpr and fpr
            end
        end
        
        tpr_t_sort = [0; sort(tpr_times);1];
        fpr_t_sort = [0; sort(fpr_times);1];
        
        % AUC computation
        auc_vs_samps(ss,r) = trapz(fpr_t_sort,tpr_t_sort); %integrate . area under curve
        ens_num(ss,r) = Nens;
        % Computign mean error
        new_ens_seq_out = sum(bsxfun(@times, ensmat_out',1:nens),2);        
        new_ens_seq_in = sum(bsxfun(@times, ensmat_in',1:nens),2);
        mean_err(ss,r) = 1 - sum(new_ens_seq_out==new_ens_seq_in)/T;
        
    end
    toc
end



