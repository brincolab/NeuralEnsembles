%% Methods parameters
pars.sr = 1; % Hz
pars.bin = 0.02; % seconds
pars.dc = 0.01; % cut-off for distances
pars.ishalo = 1; % not using halo refinement of the clusters
pars.npcs = 5;
pars.maxmem = 8*10^9; % maximal Gb of RAM 
pars.minspk = 3; % minimum 3 spikes per pattern
pars.nsur = 10; % surrogates for core-cells, should be 1000 or more.
pars.prct = 99;% percentile on the surrogate core-cell distribution
pars.ccthr = 0.1; % minimal correlation between template and pattern
pars.sampfac = 2; % sampling factor; controls how many subsamples ara drawn.
pars.thrmet = 'fit'; % method to automatically detect centroids

%% extracting esnembles from 4 experiments
load('/home/octodon/Documents/MATLAB/FastEnsDet/4exps/4_exps_6_conds_concatenated.mat')
nexps = 4;
ensData = cell(nexps,1);

%% tomando todo el registro
pars.maxmem = 10*10^9; % maximal Gb of RAM 
pars.sampfac = 1000; % sampling factor; controls how many subsamples ara drawn.
ens_out = struct();
%%
for e=1:nexps   
    tic
    [ens_out.sampCent,ens_out.Ts] = ...
    subSample_ensembles_pca(cat_rasters{e},pars.minspk,pars.sampfac,pars.dc,...
    pars.ishalo,pars.thrmet,pars.maxmem,pars.npcs,['v4_samp_1000_halo_corrected_exp_',num2str(e)]);

    ensData{e} = ens_out;
    save(['v4_checkpoint.mat'],'ensData','pars');
    toc
end

%% Plotting the variance for each set
figure
for e=1:4
    [~,~,latent] = pca(ensData{e}.sampCent');
    subplot(2,2,e)
    semilogx(latent(1:20))
end
    
    


%% re-clustering using differente npcs and dcs
npcs = 6;
prctdata = 0.01;
for e=4%:nexps
    figure
    [~,pcs] = pca(ensData{e}.sampCent','NumComponents',npcs);
    centcor = pdist2(pcs,pcs); %
    [Nens,ensId,centInd,haloInd,delta,rho,predbounds]=find_ens_by_density(centcor,prctdata,'fit',1,0);
    
    %%
    cols = othercolor('Dark28',Nens);
    selid = randi(length(haloInd),1000);
    selpcs = pcs(selid,:);
    selhalo = haloInd(selid);
    
    subplot(1,2,1)    
    %plot3(pcs(haloInd==0,1),pcs(haloInd==0,2),pcs(haloInd==0,3),'.','color',[0.8 0.8 0.8]);hold on
    plot3(selpcs(selhalo==0,1),selpcs(selhalo==0,2),selpcs(selhalo==0,3),'.','color',[0.5 0.5 0.5],'markersize',2);hold on
    for c=1:Nens
        plot3(selpcs(selhalo==c,1),selpcs(selhalo==c,2),selpcs(selhalo==c,3),'.','color',cols(c,:))
        hold on
    end
    
    subplot(1,2,2)    
    plot(rho(selid),delta(selid),'k.');hold on
    for c=1:Nens
        plot(rho(centInd==c),delta(centInd==c),'.','color',cols(c,:),'markersize',20);hold on
    end
    
end

%% Plotting results: delta vs rho
figure
for e=1:4
    subplot(2,2,e)
    plot(ensData{e}.rho,ensData{e}.delta,'.');hold on
    plot(ensData{e}.rho(ensData{e}.centInd>0),ensData{e}.delta(ensData{e}.centInd>0),'ko');hold on
end

%% plotting pcs
figure
for e=1:4
    [~,pcs] = pca(ensData{e}.sampCent');
    subplot(2,2,e)
    cols = hsv(ensData{e}.Nens);
    for c=1:ensData{e}.Nens
        plot3(pcs(ensData{e}.ensId==c,1),pcs(ensData{e}.ensId==c,2),pcs(ensData{e}.ensId==c,3),'.','color',cols(c,:))
        hold on
    end
    
end

%%
figure
for e=1:4
    [~,pcs] = pca(ensData{e}.sampCent');
    subplot(2,2,e)
    cols = hsv(ensData{e}.Nens);
    plot3(pcs(ensData{e}.haloInd==0,1),pcs(ensData{e}.haloInd==0,2),pcs(ensData{e}.haloInd==0,3),...
        '.','color',[0.5 0.5 0.5]);hold on
    for c=1:ensData{e}.Nens
        plot3(pcs(ensData{e}.haloInd==c,1),pcs(ensData{e}.haloInd==c,2),pcs(ensData{e}.haloInd==c,3),'.','color',cols(c,:))
        
        hold on
    end
    alpha(0.5)
end

