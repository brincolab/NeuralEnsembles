%% clustering indexes on retinal data
sampCents = arrayfun(@(i) load(['v3_samp_500_halo_corrected_exp_',num2str(i),'_sampCent.mat']),1:4,'uni',0);
%sampCents_mat = cellfun(@(x) cat(2,x.sampCent{:}),sampCents,'uni',0);
sampCents_mat = cellfun(@(x) cat(2,x.sampCent{randi(length(x.sampCent),6000,1)}),sampCents,'uni',0);
sampCents_mat{4} = cat(2,sampCents{4}.sampCent{randi(length(sampCents{4}.sampCent),20000,1)});
sampCents_mat{2} = cat(2,sampCents{2}.sampCent{randi(length(sampCents{2}.sampCent),90000,1)});
%% re-clustering using differente npcs and dcs
npcs = 3:10;
dc = [0.05 0.01:0.01:0.1];
nnpcs = length(npcs);
ndcs = length(dc);
scores_exps=cell(4,2);
%%
for e=1:4
    clus_scores = zeros(11,nnpcs,ndcs);
    nclus = zeros(nnpcs,ndcs);
    [~,pcs] = pca(sampCents_mat{e}','NumComponents',npcs(end));
    for p=1:nnpcs
        centcor = pdist2(pcs(:,1:npcs(p)),pcs(:,1:npcs(p)));%        
        for d=1:ndcs
            [nclus(p,d),ensId,centInd,haloInd,delta,rho,predbounds]=find_ens_by_density(centcor,dc(d),'fit',1,0);
            
            cents = find(centInd);
            clus_scores(:,p,d) = clus_scores_Dunns_DB(centcor,ensId,cents);
        end
    end
    scores_exps{e,1}=clus_scores;
    scores_exps{e,2}=nclus;
    save('cluster_scores_checkpoint.mat','scores_exps')
end
%%
figure
e=2;
for p=1:nnpcs
    subplot(3,3,p)
    %plot(dc,squeeze(clus_scores(:,p,:)))
    plot(dc,squeeze(scores_exps{e,1}(1,p,:)))
    %semilogy(dc,squeeze(clus_scores(:,p,:)))
    %plot(dc,mean(squeeze((clus_scores(1:9,p,:)))))
end
%%
e=1;
nonan_samps = sampCents_mat{e};
nonan_samps(:,sum(isnan(nonan_samps))>0) =[];
[~,pcs,latent] = pca(nonan_samps','NumComponents',10);
%%
centcor = pdist2(pcs(:,1:6),pcs(:,1:6));%
[Nens,ensId,centInd,haloInd,delta,rho,predbounds]=find_ens_by_density(centcor,0.01,'fit',1,1);
%%
%[Nens, centInd,predbounds] = decisionGraph(rho, delta, 'fit');
%% Merging cluster 2 with 10 for exp 3
%ensId(ensId==2)=10;
%centInd(centInd==2)=10;
haloInd(haloInd==2)=10;
%%
centInd2= zeros(size(centInd));
cid = find(centInd>0);
centInd2(cid)=1:Nens-2;
dist2cent = centcor(cid,:); % re compute distances from centroid to any other point
[~,clustInd] = min(dist2cent); % re generates cluster labels  
%%
ensId = clustInd;
centInd = centInd2;
%%
figure
plot3(pcs(:,1),pcs(:,2),pcs(:,3),'.','color',[0.5 0.5 0.5],'markersize',5);hold on
%%
cols = othercolor('Dark28',Nens*2);
%selid = randi(length(haloInd),1000);
%selpcs = pcs(selid,:);
%selhalo = haloInd(selid);

figure
subplot(1,2,1)

plot3(pcs(haloInd==0,1),pcs(haloInd==0,2),pcs(haloInd==0,3),'.','color',[0.5 0.5 0.5],'markersize',1);hold on
%plot3(selpcs(selhalo==0,1),selpcs(selhalo==0,2),selpcs(selhalo==0,3),'.','color',[0.5 0.5 0.5],'markersize',2);hold on
cont = 1;
for c=1:Nens
    %plot3(selpcs(selhalo==c,1),selpcs(selhalo==c,2),selpcs(selhalo==c,3),'.','color',cols(c,:))
    plot3(pcs(haloInd==c,1),pcs(haloInd==c,2),pcs(haloInd==c,3),'.','color',cols(cont,:))    
    %plot3(pcs(ensId==c,1),pcs(ensId==c,2),pcs(ensId==c,3),'.','color',cols(cont,:))
    text(pcs(centInd==c,1),pcs(centInd==c,2),pcs(centInd==c,3),...
        num2str(c),'color','k','fontsize',16)
    hold on
    cont = cont+2;
end

% subplot(1,2,2)
% %plot3(pcs(haloInd==0,1),pcs(haloInd==0,2),pcs(haloInd==0,3),'.','color',[0.5 0.5 0.5],'markersize',1);hold on
% %plot3(selpcs(selhalo==0,1),selpcs(selhalo==0,2),selpcs(selhalo==0,3),'.','color',[0.5 0.5 0.5],'markersize',2);hold on
% cont = 1;
% for c=1:Nens
%     %plot3(selpcs(selhalo==c,1),selpcs(selhalo==c,2),selpcs(selhalo==c,3),'.','color',cols(c,:))
%     %plot3(pcs(haloInd==c,1),pcs(haloInd==c,2),pcs(haloInd==c,3),'.','color',cols(c,:))    
%     plot3(pcs(ensId==c,2),pcs(ensId==c,3),pcs(ensId==c,4),'.','color',cols(cont,:))
%     hold on
%     cont = cont+2;
% end

%
subplot(1,2,2)
%plot(rho(selid),delta(selid),'k.');hold on
plot(rho(centInd==0),delta(centInd==0),'k.');hold on
cont = 1;
for c=1:Nens
    plot(rho(centInd==c),delta(centInd==c),'.','color',cols(cont,:),'markersize',20);hold on
    text(rho(centInd==c),delta(centInd==c),num2str(c));hold on
    cont = cont+2;
end

%% extracting cluster templates as average
% extracting templates for each ensemble
enslist = unique(haloInd);
enslist(1)=[];
templates = zeros(size(nonan_samps,1),length(enslist));
for i=1:length(enslist)
    templates(:,i) = nanmean(nonan_samps(:,haloInd==enslist(i)),2); % final centroids
end
%%
figure
for n=1:length(enslist)
    subplot(5,5,n)
    imagesc(nonan_samps(:,haloInd==enslist(n)))
end
subplot(5,5,Nens+1)
imagesc(templates)
%% NOTAS
% experimento 1 usa 4 componentes principales y 0.01 con umbral de 99.9
% experimento 2 usa 4 comps y 0.01 con umbral de 99.99
% exsperimento 3 usa 5 comps y dc = 0.006 y el cluster 2 es re asignado al
% clñuster 10
% experimento 4 usa 6 comps y 0.01 y 99.999