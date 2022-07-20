% exp 1 uses 6 pcs and 0.035 on the global centroid clustering
% exp 2 uses 6 pcs and 0.022 on the global centroid clustering
% exp 3 uses 5 pcs and 0.015 on the global centroid clustering and manual
% exp 4 uses 5 pcs and 0.015 on the global centroid clustering and manual
% selection
%% Testing centroid detecting with sub samples
%load('v3_samp_500_halo_corrected_exp_1_sampCent.mat');
%load('v3_samp_500_halo_corrected_exp_4_sampCent.mat');
load('v4_samp_1000_halo_corrected_exp_4_sampCent.mat');
samps = cat(2,sampCent{:});
%load('v3_samp_500_halo_corrected_exp_4_sampCent.mat');
%samps = cat(2,samps,sampCent{:});
%% 1.- Esitmating PCA for all dataset
[coeff,pcs,latent] = pca(samps','numcomponents',6);

%% 2.- Taking random subsample and computing centroids
nsamps = length(samps);
samplist = 1:nsamps;
samplistidx = samplist;
nreps = 36;
centcoord = cell(nreps,1);
nsubsamp = floor(nsamps/nreps);
% Making indexes to sample without repetition
randidx = randperm(nsamps);
randidx = reshape(randidx(1:(nreps*nsubsamp)),[nsubsamp nreps]);
%%    
parfor r=1:nreps    
    centcor = pdist2(pcs(randidx(:,r),:),pcs(randidx(:,r),:));
    [Nens,ensId,centInd,haloInd,delta,rho,predbounds]=find_ens_by_density(centcor,0.01,'fit',1,0);
    auxcentscoord = pcs(randidx(centInd>0,r),:);
    centcoord{r} = auxcentscoord;    
end
%% Checking centroids distance
globalcents = cat(1,centcoord{:});
centdist = pdist2(globalcents,globalcents);
%
plot3(globalcents(:,1),globalcents(:,2),globalcents(:,3),'.')
%%
[Nens,ensId,centInd,haloInd,delta,rho,predbounds]=find_ens_by_density(centdist,0.015,'manual',1,1);
%
cols = othercolor('Dark28',Nens*2);
% subplot(1,3,1)
% aa = cumsum(latent)*100;
% loglog(aa,'o-');hold on
% loglog(6,aa(6),'rs')
% ylim([0 100])
% ylabel('% of Variance')
% xlabel('Rank')
%%
subplot(1,2,1)
cont = 1;
for c=1:Nens
    plot3(globalcents(ensId==c,1),globalcents(ensId==c,2),globalcents(ensId==c,3),'.','color',cols(cont,:))    
    hold on;
    text(globalcents(centInd==c,1),globalcents(centInd==c,2),globalcents(centInd==c,3),...
    num2str(c),'color','k','fontsize',16)
    cont = cont+2;
end
subplot(1,2,2)
plot(rho(centInd==0),delta(centInd==0),'k.');hold on
cont = 1;
for c=1:Nens
    plot(rho(centInd==c),delta(centInd==c),'.','color',cols(cont,:),'markersize',20);hold on
    text(rho(centInd==c),delta(centInd==c),num2str(c));hold on
    cont = cont+2;
end

%% Computing centroid average coordinates
ave_cent = arrayfun(@(x) mean(globalcents(ensId==x,:)),1:Nens,'uni',0);
ave_cent = cat(1,ave_cent{:});
% Computing distance from all datapoints to centroids
dist2cent = pdist2(ave_cent,pcs);
[~,clustInd] = min(dist2cent);
[haloInd] = haloAssignv2(dist2cent, clustInd);
%%
cols = othercolor('Dark28',Nens*2);
selid = randperm(length(pcs),20000);
selhalo = haloInd(selid);
selpcs = pcs(selid,:);
figure('units','normalized','outerposition',[0 0 1 1],'PaperPositionMode', 'auto','visible','on')
subplot(1,2,1)
cont = 1;
%plot3(pcs(haloInd==0,1),pcs(haloInd==0,2),pcs(haloInd==0,3),'.','color',[0.5 0.5 0.5],'markersize',1);
hold on
for c=1:Nens
    %plot3(pcs(clustInd==c,1),pcs(clustInd==c,2),pcs(clustInd==c,3),'.','color',cols(cont,:))    
    %plot3(pcs(haloInd==c,1),pcs(haloInd==c,2),pcs(haloInd==c,3),'.','color',cols(cont,:))    
    plot3(selpcs(selhalo==c,1),selpcs(selhalo==c,2),selpcs(selhalo==c,3),'.','color',cols(cont,:))    
    hold on;
    %text(ave_cent(c,1),ave_cent(c,2),ave_cent(c,3)*1.5,...
    %num2str(c),'color','k','fontsize',14)
    cont = cont+2;
end
grid;
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')

subplot(1,2,2)
plot(rho(centInd==0),delta(centInd==0),'k.');hold on
%plot(predbounds(:,1),predbounds(:,2),'b.');hold on
cont = 1;
for c=1:Nens
    plot(rho(centInd==c),delta(centInd==c),'.','color',cols(cont,:),'markersize',30);hold on
    text(rho(centInd==c)-0.1,delta(centInd==c)+0.05,num2str(c),'fontsize',12);hold on
    cont = cont+2;
end
ylabel('\delta')
xlabel('\rho')
ylim([0 1.2])
xlim([0 150])
grid;
%% Extracting templates and plotting
% extracting templates for each ensemble
templates = zeros(size(samps,1),Nens);
for i=1:Nens
    templates(:,i) = nanmean(samps(:,haloInd==i),2); % final centroids
end
%%
figure
for n=1:Nens
    subplot(5,5,n)
    imagesc(samps(:,haloInd==n))
end
subplot(5,5,Nens+1)
imagesc(templates)
