function [Nens,ensId,centInd,delta,rho,templates,predbounds,pcs,bincor] = ensembles_pca(raster,minSpk,prctdata,isHalo,thrMet,npcs,filename)
%
% [Nens,ensId,delta,rho,Ts,sampCent,templates] = subSample_ensembles(raster,minSpk,sampFac,kernel,prctdata,thrMet,savemem,isGpu)
%
% Find ensemble on several subsamples of the raster T. Defines the size of
% the subsample Ts based on the available memory and 'sampFac', which if >1
% the raster is oversampled (i.e. 'sampFac'*T/Ts). For each subsample,
% clusteriong is performed and a set of representative population patterns
% are extracted. Finally, all the represdentative paterns extracted over
% all sampled are pooled together and clustered again, yielding the set of
% ensemble present on all the recording.
%
% INPUTS
%
% 'raster' N x T binary matrixc whre raster(i,t)=1 if neuron i spikes at
% time t and 0 otherwise.
% 'minSpk' minimum number of population spikes to consider in analysis
% 'sampFac' oversampling factor to determine the number of subsamples.
% Multiplies T/Ts. If >1, oversampling, if <1, less samples than the whole
% recording are extracted, i.e. undersampling.
% 'prctdata' percentiel of the distance distribution to compute the density
% of the clustering step.[0,1] usually 0.01-0.02
% 'isHalo' is 1 if halo cluster refinement is used and 0 otherwise
% 'thrMet' if 0 manual selection of centroids is performed. If =>1,
% 'thrMet' shuffled version of delta*rho are generated to determine the
% threshold for the cluster centroids. If <0, uses the sorted values of
% delta*rho to find the cluster centroids.
% 'maxmem' max memory to use in bytes. Should be at least 1Gb lower than
% the total RAM available
% 'isGpu' =1 if gpu is presented and will be used and 0 otherwise.
% 'npcs' is the number of principal components to use
%
% OUTPUTS
% 'Nens' is the number of cluster,i.e. ensembles, found on all raster
% 'ensId' is a vector, where ensId(t)=E, if ensemble E was active on t entry
% 'delta' are the delta values for clustering
% 'rho' the density values
% 'Ts' the length of the subsample
% 'templates' is a N x Nens matrix, with the representative template of the
% each ensemble.
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

N = size(raster,1);
raster = raster(:,sum(raster)>minSpk);
[~,pcs] = pca(raster','NumComponents',npcs); % pca with npcs num components
bincor = pdist2(pcs,pcs); % euclidean distance on principal component space
[Nens,~,centInd,ensId,delta,rho,predbounds]=find_ens_by_density(bincor,prctdata,thrMet,isHalo,0);
    

% extracting templates for each ensemble
templates = zeros(N,Nens);
for i=1:Nens
    templates(:,i) = nanmean(raster(:,ensId==i),2); % final centroids
end

% saving
%save([filename,'_dc_',num2str(prctdata),'_npcs_',num2str(npcs),'_ensembles_output.mat'],'ensId','centInd','templates');
