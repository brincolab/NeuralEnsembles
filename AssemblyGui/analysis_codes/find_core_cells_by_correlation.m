function [neuronid,idthr,ens_cel_corr,sur_cel_cor] = find_core_cells_by_correlation(raster,ens_raster,nsur,p)
%
% [neuronid,idthr] = find_core_cells_by_correlation(raster,ens_raster,nsur,p)
%
% Detectes the core neurons of a given raster. A core neuron is defined as
% a neuron that whose correlation with an ensemble is bigger than the
% correlation by chance (random spike sequence with the same firing rate).
%
% INPUTS
%
% 'raster' is a N x T binary matrix, where raster(n,t)=1 if neuron 'n'
% fired at time 't' and 0 otherwise.
% 'ens_raster' is a nens x Tens binary matrix, where ens_raster(nen,t)=1 if
% ensemble nen is active at time t and 0 otherwise.
% 'nsur' is the number of artificial data to generate chance level
% 'p' is the percentile to use a threshold of the random distribution.
%
% OUTPUTS
%
% 'neuronid' is a Nx1 vector, where each entry is either 0 or 1, where 0 is
% non core neuron and 1 a core neuron.

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020
warning('off')
[N,T] = size(raster);
[nens,Tens] = size(ens_raster);
if T~=Tens && T~=N
    error('Raster and ensemble-raster length must be equal')    
end

ens_cel_corr=1-pdist2(raster,ens_raster,'correlation'); % correlation between cell and ensemble
sur_cel_cor = zeros(N,nens,nsur);
parfor s=1:nsur
    warning('off')
    sur_ens_seq = shuffle(ens_raster,'time');
    sur_cel_cor(:,:,s)=1-pdist2(raster,sur_ens_seq,'correlation'); % correlation between cell and ensemble
    
end

idthr = prctile(sur_cel_cor,p,3);
neuronid = ens_cel_corr>idthr;
