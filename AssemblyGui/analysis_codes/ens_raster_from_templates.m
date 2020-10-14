function [ens_seq,maxcor,tmpcor] = ens_raster_from_templates(tmplts,raster,ccthr)
% ens_raster = ens_raster_from_templates(tmplts,raster,ccthr)
% 
% Generates an ensemble raster by template mastching of tmplts on the
% raster. If no template shows correlation higher than ccthr, then no
% ensemble is defined on that bin.
%
% INPUTS
% tmplts is a N x ntmplts matrix, where N is the number of neurons and
% ntmplts the number of templates
% raster is a N x T raster, where raster(n,t)=1 if neuron n fired at time T
% and 0 otherwise.
% ccthr is the threshold for correlation between template and a given bin
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020
if islogical(raster)
    raster = single(raster);
end
warning('off');
ntmplts = size(tmplts,2);
tmpcor = 1-pdist2(raster',tmplts','correlation'); % template matching
tmpcor(isnan(tmpcor))=0; % setting to 0 al nans
[maxcor,ens_seq] = max(tmpcor,[],2); % getting maximum correlation for each bin and corresponding index
ens_seq(maxcor<ccthr)=0; % removing any correlations lower than 0.1
%ens_raster = bsxfun(@eq,ens_seq,1:ntmplts)';

