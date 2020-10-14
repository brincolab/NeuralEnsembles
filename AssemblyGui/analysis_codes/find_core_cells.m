function [neuronid,idthr] = find_core_cells(raster,nsur,p)
%
% [neuronid] = core_cells(raster)
%
% Detectes the core neurons of a given raster. A core neuron is defined as
% a neuron that participates on the raster more than expected by chance.
%
% INPUTS
%
% 'raster' is a N x T binary matrix, where raster(n,t)=1 if neuron 'n'
% fired at time 't' and 0 otherwise.
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

[N,T] = size(raster);
partidx = zeros(N,nsur);
PI = mean(raster,2);
parfor s=1:nsur
    surras = false(size(raster));
    [~,randis] = sort(rand(size(raster))); % random neuron labels
    colids = repmat((1:T),[N 1]);
    linids = sub2ind([N,T],randis(:),colids(:));
    surras(linids) = raster;     
    partidx(:,s) = mean(surras,2);

end
idthr = prctile(partidx,p,2);
neuronid = PI>=idthr;