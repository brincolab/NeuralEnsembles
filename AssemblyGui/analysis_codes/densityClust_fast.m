function [numClust, clustInd, centInd, haloInd, delta,predbounds] = densityClust_fast(dist, rho, isHalo, autosel)
% DENSITYCLUST Clustering by fast search and find of density peaks.
%   SEE the following paper published in *SCIENCE* for more details:
%       Alex Rodriguez & Alessandro Laio: Clustering by fast search and find of density peaks,
%       Science 344, 1492 (2014); DOI: 10.1126/science.1242072.
%   INPUT:
%       dist: [NE, NE] distance matrix
%       rho: local density [row vector] (previously calculated using a
%       given dc)
%       isHalo: 1 denotes that the haloInd assigment is provided, otherwise 0.
%       autosel: string with the method used for automcatic centroid
%       detection
%   OUTPUT:
%       numClust: number of clusters
%       clustInd: cluster index that each point belongs to, NOTE that -1 represents no clustering assignment (haloInd points)
%       centInd:  centroid index vector
%       haloInd: haloInd row vector [0 denotes no haloInd assignment]
%
%
% MODIFIED BY RUBEN HERZOG 2017.

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

[sortrho, ordRho] = sort(rho, 'descend');

% matrix with rho in decreasing order. Column gtmat(i,j)=1 if rho(orhRho(i))>=rho(ordhRho(j)) and 0 otherwise
% the above matrix is anti-simmetric, so in the upper diagonal we have the  lower and in the lower the greater or equal than
gtmat = bsxfun(@ge,sortrho,sortrho'); 
seldist = gtmat.*dist(ordRho,ordRho);% sel(tting to 0 all the distance of lower density points
seldist = tril(seldist,-1); % lower diagonal is greater density
seldist(seldist==0)=inf; % to avoid get zero as minimal value
%seldist(ordRho,ordRho)=seldist; % sorting to original order
[delta] = min(seldist,[],2); % taking the minimal distance to any other point with higehr density
delta(ordRho) = delta;

clear seldist

%delta(ordRho(1)) = max(delta(~isinf(delta)));
delta(rho==max(rho)) = max(delta(~isinf(delta))); % most dense points is assigned to maximal delta.
delta(isinf(delta))=0;
delta = delta';

[numClust, centInd,predbounds] = decisionGraph(rho, delta, autosel);

% each point is assigned to its closest centroid
if numClust==1
    clustInd = ones(length(delta),1);
else
    clusidx = find(centInd); % indices of centroids    
    dist2cent = dist(clusidx,:); % distance from centroid to any other point    
    [~,clustInd] = min(dist2cent);
    clustInd(clusidx) = centInd(clusidx);
    % setting to 0 clusters with 1 member    
    cluscounts = histc(clustInd,1:numClust); % number of members on each cluster
    one_mem_clus = find(cluscounts==1 | cluscounts==0); % finds indices of one or zero member clusters
    if ~isempty(one_mem_clus)
        numClust = numClust - length(one_mem_clus);
        cluslab = centInd(clusidx); % label of centroid
        clusidx(ismember(cluslab,one_mem_clus))=[]; % removing centroids with just one member
        centInd = zeros('like',centInd);
        centInd(clusidx)=1:numClust;
        dist2cent = dist(clusidx,:); % re compute distances from centroid to any other point        
        [~,clustInd] = min(dist2cent); % re generates cluster labels       
        clustInd(clusidx) = centInd(clusidx);
    end
    
end


%haloInd = haloAssign(dist, clustInd, dc, rho, isHalo);
try
    haloInd = haloAssign(dist, clustInd,centInd,isHalo);
catch
    haloInd = clustInd;
end

end


