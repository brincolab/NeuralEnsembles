%function [haloInd] = haloAssign(dist, clustInd, dc, rho, isHalo)
function [haloInd] = haloAssignv2(dist2cent, clustInd)
haloInd = clustInd; % 0 denotes no halo assignment
nens = max(clustInd);
sameclusmat = bsxfun(@eq, (1:nens)',clustInd); % matrix with points that belong to the same cluster
dist2cluscent = dist2cent.*sameclusmat; % makes zero all distances to other centroids
dist2cluscent(isnan(dist2cluscent))=0;

meandist2cent = bsxfun(@gdivide,sum(dist2cluscent,2),sum(sameclusmat,2)); % mean distance to centroids
lemeandist2cent = bsxfun(@le,dist2cent,meandist2cent); % lower or equal to the average distance to centroid
remids = ~sum(lemeandist2cent); % ids to remove
haloInd(remids) = 0;
end
