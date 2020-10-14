%function [haloInd] = haloAssign(dist, clustInd, dc, rho, isHalo)
function [haloInd] = haloAssign(dist, clustInd,centInd,isHalo)
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020
haloInd = clustInd; % 0 denotes no halo assignment
if isHalo == 1
        
    % neqAndCloseRho = bsxfun(@times,and(~bsxfun(@eq, clustInd,clustInd'),triu(dist<dc,1)),rho);
    % the above matrix contains the respective rho values of the points
    % that belong to different clusters but closer than dc.
    
    % this apporach only removes from the clusters the points in between
    % both of them
    % function for the first approach function [haloInd] = haloAssign(dist, clustInd, numClust, dc, rho, isHalo)
    %gg=arrayfun(@(x) (max(neqAndCloseRho(:,clustInd==x),[],2)<=rho') ,1:numClust,'uni',0);
    %haloId = sum(cat(2,gg{:}),2)<numClust;
    
    % this approach removes all the points with density lower than the
    % average density on the border of the clusters. Thus, removes both the
    % points in between clusters and those surrounding them but far away
    % from the centroid.
    % maxBorderRho = max(neqAndCloseRho(:));
%     maxBorderRho = mean(neqAndCloseRho(neqAndCloseRho>0));
%     haloId = rho<=maxBorderRho;
%         
%     haloInd(haloId) = 0;
    
    
    % this approach removes points that are farther than the average
    % intracluster distance to the centroid    
    sameclusmat = bsxfun(@eq, clustInd,clustInd'); % matrix with points that belong to the same cluster
    sameclus_cent = sameclusmat(centInd>0,:); % selecting only centroid rows
    dist2cent = dist(centInd>0,:); % distance from all points to the centroids
    dist2cluscent = dist2cent.*sameclus_cent; % makes zero all distances to other centroids
    dist2cluscent(isnan(dist2cluscent))=0;
    
    meandist2cent = bsxfun(@gdivide,sum(dist2cluscent,2),sum(sameclus_cent,2)); % mean distance to centroids
    lemeandist2cent = bsxfun(@le,dist2cent,meandist2cent); % lower or equal to the average distance to centroid
    remids = ~sum(lemeandist2cent); % ids to remove
    haloInd(remids) = 0;
    
    
    
end



%
%% 2.- computes average rho between both points satifying above condition
%avgRho = (bsxfun(@plus,rho,rho')./2).*neqAndClose;
end
