function [dc, rho] = paraSetv2(dist, percNeigh)
% using the modified version of Yger et al 2016, related to spike sorting.

[NE, ~] = size(dist);
dc = round(percNeigh.*NE); % percentage of nearest neighbors, S in the cited paper.
dist_sorted = sort(dist);
rho = 1./mean(dist_sorted(2:dc+1,:));
rho(isnan(rho))=0;

end
