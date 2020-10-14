function [dc, rho] = paraSetv2(dist, percNeigh)
% using the modified version of Yger et al 2016, related to spike sorting.
%
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

[NE, ~] = size(dist);
dc = round(percNeigh.*NE); % percentage of nearest neighbors, S in the cited paper.
dist_sorted = sort(dist);
rho = 1./mean(dist_sorted(2:dc+1,:)); % Yger 2016
rho(isnan(rho))=0;

end
