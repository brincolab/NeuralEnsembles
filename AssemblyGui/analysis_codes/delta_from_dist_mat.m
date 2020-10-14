function delta = delta_from_dist_mat(dist, rho)
% computed delta for the rho vs delta plot. 
% based on clustering by density peaks
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020
[~, ordRho] = sort(rho, 'descend');

% matrix with rho in decreasing order. Column gtmat(i,j)=1 if rho(orhRho(i))>=rho(ordhRho(j)) and 0 otherwise
% the above matrix is anti-simmetric, so in the upper diagonal we have the  lower and in the lower the greater or equal than
gtmat = bsxfun(@ge,rho(ordRho),rho(ordRho)'); 
seldist = gtmat.*dist(ordRho,ordRho);% setting to 0 all the distance of lower density points
seldist = tril(seldist,-1); % lower diagonal is greater density
seldist(seldist==0)=inf; % to avoid get zero as minimal value
seldist(ordRho,ordRho)=seldist; % sorting to original order
[delta] = min(seldist,[],2); % taking the minimal distance to any other point with higehr density

clear seldist

%delta(ordRho(1)) = max(delta(~isinf(delta)));
delta(rho==max(rho)) = max(delta(~isinf(delta)));
delta(isinf(delta))=0;
delta = delta';