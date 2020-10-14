function [Nens,ensId,centInd,haloInd,delta,rho,predbounds]=find_ens_by_density(dist,percNeigh,clusmet,ishalo,plt)
% [Nens,ensId,centInd,delta]=find_ens_by_density(dist,percNeigh,sdfact)
%
% Search for ensembles on the 'dist' distance matrix using the density
% peaks clustering methods by Rodriguez and Laio 2014.
%
% INPUTS
%
% 'dist' is a nbins x nbins matrix of distances between the significant bins
% found on data
% 'percNeigh' average percentage of neighbours, ranging from [0, 1]
% 'clusmet' is the method used for centroid search
% 'ishalo' is 1 if halo cluster refinement is used and 0 otherwise
% 'plt' is operator for plotting. If 1 plots, nothing otherwise.
%
% OUTPUTS
%
% 'Nens' is the number of ensembles found
% 'ensId' is the corresponding ensemble of each bin of dist
% 'centInd' is the moste representative bin of the ensemble
% 'delta' is the minimum distance between point it and any other point with
% higher density
% 'rho' is the local density of point i
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020


[~, rho] = paraSetv2(dist, percNeigh );
[Nens, ensId, centInd, haloInd, delta,predbounds] = densityClust_fast(dist, rho,ishalo,clusmet);

if plt~=0
    figure;
    plot(rho, delta, 's', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'b');hold on
    
    pcols = hsv(Nens);
    for i=1:Nens
        plot(rho(centInd==i), delta(centInd==i),'o', 'markersize',14,'color',pcols(i,:));hold on
    end
    title('Decision Graph');
    xlabel('\rho');
    ylabel('\delta');
    
end