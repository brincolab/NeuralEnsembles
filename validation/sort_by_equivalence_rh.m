function [ensmat_in_sort, enscells_in_sort, ens_id] = sort_by_equivalence_rh(ensmat_gt,ensmat_est,enscells_est,metric)
% Sort ensembles detected by its equivalences with generated ensembles sequences.
% INPUT:
%   ensmat_gt: <logical> nens_1-by-T matrix with the activity of each generated ensemble.
%   ensmat_est: <logical> nens_2-by-T matrix with the activity of each estimated ensemble.
%   enscells_est: <logical> N-by-nens_2 matrix with the core cells of each estimated ensemble.
%   metric: <string> with the metric used to find equivalences. See pdist2
%   for metrics
% OUTPUT:
%   ensmat_in_sort: <logical> nens-by-T-by-M matrix with the activity of each ensemble
%     detected sorted with respect to ensembles generated.
%   enscells_in_sort: <logical> N-by-nens_1 matrix with the core cells of each ensemble
%     detected sorted with respect to ensembles generated.

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

[nens1,T] = size(ensmat_gt);
[N,nens2] = size(enscells_est);
% min_nens = min(nens1,nens2);

% similarity matrix between GT ensemble sequence and estimated using
% 'metric'
seq_corr = 1-pdist2(ensmat_gt*1,ensmat_est*1,metric); 
% computing the maximum similarity and keeping the indexes of the
% detected ensembles that show maximal similarity with GT ensembles.

if nens1>nens2 % underdetection
    [~,ens_id] = max(seq_corr,[],1);
    ensmat_in_sort = false(nens1,T);
    enscells_in_sort = false(N,nens1);
    ensmat_in_sort(ens_id,:) = ensmat_est;
    enscells_in_sort(:,ens_id) = enscells_est;
    
else % equal or overdetection
    [~,ens_id] = max(seq_corr,[],2);
    ensmat_in_sort = ensmat_est(ens_id,:);
    enscells_in_sort = enscells_est(:,ens_id);
end

