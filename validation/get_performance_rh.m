% function performance = get_performance_rh(ensmat_gt,enscells_gt,ensmat_det,enscells_det)
function perf = get_performance_rh(ensmat_gt,enscells_gt,ensmat_det,enscells_det)
% Computes the performance of detection algorithm based on ROC AUC and
% Jaccard index.
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

nens = size(enscells_gt,2);
ensmat_global_gt = logical(sum(ensmat_gt,1));
ensmat_global_det = logical(sum(ensmat_det,1));
isens = sum(ensmat_det,2)>0; %Vector indicating which ensembles were detected

perf.fpr_indiv = nan(nens,1);
perf.tpr_indiv = perf.fpr_indiv;
perf.fpr_core = perf.fpr_indiv;
perf.tpr_core = perf.fpr_indiv;
perf.ham_indiv = perf.fpr_indiv;
perf.ham_core = perf.fpr_indiv;
perf.jac_indiv = perf.fpr_indiv;
perf.jac_core = perf.fpr_indiv;
perf.cor_indiv = perf.fpr_indiv;
perf.cor_core = perf.fpr_indiv;

% Calculate TPR-FPR for global activation sequence
[perf.fpr_global, perf.tpr_global] = fpr_tpr(ensmat_global_gt', ensmat_global_det');
% Calculate TPR-FPR for individual activation sequence
[perf.fpr_indiv(isens), perf.tpr_indiv(isens)] = fpr_tpr(ensmat_gt(isens,:)', ensmat_det(isens,:)');
% Calculate TPR-FPR for core-cells
[perf.fpr_core(isens), perf.tpr_core(isens)] = fpr_tpr(enscells_gt(:,isens), enscells_det(:,isens));
% Calculate AUC ROC
perf.auc_gs = trapz([0;sort(perf.fpr_global);1],[0;sort(perf.tpr_global);1]);
perf.auc_is = trapz([0;sort(perf.fpr_indiv);1],[0;sort(perf.tpr_indiv);1]);
perf.auc_cc = trapz([0;sort(perf.fpr_core);1],[0;sort(perf.tpr_core);1]);

% Calculate hamming distance
perf.ham_glob = 1-pdist2(ensmat_global_gt,ensmat_global_det,'ham');
perf.ham_indiv(isens) = diag(1-pdist2(ensmat_gt(isens,:),ensmat_det(isens,:),'ham'));
perf.ham_core(isens) = diag(1-pdist2(enscells_gt(:,isens)', enscells_det(:,isens)','ham'));

% Calculate jaccard distance
perf.jac_glob = 1-pdist2(ensmat_global_gt,ensmat_global_det,'jac');
perf.jac_indiv(isens) = diag(1-pdist2(ensmat_gt(isens,:),ensmat_det(isens,:),'jac'));
perf.jac_core(isens) = diag(1-pdist2(enscells_gt(:,isens)', enscells_det(:,isens)','jac'));

% Calculate correlation distance
perf.cor_glob = 1-pdist2(ensmat_global_gt,ensmat_global_det,'cor');
perf.cor_indiv(isens) = diag(1-pdist2(ensmat_gt(isens,:),ensmat_det(isens,:),'cor'));
perf.cor_core(isens) = diag(1-pdist2(enscells_gt(:,isens)', enscells_det(:,isens)','cor'));
