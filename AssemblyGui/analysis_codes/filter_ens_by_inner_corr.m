function [ens_corr,corr_thr,selection] = filter_ens_by_inner_corr(raster,core_cells,sdfact)

C_cells = corrcoef(raster');
C_cells = triu(C_cells,1);
C_cells(C_cells==0)=nan;
noncore_corr = nanmean(nanmean(C_cells(sum(core_cells,2)==0,sum(core_cells,2)==0)));
noncore_std = nanstd(nanstd(C_cells(sum(core_cells,2)==0,sum(core_cells,2)==0)));

if isnan(noncore_corr)
    noncore_corr = nanmean(nanmean(C_cells));
    noncore_std = nanstd(nanstd(C_cells));
end
ens_corr = arrayfun(@(x) nanmean(nanmean(C_cells(core_cells(:,x),core_cells(:,x)))), 1:size(core_cells,2));
corr_thr = noncore_corr+sdfact*noncore_std;
selection = (ens_corr>corr_thr);