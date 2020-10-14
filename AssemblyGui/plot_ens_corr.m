function plot_ens_corr(ax, ens_corr, corr_thr, ens_cols)

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020
cla(ax);
hold(ax,'on');

nens = length(ens_corr);

for e=1:nens
    bar(ax,e,ens_corr(e),'facecolor',ens_cols(e,:),'edgecolor',ens_cols(e,:));
    hold(ax,'on');
end

plot(ax,[0.5 nens+0.5],[corr_thr corr_thr],'r--');
hold(ax,'on');

xlim(ax,[0.5 nens+0.5])
ylim(ax,[0 max([ens_corr(:);corr_thr])*1.1])
xlabel(ax,'Ensemble Id.')
ylabel(ax,'Core-Cells Mean Correlation')

hold(ax,'off');