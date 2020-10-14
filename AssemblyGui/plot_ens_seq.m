function plot_ens_seq(ax,tt,ens_labs,ens_cols,sellabs)

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

nens = max(ens_labs);
cla(ax);
try
    if nens==0
        plot(ax,[tt(1) tt(end)],[0 0],'k-');
        return
    else
        hold(ax,'on');
        
        plot(ax,tt,ens_labs,'k--');
        hold(ax,'on');
        
        plot(ax,tt(ens_labs==0),(ens_labs(ens_labs==0)),'ko');
        hold(ax,'on');
        for e=1:nens
            plot(ax,tt(ens_labs==e),ens_labs(ens_labs==e),'.','markersize',25,'color',ens_cols(e,:));
            hold(ax,'on');
        end
        ylim(ax,[0 nens*1.1])
        
    end
catch
    plot(ax,[tt(1) tt(end)],[0 0],'k-');
end
%xlim(ax,[tt(1) tt(end)]);
xlabel(ax,'Time (s)')
ylabel(ax,'Ensemble')
set(ax,'ytick',1:nens,'yticklabel',sellabs);
hold(ax,'off');