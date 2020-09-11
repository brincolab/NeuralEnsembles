function plot_pca(ax,pcs,ens_labs,ens_cols)
cla(ax);
np = size(pcs,2);
if isempty(ens_labs) || isempty(ens_cols) % plots raster with no color
    if np>2
        plot3(ax,pcs(:,1),pcs(:,2),pcs(:,3),'k.');
    else
        plot(ax,pcs(:,1),pcs(:,2),'k.');
    end
    
    hold(ax,'on')
    
else
    nens = max(ens_labs);
    if np>2
        plot3(ax,pcs(ens_labs==0,1),pcs(ens_labs==0,2),pcs(ens_labs==0,3),'k.');
        hold(ax,'on')
        for e=1:nens
            plot3(ax,pcs(ens_labs==e,1),pcs(ens_labs==e,2),pcs(ens_labs==e,3),'marker','.','color',ens_cols(e,:),'linestyle','none'); %
        end
    else
        plot(ax,pcs(ens_labs==0,1),pcs(ens_labs==0,2),'k.');
        hold(ax,'on')
        for e=1:nens
            plot(ax,pcs(ens_labs==e,1),pcs(ens_labs==e,2),'marker','.','color',ens_cols(e,:),'linestyle','none'); %
            
        end
    end
    hold(ax,'on')
end

xlabel(ax,'PC1')
ylabel(ax,'PC2')
zlabel(ax,'PC3')
hold(ax,'off');
