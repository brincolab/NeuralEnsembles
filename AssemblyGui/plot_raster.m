function plot_raster(ax,tt,raster,ens_labs,ens_cols)
% plots raster and colors each bin according to the identity given on
% ens_labs vector. ens cols is a nens x 3 matrix with colors for each
% ensemble. If ens_labs and ens_cols are empty, plots only the raster.
% ax is the axes to plot

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020
cla(ax);
[N,T] = size(raster);

try
    if isempty(ens_labs) || isempty(ens_cols) % plots raster with no color
        for n=1:N
            plot(ax,tt(raster(n,:)==1),raster(n,raster(n,:)==1)*n,'k.');
            hold(ax,'on')
        end
    else
        nens = max(ens_labs);
        for n=1:N
            plot(ax,tt(ens_labs==0 & raster(n,:)==1),raster(n,ens_labs==0 & raster(n,:)==1)*n,'k.'); % plotting bins with no ensemble activity
            hold(ax,'on')
            
            for e=1:nens
                plot(ax,tt(ens_labs==e),raster(n,ens_labs==e)*n,'.','color',ens_cols(e,:)); % plotting bins WITH ensemble activity
                hold(ax,'on')
            end
        end
    end
catch
    plot(ax,[0 tt(end)],[0 0],'k-');
end

%xlim(ax,[0 tt(end)])
ylim(ax,[0 N+1])
ylabel(ax,'Neuron')
set(ax,'xticklabel','')
hold(ax,'off');

