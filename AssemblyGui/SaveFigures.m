function SaveFigures(axes_handles)
% Saves the figure on the 'axes_handles' axes and opens a dialog box for
% choosing the place to save
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

[FileName,PathName] = uiputfile('*.png', 'Choose Place and Name to save figure');


fig_names = {'PCA','delta_rho','raster','ens_seq','ens_cel_corr','core_cells','inner_corr','npcs};

for i=1:length(axes_handles)
    f2=figure('paperunits','normalized', 'paperposition',[0.2 0.2 1 0.7],'visible', 'off');
    hax_new = copyobj(axes_handles(i), f2);
    set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
    if i==5 || i==6
        colormap(hax_new,redblue)
    end
    print(f2,'-dpng',[PathName,FileName(1:end-4),'_',fig_names{i},'.png'],'-r300');
    close(f2); %and close it
end


