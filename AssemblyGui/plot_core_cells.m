function plot_core_cells(ax,core_cells,clims)
cla(ax);
hold(ax,'on');

[N,nens] = size(core_cells);

axes(ax)

% images
imagesc(core_cells, clims);
colormap(ax,'redblue')
xlim(ax,[0 nens+1])
ylim(ax,[0 N+1])

% lines
init_x = 0.5;
for e=1:nens+1
    plot(ax,[init_x init_x], [0 N+1], 'k-')
    hold(ax,'on');
    init_x = init_x + 1;
end
set(ax,'yticklabel','','ytick',[])
ylabel(ax,'')
box(ax,'on')

hold(ax,'off');