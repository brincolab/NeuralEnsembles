function plot_eigs(ax,eigs,seleig)

cla(ax);
hold(ax,'on');

plot(ax,1:length(eigs),eigs,'o--')
%cumeigs = cumsum(eigs);
%plot(ax,1:length(eigs),cumeigs,'o--')
plot(ax,seleig,eigs(seleig),'rs','markersize',20);
set(ax,'xscale','log','yscale','linear','xlim',[1 length(eigs)]);
xlabel(ax,'PCs')
ylabel(ax,'% var');

hold(ax,'off');