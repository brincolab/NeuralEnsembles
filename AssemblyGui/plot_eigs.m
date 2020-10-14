function plot_eigs(ax,eigs,seleig)

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

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