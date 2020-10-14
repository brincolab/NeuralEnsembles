function plot_delta_rho(ax,rho,delta,cents,predbounds,ens_cols)
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

[~,pred_id] = sort(predbounds(:,1));

%pred bound
cla(ax);
plot(ax,predbounds(pred_id,1),predbounds(pred_id,2),'b--','linewidth',2);
hold(ax,'on');

% delta rho
% maxrho = 10;
% rho(rho>maxrho) =maxrho;
plot(ax,rho,delta,'k.');
hold(ax,'on');
nens = sum(cents>0);
for e=1:nens
    plot(ax,rho(cents==e),delta(cents==e),'.','markersize',25,'color',ens_cols(e,:));
    text(rho(cents==e)*1.01,delta(cents==e)*1.01,num2str(e), 'Parent', ax,'fontsize',12);
    hold(ax,'on');
end

xlim(ax,[0 max(rho(~isinf(rho))).*1.1])
ylim(ax,[0 max(delta(~isinf(delta)))*1.1]);
xlabel(ax,'\rho')
ylabel(ax,'\delta')
hold(ax,'off');