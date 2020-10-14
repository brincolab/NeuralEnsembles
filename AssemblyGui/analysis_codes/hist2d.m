function [x_hist,y_hist,thr] = hist2d(xdata,ydata,nbins)
%%
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020
remids = isnan(xdata) | isnan(ydata);
xdata(remids)=[];
ydata(remids)=[];
[sortx,idx] = sort(xdata);
sorty = ydata(idx);
%bins = linspace(min(xdata),max(xdata),nbins);
bins = logspace(log(min(xdata)),log(max(xdata)),nbins);
[~,ids] = histc(sortx,bins);

x_hist=arrayfun(@(x) median(sortx(ids==x)),1:nbins);
y_hist=arrayfun(@(x) median(sorty(ids==x)),1:nbins);
x_lim=arrayfun(@(x) mad(sortx(ids==x),1),1:nbins);
y_lim=arrayfun(@(x) mad(sorty(ids==x),1),1:nbins);
xthr = x_hist+3*x_lim;
ythr = y_hist+3*y_lim;
remids = isnan(xthr) | isnan(ythr);
xthr(remids)=[];
ythr(remids)=[];

thr=interp1(xthr,ythr,xdata);
thr(isnan(thr))=0;

