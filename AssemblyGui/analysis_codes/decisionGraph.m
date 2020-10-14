function [numClust, centInd,predbounds] = decisionGraph(rho, delta, met)
% DECISIONGRAPH Decision graph for choosing the cluster centroids.
%   INPUT:
%       rho: local density [row vector]
%       delta: minimum distance between each point and any other point with higher density [row vector]
%       isManualSelect: 1 denote that all the cluster centroids are selected manually, otherwise 0
%  OUTPUT:
%       numClust: number of clusters
%       centInd:  centroid index vector
% Modified by RUBEN HERZOG 2017 in order to automatically detect the
% centroids based on diferent methods.
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

NE = length(rho);
numClust = 0;
centInd = zeros(1, NE);
predbounds = [];

switch met
    case 'manual'
        fprintf('Manually select a proper rectangle to determine all the cluster centres (use Decision Graph)!\n');
        fprintf('The only points of relatively high *rho* and high  *delta* are the cluster centers!\n');
        plot(rho, delta, 's', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'b');
        title('Decision Graph', 'FontSize', 17);
        xlabel('\rho');
        ylabel('\delta');
        
        rectangle = getrect;
        minRho = rectangle(1);
        minDelta = rectangle(2);
        
        for i = 1 : NE
            if (rho(i) > minRho) && (delta(i) > minDelta)
                numClust = numClust + 1;
                centInd(i) = numClust;
            end
        end
        
    case 'sur' % automaticaly detects centroids based on surrogate delta*rho
        nsur = 1000;
        surdot = zeros(nsur,NE);
        for i=1:nsur
            surdot(i,:) = delta(randperm(NE)).*rho(randperm(NE));
        end
        thr_sur = prctile(surdot(:),99); % p <0.01 % threshold on the product delta*rho
        thr_delta = mean(delta) + 1*std(delta); % threshold on delta values
        thr_rho = mean(rho) + 1*std(rho); % threshold on rho values
        centid = find((delta.*rho)>=thr_sur & delta>=thr_delta & rho>=thr_rho);
        numClust = length(centid);
        centInd(centid) = 1:numClust;
        
    case 'cutoff'; % uses cutoff of the delta*rho spectrum
        % Computes the derivative of the log of delta*rho and sets as the
        % number of clusters the maximum difference between consecutive
        % values.        
        [delta_rho,sort_idx] = sort(delta.*rho,'descend'); %sorts the product delta*rho
        maxlim = NE;
        % assumes that there are no more than 100 centroids
        if NE>100
            maxlim = 100;
        end
        [~,maxid] = max(abs(diff(log(delta_rho(1:maxlim))))); % absolute difference of the log of the product delta*rho        
        centid = sort_idx(1:maxid);
        numClust = length(centid);
        centInd(centid) = 1:numClust;
        
    case 'fit'
        % fitting a line to log(delta) vs log(rho) and using 99.9
        % upper confidence interval of fit as threshold        
        
        % Set up fittype and options.        
        ft = fittype( 'poly1' );        
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Robust = 'bisquare';
        mindelta = 10^-4; % min delta to be considers, improves the fit
        nzind = find(delta>mindelta & rho>0); % uses only points with delta>mindelta and rho>0
        nzdelta = double(delta(nzind));% just to make sure they are double
        nzrho = double(rho(nzind));
        
        % fitting in log space        
        [fitresult] = fit(log(nzrho)',log(nzdelta)', ft, opts ); 
        predbounds = predint(fitresult,log(nzrho),0.9999,'observation','on'); % prediction bound computation
        
        % selecting as centroids points above the confidence bound,
        % whose rhos values are at least higher than the 10-th perncentile
        % of rho distribution and points whose delta value are higher than
        % the 70-th percentile of the delta distribution
        %auxid = (nzdelta>exp(predbounds(:,2))' & nzrho>prctile(nzrho,10) & nzdelta>prctile(nzdelta,90));
        auxid = (nzdelta>exp(predbounds(:,2))');
        predbounds = cat(2,nzrho',exp(predbounds(:,2)));
        centid = nzind(auxid); % indices on the original basis
        numClust = length(centid);
        centInd(centid) = 1:numClust;
        
%         % %%
%         plot(rho,delta,'bo');hold on;
%         plot(rho(centid),delta(centid),'ks');
%         hold on;plot(nzrho,exp(predbounds(:,2)),'r.')
    case 'hist'
        % histogram method of rho and delta
        nbins = 10; % number of bins for the histogram
        nzind = find(delta>0 | rho>0); % uses only non zero delta and rho dataoints
        nzdelta = double(delta(nzind));
        nzrho = double(rho(nzind));
        [~,~,thr] = hist2d(nzrho,nzdelta,nbins);
        auxid = nzdelta>thr & nzrho>prctile(nzrho,10) & nzdelta>prctile(nzdelta,50);
        centid = nzind(auxid); % indices on the original basis                   
        numClust = length(centid);
        centInd(centid) = 1:numClust;
        predbounds = cat(2,nzrho',thr');
%         plot(rho,delta,'o');hold on;
%         plot(rho(centid),delta(centid),'ks');hold on        
%         plot(nzrho,thr,'.');
        
        
end

