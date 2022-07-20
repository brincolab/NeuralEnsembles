%% Making the ensemble raster
load('/home/octodon/Documents/MATLAB/FastEnsDet/4exps/4_exps_6_conds_concatenated.mat')
nexps = 4;
templates =  arrayfun(@(x) load(['v4_templates_exp_',num2str(x),'.mat']),1:nexps,'uni',0);
templates = {templates{1}.templates,templates{2}.templates,templates{3}.templates,templates{4}.templates};
%%
pars.ccthr = 0.3;
[ens_seq,maxcor,tmpcor] = cellfun(@(x,y) ens_raster_from_templates(x,y,pars.ccthr),templates',cat_rasters,'uni',0);

%% checking correlation distribution
bins = 0.005:0.005:1;
figure
for e=1:nexps
    subplot(2,2,e)
    %histc(tmpcor{e}(:)
    Nens = size(tmpcor{e},2);
    for n=1:Nens
        [nel] = histcounts(tmpcor{e}(:,n),bins);
        plot(bins(2:end),nel);hold on
    end
end
%% checking the max correlation distribution
figure
for e=1:nexps
    subplot(2,2,e)
    %[nel] = histcounts(maxcor{e}(:,n),bins);
    %[nel,base] = histcounts(maxcor{e});
    %plot(base(2:end),nel);hold on
    [nel] = histcounts(maxcor{e},bins,'normalization','probability');
    plot(bins(2:end),nel);hold on
    plot([mean(maxcor{e}(maxcor{e}>0)) mean(maxcor{e}(maxcor{e}>0))],...
        [0 0.02],'r--');hold on
    plot([mean(maxcor{e}(maxcor{e}>0)) mean(maxcor{e}(maxcor{e}>0))]- 1.*std(maxcor{e}(maxcor{e}>0)),...
        [0 0.02],'r--');hold on
    ylim([0 0.02])
    
end
%%
nconds = 6;
condnames = {'iSA','iPA','WN','NM', 'fPA','fSA'};
for e=1%:4
    Nens = size(tmpcor{e},2);
    cols = othercolor('Dark28',Nens*2);    
    condend = cumsum(raster_lengths(e,:));
    condend = [1 condend];
    figure
    for c=1%:nconds
        cont = 1;
        %subplot(2,3,c,'align')
        
        for n=1:Nens
            try
            %xx = find(ens_seq{e}(condend(c+1)-2000:condend(c+1))==n);
            xx = find(ens_seq{e}(condend(c):condend(c+1))==n);
            plot(xx,n,'.','color',cols(cont,:));hold on
            end
            
            cont = cont+2;
        end
        ylim([0 Nens+1])
        title(condnames{c})
    end
end

%%
nconds = 6;
condnames = {'iSA','iPA','WN','NM', 'fPA','fSA'};
for e=1%:4
    Nens = size(tmpcor{e},2);
    cols = othercolor('Dark28',Nens*2);    
    condend = cumsum(raster_lengths(e,:));
    condend = [1 condend];
    figure
    for c=1%:nconds
        cont = 1;
        %subplot(2,3,c,'align')
        
        for n=1:Nens
            try
            xx = condend(c):condend(c+1);
            %xx = find(ens_seq{e}(condend(c):condend(c+1))==n);
            plot(xx,tmpcor{e}(xx,n)+n,'color',cols(cont,:));hold on
            end
            
            cont = cont+2;
        end
        ylim([-1 Nens+2])
        title(condnames{c})
    end
end