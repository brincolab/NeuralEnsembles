%% Computing the ensemble sequences
load('/home/octodon/Documents/MATLAB/raster2Ens/Analyzing4Exp_FullRecording/v4_4_exps_ens_seq.mat')
load('/home/octodon/Documents/MATLAB/FastEnsDet/4exps/4_exps_6_conds_concatenated.mat')
nexps = 4;
condnames = {'iSA','iPA','WN','NM', 'fPA','fSA'};
%% removing 0s from the sequence
ens_seq_nz = cellfun(@(x) x(x>0),ens_seq,'uni',0);

%% ensemble sequence for stimuli
ras_ind = cumsum(cat(2,ones(4,1),raster_lengths),2);
ens_seq_stim = cell(4,6);
for e=1:4
    for c=1:6
        ens_seq_stim{e,c} = ens_seq{e}(ras_ind(e,c):ras_ind(e,c+1)-1);
    end
end
%% removing zero
ens_seq_nz_stim = cellfun(@(x) x(x>0),ens_seq_stim,'uni',0);

%% computing transition matrices
transmat_full = cellfun(@(x) transmat_from_seq(x),ens_seq_nz,'uni',0);
transmat_stim = cellfun(@(x) transmat_from_seq(x),ens_seq_nz_stim,'uni',0);
%% computing correlation between transmats
tm_cor = zeros(6,6,4);
for e=1:4
    aux = zeros(6);
    for c=1:6
        x = transmat_stim{e,c};
        x(isnan(x))=0;
        for c1=c+1:6            
            y = transmat_stim{e,c1};
            y(isnan(y))=0;
            C = corrcoef(x(:),y(:));
            aux(c,c1)=C(1,2);
        end
    end
    tm_cor(:,:,e) = aux;
end
%% images
figure
for e=1:4
    subplot(2,2,e)
    imagesc(tm_cor(:,:,e))
    colorbar
    axis square
end
            
        
%% max for each condition
max_transprob = max(max(cellfun(@(x) max(x(:)),transmat_stim),[],2),cellfun(@(x) max(x(:)),transmat_full));
min_transprob = min(min(cellfun(@(x) min(x(x>0)),transmat_stim),[],2),cellfun(@(x) min(x(x>0)),transmat_full));
%% plotting transitions matrices for stimuli and full
e=4;
for c=1:6
    subaxis(1,7,c)
    imagesc(log(transmat_stim{e,c}),log([min_transprob(e) max_transprob(e)]))
    axis square
    set(gca,'xticklabel','','yticklabel','');
end
subaxis(1,7,7)
imagesc(log(transmat_full{e}),log([min_transprob(e) max_transprob(e)]))
axis square
set(gca,'xticklabel','','yticklabel','');
%% surrogate transitions
nsur = 9600;
nblocs = 12;
blocksz = floor(nsur/nblocs);

sur_transmatall = cell(4,nblocs);
for e=1:4
    seq = ens_seq_nz{e};
    parfor n=1:nblocs
        [~,ids] = (sort(randi(2000,[length(seq) blocksz],'uint16')));
        repaux = repmat((seq),[1 blocksz]);
        repaux = repaux(ids); % shuffling each column
        repaux = num2cell(repaux,1);
        sur_transmat = cellfun(@transmat_from_seq,repaux,'uni',0);
        sur_transmat = cat(3,sur_transmat{:});
        sur_transmatall{e,n} = sur_transmat;
    end
end

% conds
sur_transmatall_conds = cell(4,6,nblocs);
for e=1:4
    for c=1:6
        seq = ens_seq_nz_stim{e,c};
        parfor n=1:nblocs
            [~,ids] = (sort(randi(2000,[length(seq) blocksz],'uint16')));
            repaux = repmat((seq),[1 blocksz]);
            repaux = repaux(ids); % shuffling each column
            repaux = num2cell(repaux,1);
            sur_transmat = cellfun(@transmat_from_seq,repaux,'uni',0);
            sur_transmat = cat(3,sur_transmat{:});
            sur_transmatall_conds{e,c,n} = sur_transmat;
        end
    end
end

%load('D:\My Documents\MATLAB\FastEnsDet\figures\4exp6conds\sur_transmat_1000_full_conds.mat')

%% concatenating data
cat_sur_transmatall = cell(4,1);
cat_sur_transmatall_conds = cell(4,6);

for e=1:4
    cat_sur_transmatall{e} = cat(3,sur_transmatall{e,:});
    for c=1:6
        cat_sur_transmatall_conds{e,c} = cat(3,sur_transmatall_conds{e,c,:});
    end
end

%% thresholds
prct = 90; % p <0.01

thr_all = cellfun(@(x) prctile(x,prct,3),cat_sur_transmatall,'uni',0);
thr_conds = cellfun(@(x) prctile(x,prct,3),cat_sur_transmatall_conds,'uni',0);

%% mapping the observed transition value to theshuffled distribution
tmat_pvals_stim = cell(nexps,nconds);
mean_sur_trans = tmat_pvals_stim;
for e=1:nexps
    nens = length(transmat_full{e});
    for c=1:nconds    
        mean_sur_val = mean(cat_sur_transmatall_conds{e,c},3);
        pvalmat = zeros(nens);  
        sur_pvalmat = pvalmat;
        for s=1:nens
            for s1=1:nens
                %pvalfun = @(x,y) reshape(mean(bsxfun(@le,x(:),y(:).'))*100,size(x));
                %pvalmat(s,s1) = pvalfun(cat_sur_transmatall_conds{e,c}(s,s1),transmat_stim{e,c}(s,s1));
                nless = sum(cat_sur_transmatall_conds{e,c}(s,s1,:)<transmat_stim{e,c}(s,s1));
                nequal = sum(cat_sur_transmatall_conds{e,c}(s,s1,:)==transmat_stim{e,c}(s,s1));
                pvalmat(s,s1) = 100 * (nless +0.5*nequal)/nsur;
                
                nless = sum(cat_sur_transmatall_conds{e,c}(s,s1,:)<mean_sur_val(s,s1));
                nequal = sum(cat_sur_transmatall_conds{e,c}(s,s1,:)==mean_sur_val(s,s1));
                sur_pvalmat(s,s1) = 100 * (nless +0.5*nequal)/nsur;
                
            end
        end
        tmat_pvals_stim{e,c} = pvalmat;
        mean_sur_trans{e,c} = sur_pvalmat;
    end
end
    
%%
bins = linspace(0,100,60);

nel_sur = histc(sur_pvalmat(:),bins);
nel_obs = histc(pvalmat(:),bins);
figure
subplot(211)
plot(bins,nel_obs,'r');hold on
plot(bins,nel_sur,'k');hold on

subplot(212)
[sort_obs,sortid] = sort(pvalmat(:),'descend');
plot(sort_obs,'ro-');hold on
plot(sur_pvalmat(sortid),'ko-');hold on
set(gca,'yscale','log')
%% binary transition matrices
bin_transmat_stim = cellfun(@(x,y) x>y,transmat_stim,thr_conds,'uni',0);
bin_transmat_full = cellfun(@(x,y) x>y,transmat_full,thr_all,'uni',0);
%% invariant transitions
nconds = 6;
invar_transmat_stim = cell(nexps,1);
for e=1:nexps
    invar_transmat_stim{e} = sum(cat(3,bin_transmat_stim{e,:}),3)==nconds;
end

%% strongly connected components
sparse_stim =cellfun(@sparse,bin_transmat_stim,'uni',0);
sparse_full =cellfun(@sparse,bin_transmat_full,'uni',0);
[~,sc_stim] = cellfun(@graphconncomp,sparse_stim,'uni',0);
[~,sc_full] = cellfun(@graphconncomp,sparse_full,'uni',0);
%% binarized and weighted in and out degree
w_in_deg_full = cellfun(@(x) sum(x),transmat_full,'uni',0);
w_in_deg_stim = cellfun(@(x) sum(x),transmat_stim,'uni',0);
b_in_deg_full = cellfun(@(x) sum(x),bin_transmat_full,'uni',0);
b_out_deg_full = cellfun(@(x) sum(x,2),bin_transmat_full,'uni',0);
b_in_deg_stim = cellfun(@(x) sum(x),bin_transmat_stim,'uni',0);
b_out_deg_stim = cellfun(@(x) sum(x,2),bin_transmat_stim,'uni',0);
%% Plotting raw and binarized data

for e=1:4
    figure
    for c=1:6
        subaxis(2,8,c)
        imagesc(log(transmat_stim{e,c}),log([min_transprob(e) max_transprob(e)]))
        axis square
        set(gca,'xticklabel','','yticklabel','');
        title(condnames{c})
        
        subaxis(2,8,c+8)
        imagesc(bin_transmat_stim{e,c})
        axis square
        set(gca,'xticklabel','','yticklabel','');
    end
    subaxis(2,8,7)
    imagesc(log(transmat_full{e}),log([min_transprob(e) max_transprob(e)]))
    axis square
    set(gca,'xticklabel','','yticklabel','');
    title('Full Recording')
    
    subaxis(2,8,8)
    imagesc(invar_transmat_stim{e}.*transmat_full{e})
    axis square
    set(gca,'xticklabel','','yticklabel','');
    title('Invariant Transition')
    
    subaxis(2,8,15)
    imagesc(bin_transmat_full{e})
    axis square
    set(gca,'xticklabel','','yticklabel','');
    
    subaxis(2,8,16)
    imagesc(invar_transmat_stim{e})
    axis square
    set(gca,'xticklabel','','yticklabel','');
    
    
end



%% plotting the graphs
%  circular coordinates
condid = [1 2 3 5 6 7];
figfold = 'D:\My Documents\MATLAB\FastEnsDet\figures\4exp6conds\';
for e=1%:4
    nnodes = length(bin_transmat_full{e});
    th = linspace(0,2*pi,nnodes);
    xpos = cos(th);
    ypos = sin(th);
    figure('units','normalized','outerposition',[0 0 1 1],'PaperPositionMode', 'auto','visible','on')
    for c=1:6
        subaxis(2,4,condid(c))
        plot_directed_graph(bin_transmat_stim{e,c},xpos,ypos,10,[]);
        title(condnames{c})     
        set(gca,'xtick',[],'ytick',[],'xticklabel','','yticklabel','','color',[240 248 255]./255);
        box on;
    end
    subaxis(2,4,4)
    plot_directed_graph(bin_transmat_full{e},xpos,ypos,10,[]);
    title('Whole')  
    set(gca,'xtick',[],'ytick',[],'xticklabel','','yticklabel','','color',[240 248 255]./255);
    box on;
    
    subaxis(2,4,8)
    plot_directed_graph(invar_transmat_stim{e},xpos,ypos,10,[]);
    title('Invariant')  
    set(gca,'xtick',[],'ytick',[],'xticklabel','','yticklabel','','color',[240 248 255]./255);
    box on;
    %print(gcf,'-dpng',[figfold,'Ensembles_trans_network_exp_',num2str(e),'.png'],'-r300')
    %close gcf;
  
end

%% graphs with sontrgly components 
for e=1:4
    nnodes = length(bin_transmat_full{e});
    th = linspace(0,2*pi,nnodes);
    xpos = cos(th);
    ypos = sin(th);
    figure('units','normalized','outerposition',[0 0 1 1],'PaperPositionMode', 'auto','visible','off')
    for c=1:6
        subaxis(2,4,condid(c))
        plot_directed_graph(bin_transmat_stim{e,c},xpos,ypos,10,sc_stim{e,c});
        title(condnames{c})     
        set(gca,'xtick',[],'ytick',[],'xticklabel','','yticklabel','','color',[240 248 255]./255);
        box on;
    end
    subaxis(2,4,4)
    plot_directed_graph(bin_transmat_full{e},xpos,ypos,10,sc_full{e});
    title('Whole')  
    set(gca,'xtick',[],'ytick',[],'xticklabel','','yticklabel','','color',[240 248 255]./255);
    box on;
    
    subaxis(2,4,8)
    plot_directed_graph(invar_transmat_stim{e},xpos,ypos,10,sc_full{e});
    title('Invariant')  
    set(gca,'xtick',[],'ytick',[],'xticklabel','','yticklabel','','color',[240 248 255]./255);
    box on;
    %print(gcf,'-dpng',[figfold,'Ensembles_trans_network_Components_exp_',num2str(e),'.png'],'-r300')
    %close gcf;
  
end

%% plotting ensemble rate vs in and out degree
ens_rate_full = cellfun(@(x) mean(x,2),ens_raster,'uni',0);

figure
for c=1:6
    subplot(3,6,c)
    for e=1:4
        plot(ens_rate{e}(:,c),w_in_deg_stim{e,c},'ko');hold on
    end
    
    subplot(3,6,c+6)
    for e=1:4
        plot(ens_rate{e}(:,c),b_in_deg_stim{e,c},'ko');hold on
    end
    subplot(3,6,c+12)
    for e=1:4
        plot(ens_rate{e}(:,c),b_out_deg_stim{e,c},'ko');hold on
    end
end
    
















