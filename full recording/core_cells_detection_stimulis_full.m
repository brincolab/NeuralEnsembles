%% Computign core-cells for each stimulus and for full recording
load('/home/octodon/Documents/MATLAB/FastEnsDet/4exps/4_exps_6_conds_concatenated.mat')
load('/home/octodon/Documents/MATLAB/raster2Ens/Analyzing4Exp_FullRecording/v4_4_exps_ens_seq.mat')

%% raster and ensemble sequence for stimuli
ras_ind = cumsum(cat(2,ones(4,1),raster_lengths),2);
ens_seq_stim = cell(4,6);
ras_stim = ens_seq_stim;
for e=1:4
    for c=1:6
        ens_seq_stim{e,c} = ens_seq{e}(ras_ind(e,c):ras_ind(e,c+1)-1);
        ras_stim{e,c} = cat_rasters{e}(:,ras_ind(e,c):ras_ind(e,c+1)-1);
    end
end
%% Ensemble raster
ensras_stim = cellfun(@(x) bsxfun(@eq,x',(1:max(x))'),ens_seq_stim,'uni',0);
ensras_full = cellfun(@(x) bsxfun(@eq,x',(1:max(x))'),ens_seq,'uni',0);
%% Core cell computation
nsur = 1000;
prct = 99.9;
[stim_neuronid,stim_idthr,stim_ens_cel_corr] = cellfun(@(x,y) ...
    find_core_cells_by_correlation(x*1,y*1,nsur,prct),ras_stim,ensras_stim,'uni',0);
[full_neuronid,full_idthr,full_ens_cel_corr] = cellfun(@(x,y) ...
    find_core_cells_by_correlation(x*1,y*1,nsur,prct),cat_rasters,ensras_full,'uni',0);