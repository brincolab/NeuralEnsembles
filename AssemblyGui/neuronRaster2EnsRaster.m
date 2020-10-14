function [ens_seq,maxcor,templates,core_cells,clust_out,pars] = neuronRaster2EnsRaster(spktimes,filename,pars)
%
% [ens_seq,templates,clust_out] = neuronRaster2EnsRaster(spktimes,pars)
%
% Generates ensemble raster (N_ens x T) from neuron raster (N_neu x T),
% whose timestamps are contained on the cell array spktimes.
% The ensemble sequences is on ens_seq vector.
%
% INPUTS
% 'spktimes' is a N_neur x 1 cell array, where spktimes(n) containes the
% timestamps of the n-th neuron on time points (in sampling rate
% resolution)
% 'filename' is a string with a name for the file to write files to disk
% 'pars' is a structure with the following fields
%   'bin' the bin size for the raster in seconds
%   'sr' is the sampling rate in Hz.
%   'minspk' is the size for a population spike
%   'sampfac' is the number of subsamples
%   'dc' is the distance cut-off for the clustering procedure
%   'ishalo' is 1 if halo cluster refinement is osed and 0 otherwise
%   'thrmet' is the method used for detecting centroids. Methods can be:
%       'manual' NOT RECOMMENDED. Ask for manual selection of centroids on each iteration
%        'sur' computes a surrogate distribution for the product delta*rho
%        'cutoff' computes the cutoff on the spectrum of delta*rho
%        'fit' fit a power-law to the rho vs delta plot and uses prediction
%        bounds and threshold. RECOMMENDED!
%        'hist' computes the 2d histogram of delta and rho and uses median + mad as threshold
%   'maxmem' max memory to use in bytes. Should be at least 1Gb lower than
%   the total RAM available.
%   'npcs' is the number of principal components used in clustering
%   'ccthr' is the minimal correlation for template matching
%   'nsur' number of surrogate for core-cell detection
%   'prct' percentile of the surrogate distribution used as threshold for
%   core-cell detection.
%
% OUTPUTS
% 'ens_seq' T x 1 vector where ens_seq(t)=e, if ensemble 'e' is active on
% time 't'
% 'maxcorr' is a T x 1 where maxcorr(t)=cval, where cval is the highest
% correlation of all the templates and the bin 't'.
% 'templates' is a N_neur x N_ens, where templates(n,e) is a vector with
% the weight of neuron 'n' on the ensemble 'e'.
% 'core_cells' is a N_ens x 1 cell array with the indices of core neurons
% for each ensemble.
% 'clust_out' is a the output of the global clustering step with the following fields:
%   'rho' is the density of each datapoint 
%   'delta' is the minimum distance to the next higher density point
%   'ensId' is a vector where ensId(s)=c, means that sample 's' belongs to
%   cluster 'c'
%   'sampcent' are the templates used on clustering
%   'cents' are the indices of the centroids
%   'predbounds' is the treshold on the delta vs rho plot
% 'pars' is the set of parameters used in the analysis.
%
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

% Checking input parameters. If any input parameter is absent, will use
% default
[pars] = ensRast_checkPars(pars); % checks the input parameters and set the to default if any is absent
disp('Analysis will run with the following parameters')
disp(pars);% displays parameters used

% Generating raster
[raster] = ts2binaryraster(spktimes,pars.bin,pars.sr);

% Subsampling raster and clustering
[~,ensId,centInd,delta,rho,~,sampCent,templates,predbounds] = ...
    subSample_ensembles_pca(raster,pars.minspk,pars.sampfac,pars.dc,...
    pars.ishalo,pars.thrmet,pars.maxmem,pars.npcs,filename);

% Template Matching
[ens_seq,maxcor] = ens_raster_from_templates(templates,raster,pars.ccthr);

% Core-cells detection
core_cells = arrayfun(@(x) find(find_core_cells(raster(:,ens_seq==x),pars.nsur,pars.prct)),1:max(ens_seq),'uni',0);

% Generating outputs
clust_out.rho = rho;
clust_out.delta = delta;
clust_out.ensId = ensId;
clust_out.sampcent = sampCent;
clust_out.cents = centInd;
clust_out.predbounds = predbounds;


function [pars] = ensRast_checkPars(pars)
% checks the inputs parameters for the neuronRaster2EnsRaster function
if ~isfield(pars,'bin') || isempty(pars.bin) % checking bin size
    pars.bin = 0.02; % default
    warning(['Default bin size = ',num2str(pars.bin),' s is used'],'off')
end
if ~isfield(pars,'sr') || isempty(pars.sr) % checking bin size
    pars.sr = 20000; % default
    warning(['Default sampling rate = ',num2str(pars.sr),' Hz is used'])
end
if ~isfield(pars,'minspk') || isempty(pars.minspk) % checking min pop spike size
    pars.minspk = 3; % default
    warning(['Default minspk = ',num2str(pars.minspk),' spikes is used'])
end
if ~isfield(pars,'sampfac') || isempty(pars.sampfac) % checking number of samples
    pars.sampfac = 100; % default
    warning(['Default number of samples = ',num2str(pars.sampfac),' is used'])
end
if ~isfield(pars,'thrmet') || isempty(pars.thrmet) % checking centroid detection method
    pars.thrmet = 'fit'; % default
    warning(['Default ',(pars.thrmet),' centroid searching method is used'])
end
if ~isfield(pars,'maxmem') || isempty(pars.maxmem) % checking memory to keep unused
    pars.maxmem = 4*10^9; % default 4Gb
    warning(['Default maximal memory use of ',num2str(pars.maxmem),' bytes is used'])
end
if ~isfield(pars,'npcs') || isempty(pars.npcs) % checking number pf principal components
    pars.npcs = 6; % default
    warning(['Default number of principal components = ',num2str(pars.npcs),' is used'])
end
if ~isfield(pars,'ccthr') || isempty(pars.ccthr) % checking minimal correlation for template matching step
    pars.ccthr = 0.1; % default
    warning(['Default minimal correlation for template matching = ',num2str(pars.ccthr),' is used'])
end
if ~isfield(pars,'nsur') || isempty(pars.nsur) % checking minimal correlation for template matching step
    pars.nsur = 1000; % default
    warning(['Default number of surrogates = ',num2str(pars.nsur),' is used'])
end
if ~isfield(pars,'prct') || isempty(pars.prct) % checking minimal correlation for template matching step
    pars.prct = 99; % default
    warning(['Default percentile on surrogate distribution = ',num2str(pars.prct),...
        ' is used as threshold for core-cells detection'])
end
if ~isfield(pars,'dc') || isempty(pars.dc) % checking minimal correlation for template matching step
    pars.dc = 0.01; % default
    warning(['Default cut-off parameter = ',num2str(pars.prct),' is used '])
end
if ~isfield(pars,'ishalo') || isempty(pars.ishalo) % checking minimal correlation for template matching step
    pars.ishalo = 0; % default
    warning(['Not using halo cluster refinement by default'])
end