function [raster] = ts2binaryraster(spks,bin,sr)
% [raster,poprate] = ts2raster_poprate(spks,bin,sr,op)
%
% converts timestamps (timepoints) of N x 1 (or transpose) cell array
% 'spks' samples at 'sr' sampling rate into a
% binary (logical) raster 'raster' with bin size 'bin' s.

N = length(spks); % Number of units
maxt = cellfun(@(x) max(x),spks,'uniformoutput',0); % maximum time for all the units
nonempt = cellfun(@(x) ~isempty(x),spks);
maxt = max(cell2mat(maxt(nonempt)));

bin = bin*sr; % bins in timepoints
tt = [0:bin:maxt-bin]; % timebase for raster

raster = false(N,length(tt));
for n=1:N        
    raster(n,:) = (histc(spks{n},tt));    
end
raster = logical(raster); % binarizing



