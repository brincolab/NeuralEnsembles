function [ensmat_in, enscells_in, nens, time_comp] = get_carrillo_ens(raster, param)
% Gets the activation sequence and core-cells by Carrillo's algorithm.
% INPUT:
%   raster: <logical> N-by-T matrix with the raster. N is the number of total neurons
%     and T is the number of total frames.
%   param: [optional] <struct> with the parameters to call the detection algorithm.
% OUTPUT:
%   ensmat_in: <logical> nens-by-T matrix with the activity of each ensemble.
%   enscells_in: <logical> N-by-nens matrix with the core cells of each ensemble.
%   nens: <integer> with the number of ensembles detected.
%   time_comp: <float> with the time of computing in detection algorithm.

% If param variable don't exist, then it's created by default

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020
if ~exist('param','var')
%     param = struct('pks',4,'ticut',[],'jcut',[],'state_cut',[]); % automatic parameters based on surrogate data
%     param = struct('pks',10,'ticut',0.3,'jcut',0.08,'state_cut',10); % seems optimal between computer performance and precision with synth data
    param = struct('pks',3); % Default parameters 
end

% Detect ensembles and core-cells with Carrillo's algorithm
tic;
[core_svd, state_pks_full] = findSVDensemble2(raster, param);
time_comp = toc;

nens = size(core_svd,1); %Number of ensembles detected
N = size(raster,1); %Number of total neurons
T = size(raster,2); %Number of total frames (bins)

% Convert format from core_svd to enscells_in
enscells_in = false(N, nens); 
for i = 1:nens
    Ncellsperens = size(core_svd{i},1);
    for j = 1:Ncellsperens
        cell_index = core_svd{i}(j);
        enscells_in(cell_index,i) = true;
    end
end

%Convert format from state_pks_full to ensmat_in
ensmat_in = false(nens, T); 
for i = 1:T
    if (state_pks_full(i) ~= 0)
       ensmat_in(state_pks_full(i),i) = true;
    end
end

end