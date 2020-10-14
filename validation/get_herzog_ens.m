function [ensmat_in, enscells_in, nens, time_comp] = get_herzog_ens(raster, param)
% Gets the activation sequence and core-cells by Herzog's algorithm.
% INPUT:
%   raster: <logical> N-by-T matrix with the raster, where N is the number of total neurons
%     and T is the number of total frames.
%   param: [optional] <struct> with the parameters to call the detection algorithm.
% OUTPUT:
%   ensmat_in: <logical> nens-by-T matrix with the activity of each ensemble.
%   enscells_in: <logical> N-by-nens matrix with the core cells of each ensemble.
%   nens: <integer> with the number of ensembles detected.
%   time_comp: <float> with the time of computing in detection algorithm.

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

% If param variable don't exist, then it's created by default
if(~exist('param','var'))
    param = struct('dc',[],'npcs',[],'minspk',[],'nsur',[],'prct',[],'cent_thr',[],'inner_corr',[],'minsize',[]);
end

% Detect ensembles and core-cells with Herzog's algorithm
tic;
[ensmat_in, enscells_in] =  raster2ens_by_density(raster, param);
time_comp = toc;

nens = size(enscells_in,2); %Number of ensembles detected

end