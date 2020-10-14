function [ensmat_in,enscells_in,raster,frates] = generate_data(N,T,nens,ncellsperens,fr,timesperens,show_plot)
% Generate a synthetical raster with pre-defined ensembles.
% INPUT:
%   N : <integer> with the number of neurons in the population.
%   T: <integer> with the number of total frames.
%   nens: <integer> with the number of ensembles.
%   ncellsperens: <integer> Nens-by-1 vector with the number of core cells per ensemble.
%       If it's a scalar then all ensembles have the same number of core cells. If it's equal to 0, 
%       the number is generated randomly. 
%   fr: <float> with the "firing rate" (actually represents a variance in a pdf that generates 
%       the raster).    
%   timesperens: <float> Nens-by-1 vector with the proportional time that each ensemble is
%       actived. If it's a scalar then represents the proportional time that no ensemble
%       is active.
% OUTPUT:
%   ensmat_in: <logical> Nens-by-T matrix with the activity of each ensemble.
%   enscells_in: <logical> N-by-Nens matrix with the core cells of each ensemble.
%   raster: <logical> N-by-T matrix with the raster.
%   frates: <float> N-by-1 vector with the firing rate of each neuron.

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

%% Formatting parameters
%Converts Ncellsperens in a vector
if isscalar(ncellsperens)
    if ncellsperens > 0
        ncellsperens_vec = ncellsperens*ones(nens,1);
    else
        random = randi([10,floor(N*fr)]);
        ncellsperens_vec = random*ones(nens,1);
    end
else
    ncellsperens_vec = ncellsperens;
end

% Converts timesperens in a vector
if isscalar(timesperens)
    proportion = (1-timesperens)/nens;
    timesperens_vec = proportion*ones(nens,1);
else
    timesperens_vec = timesperens;
end

% Generate raster with ensembles
[ensmat_in,enscells_in,raster,frates] = MakeEnsembles_fix_rate(N,fr,T,nens,ncellsperens_vec,timesperens_vec);

% Plot sequence ensemble activation and core-cells 
if exist('show_plot','var') && show_plot
    plot_ens(ensmat_in, enscells_in, ncellsperens, 'Generated');
end

end