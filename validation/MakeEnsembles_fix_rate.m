function [ensmat_in,enscells_in,raster,frates] = MakeEnsembles_fix_rate(ncells,fr,T,nens,ncellsperens,ntimesperens)

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020
if sum(ntimesperens)>1
    error('sum(ntimesperens) must be lower than 1')
    
elseif any(ncellsperens>ncells)
    error('ncellsperens be lower than ncells')
    
% elseif strcmp(method,'pearson')==0 && strcmp(method,'cos')==0 && strcmp(method,'propeq')==0
%     error('Wrong similarity method')
% 
% elseif numel(fr)~=ncells && numel(fr)>1
%     error('For a fr vector length(fr) must be equal to ncells')
% 
% elseif any(fr>1) % max allowed firing rate is 0.6
%     error('fr can''t be greater than 0.9')
        
end

% 1.- Generating the firing rates vector
tt = 1:T; % timebase of the raster
if numel(fr)>1 % fr is a vector
    frates = reshape(fr,[ncells 1]); % just to be sure that has the right shape    

else % fr is scalar
    frates = 1;
    while length(frates)<ncells % looping while the number of positive entries of frates is lower than ncells
        frates = fr*randn(ncells*10,1); % normally distributed firing rates with mean avefr
        frates(frates<0.001 | frates>0.5 ) =[]; % removing negative firing rates and very high firing rates
    end
    frates = frates(randperm(length(frates),ncells)); % extracting random values from the positive part of the distribution
end

% 2.- Generating ensembles (cells participating), activation time and plugging ensembles on raster
raster = false(ncells,T);
enscells_in = false(ncells,nens); % binary matrix with (i,j)=1 if neuron i participates on ensemble j
ensmat_in = false(nens,T); % binary matrix with (i,j)=1 if ensemble i was active on the significant bin 
tt2 = tt;
tmax = T;
for i=1:nens
    % 2.1 enerating the ensemble cells
    cells = randperm(ncells,ncellsperens(i)); % ensemble cells indices
    enscells_in(cells,i) = true;
    
    % 2.2 generating the ensemble activation times
    enspos = randperm(tmax,floor(ntimesperens(i)*T)); % generating random indices for the ensembles activation
    ensmat_in(i,tt2(enspos)) = true; % assigning to the activation binary matrix
    
    % 2.3 Plugging the ensemble on the raster
    raster(:,tt2(enspos)) = repmat(enscells_in(:,i),[1 floor(ntimesperens(i)*T)]); % plugging the ensemble
    
    tt2(enspos)=[]; % removing those indices from the timebase to avoid repetitions
    tmax = length(tt2); % maximum number of indices after removal
end

% 3.- Checking excess or deficit of spikes
fr_ens = sum(raster,2);
E_fr = floor(frates.*T);
fr_dif = fr_ens - E_fr; % negative is excess and positive deficit
rem_spk = find(fr_dif>0); % positive means spike removal
add_spk = find(fr_dif<0); % negative means spike addition

% 3.2 Adding spikes
for n=1:length(add_spk)
    no_spks = find(raster(add_spk(n),:)==0); % bins without spikes
    add_times = randperm(length(no_spks),abs(fr_dif(add_spk(n)))); % abs(fr_dif) random indices
    raster(add_spk(n),no_spks(add_times))=true;
end

% 3.3 removing spikes
for n=1:length(rem_spk)
    yes_spks = find(raster(rem_spk(n),:)); % bins with spikes
    rem_times = randperm(length(yes_spks),abs(fr_dif(rem_spk(n)))); % abs(fr_dif) random indices
    raster(rem_spk(n),yes_spks(rem_times))=false;
end

end
