% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

clear;clc;close all;
%% Parameters for raster
T = 500; %Number of frames
fr = 0.15; %Firing rate of neurons
nens = 12; %Number of ensembles
time_inactivity = 0.2; %Percentage of time without ensemble
ncellsperens = 35; %Number of neurons per ensemble
N_min = 100; %Minimum value of neuron population
N_max = 200; %Maximum value of neuron population
N_size = 3; %Number of different values of neuron population
N = floor(linspace(N_min,N_max,N_size))'; %Vector with number of neurons population (variable)
R = 4; %Number of repetitions

%% Pre-allocated memory
performance_lc = cell(N_size,R); %Matrix with performance metrics using TPR and FPR
performance_rh= performance_lc; %Matrix with performance metrics using Jaccard distance

nens_alg_lc = zeros(N_size,R); %Vector with number of ensembles detected by algorithms
nens_alg_rh = nens_alg_lc;

times_comp_lc = zeros(N_size,R); %Vector with computing time of algorithms
times_comp_rh =times_comp_lc ; 
%%
raster_s = cell(N_size,R);
t_init = tic;
for n = 1:N_size        
    t_0 = tic;
    nsel = N(n);
    for r = 1:R
        %% Generate raster
        [ensmat_in,enscells_in,raster,frates] = generate_data(nsel,T,nens,ncellsperens,fr,time_inactivity);
        raster_s{n,r} = raster;
        %% Ensembles detection and sorting
        % LC
        [ensmat_in_alg, enscells_in_alg, nens_alg_lc(n,r), times_comp_lc(n,r)] = get_carrillo_ens(raster);        
        [ensmat_lc,enscells_lc] = sort_by_equivalence_rh(ensmat_in,ensmat_in_alg,enscells_in_alg,'corr');        
        % RH
        [ensmat_in_alg, enscells_in_alg, nens_alg_rh(n,r), times_comp_rh(n,r)] = get_herzog_ens(raster);                
        [ensmat_rh,enscells_rh] = sort_by_equivalence_rh(ensmat_in,ensmat_in_alg,enscells_in_alg,'corr');
        %% Get performance
        % LC
        performance_lc{n,r} = get_performance_rh(single(ensmat_in),single(enscells_in),single(ensmat_lc),single(enscells_lc));
        % RH
        performance_rh{n,r} = get_performance_rh(single(ensmat_in),single(enscells_in),single(ensmat_rh),single(enscells_rh));

    end
    t_i = toc(t_0);
    disp([' All reps done in ', num2str(t_i),' seconds'])

end
t_elapsed = toc(t_init);
disp(['Total time elapsed: ', num2str(t_elapsed),' seconds'])

%% Extracting performance
nmets = 2;
global_seq_sim_lc = zeros(N_size,R);
indivi_seq_sim_lc = zeros(nens,N_size,R);
core_sim_lc = indivi_seq_sim_lc;

global_seq_sim_rh = zeros(N_size,R);
indivi_seq_sim_rh = zeros(nens,N_size,R);
core_sim_rh = indivi_seq_sim_rh;
met = 'cor';

for n = 1:N_size
    for r=1:R
        global_seq_sim_lc(n,r) = performance_lc{n,r}.([met,'_glob']);
        indivi_seq_sim_lc(:,n,r) = performance_lc{n,r}.([met,'_indiv']);
        core_sim_lc(:,n,r) = performance_lc{n,r}.([met,'_core']);
        
        global_seq_sim_rh(n,r) = performance_rh{n,r}.([met,'_glob']);
        indivi_seq_sim_rh(:,n,r) = performance_rh{n,r}.([met,'_indiv']);
        core_sim_rh(:,n,r) = performance_rh{n,r}.([met,'_core']);
    end
end

% Plotting performance curves
% Plotting only correlation measure
legends = {'Carrillo', 'Herzog'};
figure
subplot(3,2,1)
shadedErrorBar(N,mean(times_comp_lc,2),std(times_comp_lc,0,2),'lineprops',{'o-'});hold on
shadedErrorBar(N,mean(times_comp_rh,2),std(times_comp_rh,0,2),'lineprops',{'o-'});hold on
ylabel('Processing Time (sec)')


subplot(3,2,2)
shadedErrorBar(N,mean(nens_alg_lc,2),std(nens_alg_lc,0,2),'lineprops',{'o-'});hold on
shadedErrorBar(N,mean(nens_alg_rh,2),std(nens_alg_rh,0,2),'lineprops',{'o-'});hold on
plot([N(1) N(end)],[nens nens],'r--');hold on
ylabel('Detected Ensembles')


subplot(3,2,3)
shadedErrorBar(N,nanmean(global_seq_sim_lc,2),nanstd(global_seq_sim_lc,0,2),'lineprops',{'o-'});hold on
shadedErrorBar(N,nanmean(global_seq_sim_rh,2),nanstd(global_seq_sim_rh,0,2),'lineprops',{'o-'});hold on
ylabel('Global Sequence Correlation')
ylim([0 1])


subplot(3,2,4)
shadedErrorBar(N,mean(nanmean(indivi_seq_sim_lc),3),std(nanmean(indivi_seq_sim_lc),0,3),'lineprops',{'o-'});hold on
shadedErrorBar(N,mean(nanmean(indivi_seq_sim_rh),3),std(nanmean(indivi_seq_sim_rh),0,3),'lineprops',{'o-'});hold on
ylabel({'Mean Ensemble','Sequence Similarity'})
xlabel('Network Size')
ylim([0 1])

subplot(3,2,5)
shadedErrorBar(N,mean(nanmean(core_sim_lc),3),std(nanmean(core_sim_lc),0,3),'lineprops',{'o-'});hold on
shadedErrorBar(N,mean(nanmean(core_sim_rh),3),std(nanmean(core_sim_rh),0,3),'lineprops',{'o-'});hold on
ylabel({'Mean Ensemble','Core-Cells Similarity'})
xlabel('Network Size')
ylim([0 1])

%%
subplot(3,2,6)
shadedErrorBar(N,mean(times_comp_lc,2),std(times_comp_lc,0,2),'lineprops',{'o-'});hold on
shadedErrorBar(N,mean(times_comp_rh,2),std(times_comp_rh,0,2),'lineprops',{'o-'});hold on
ylabel('Processing Time (sec)')
xlabel('Network Size')
legend(legends)
