% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

clear;clc;close all;
%% Parameters for raster
N = 100; %Number of neurons
fr = 0.15; %Firing rate of neurons
nens = 6; %Number of ensembles
time_inactivity = 0.2; %Percentage of time without ensemble
ncellsperens = 20; %Number of neurons per ensemble
% T_min = 1000; %Minimum value of frames
% T_max = 5000; %Maximum value of frames
% T_size = 5; %Number of different values of frames
% T = floor(linspace(T_min,T_max,T_size)); %Vector with number of frames (variable)
% T = [500 1000 2500 5000 7500 10000 25000 50000];
T = [500 1000 2500];
T_size = length(T);
R = 4; %Number of repetitions

%% Pre-allocated memory
performance_lc = cell(T_size,R); %Matrix with performance metrics using TPR and FPR
performance_rh= performance_lc; %Matrix with performance metrics using Jaccard distance

nens_alg_lc = zeros(T_size,R); %Vector with number of ensembles detected by algorithms
nens_alg_rh = nens_alg_lc;

times_comp_lc = zeros(T_size,R); %Vector with computing time of algorithms
times_comp_rh =times_comp_lc ; 

%%
t_init = tic;
for t = 1:T_size        
    t_0 = tic;    
    tsel = T(t);
    for r = 1:R
        %% Generate raster
        [ensmat_in,enscells_in,raster,frates] = generate_data(N,tsel,nens,ncellsperens,fr,time_inactivity);        
        %% Ensembles detection and sorting
        raster = single(raster);
        % LC
        [ensmat_in_alg, enscells_in_alg, nens_alg_lc(t,r), times_comp_lc(t,r)] = get_carrillo_ens(raster);        
        [ensmat_lc,enscells_lc] = sort_by_equivalence_rh(ensmat_in,ensmat_in_alg,enscells_in_alg,'corr');        
        % RH
        [ensmat_in_alg, enscells_in_alg, nens_alg_rh(t,r), times_comp_rh(t,r)] = get_herzog_ens(raster);                
        [ensmat_rh,enscells_rh] = sort_by_equivalence_rh(ensmat_in,ensmat_in_alg,enscells_in_alg,'corr');
        %% Get performance
        % LC
        performance_lc{t,r} = get_performance_rh(single(ensmat_in),single(enscells_in),single(ensmat_lc),single(enscells_lc));
        % RH
        performance_rh{t,r} = get_performance_rh(single(ensmat_in),single(enscells_in),single(ensmat_rh),single(enscells_rh));

    end
    t_i = toc(t_0);
    disp([' All reps done in ', num2str(t_i),' seconds'])
%     save('checkpoint_ensembles_v3.mat','performance_lc','performance_rh',...
%         'nens_alg_lc','nens_alg_rh','times_comp_lc','times_comp_rh','T')

end
t_elapsed = toc(t_init);
disp(['Total time elapsed: ', num2str(t_elapsed),' seconds'])


%% Plot performances curves
legends = {'SVD-based', 'Density-based'};
x_label = 'T (bins)';nmets = 2;
global_seq_sim_lc = zeros(T_size,R);
indivi_seq_sim_lc = zeros(nens,T_size,R);
core_sim_lc = indivi_seq_sim_lc;

global_seq_sim_rh = zeros(T_size,R);
indivi_seq_sim_rh = zeros(nens,T_size,R);
core_sim_rh = indivi_seq_sim_rh;
met = 'cor';

for t = 1:T_size
    for r=1:R
        global_seq_sim_lc(t,r) = performance_lc{t,r}.([met,'_glob']);
        indivi_seq_sim_lc(:,t,r) = performance_lc{t,r}.([met,'_indiv']);
        core_sim_lc(:,t,r) = performance_lc{t,r}.([met,'_core']);
        
        global_seq_sim_rh(t,r) = performance_rh{t,r}.([met,'_glob']);
        indivi_seq_sim_rh(:,t,r) = performance_rh{t,r}.([met,'_indiv']);
        core_sim_rh(:,t,r) = performance_rh{t,r}.([met,'_core']);
    end
end

% Plotting performance curves
% Plotting only correlation measure
figure
subplot(3,2,1)
h1=shadedErrorBar(T,mean(times_comp_lc,2),std(times_comp_lc,0,2),'lineprops',{'o-'});hold on
h2=shadedErrorBar(T,mean(times_comp_rh,2),std(times_comp_rh,0,2),'lineprops',{'o-'});hold on
ylabel('Processing Time (sec)')
legend([h1.mainLine h2.mainLine],legends,'location','northwest')


subplot(3,2,2)
shadedErrorBar(T,mean(nens_alg_lc,2),std(nens_alg_lc,0,2),'lineprops',{'o-'});hold on
shadedErrorBar(T,mean(nens_alg_rh,2),std(nens_alg_rh,0,2),'lineprops',{'o-'});hold on
plot([T(1) T(end)],[nens nens],'r--');hold on
ylabel('Detected Ensembles')


subplot(3,2,3)
shadedErrorBar(T,nanmean(global_seq_sim_lc,2),nanstd(global_seq_sim_lc,0,2),'lineprops',{'o-'});hold on
shadedErrorBar(T,nanmean(global_seq_sim_rh,2),nanstd(global_seq_sim_rh,0,2),'lineprops',{'o-'});hold on
ylabel('Global Sequence Correlation')
ylim([0 1])


subplot(3,2,4)
shadedErrorBar(T,mean(nanmean(indivi_seq_sim_lc),3),std(nanmean(indivi_seq_sim_lc),0,3),'lineprops',{'o-'});hold on
shadedErrorBar(T,mean(nanmean(indivi_seq_sim_rh),3),std(nanmean(indivi_seq_sim_rh),0,3),'lineprops',{'o-'});hold on
ylabel({'Mean Ensemble','Sequence Similarity'})
xlabel(x_label)
ylim([0 1])

subplot(3,2,5)
shadedErrorBar(T,mean(nanmean(core_sim_lc),3),std(nanmean(core_sim_lc),0,3),'lineprops',{'o-'});hold on
shadedErrorBar(T,mean(nanmean(core_sim_rh),3),std(nanmean(core_sim_rh),0,3),'lineprops',{'o-'});hold on
ylabel({'Mean Ensemble','Core-Cells Similarity'})
xlabel(x_label)
ylim([0 1])

