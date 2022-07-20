%% Computing entropy of sequences

[seq_ent,seq_maxent,dkl] = cellfun(@(x) compute_seq_ent(x),ens_seq_stim);
%% entropy of sequences of increasing word length
maxlen = 24;
minlen = 0;
[words,counts] = count_sequences(ens_seq_stim{1,4},minlen,maxlen);
[words_nz,counts_nz] = count_sequences(ens_seq_nz_stim{1,4},minlen,maxlen);
%%  computing with onlt observed states
nwords = maxlen;
w_ent = zeros(nwords,1);
w_dkl = w_ent;
w_ent_nz = w_ent;
w_dkl_nz = w_dkl;
maxent = w_dkl;
maxent_nz = maxent;
for w=1:nwords
%     nsyms = length(transmat_full{1})+1;
%     nsyms_w = nsyms^w;
    T = length(ens_seq_stim{1,4});
    T_w = T-w+1;
    T_nz = length(ens_seq_nz_stim{1,4});
    T_nz_w = T_nz -w +1;
    symbprobs = counts{w}./T_w; % probabilities of symbols
    symbprobs_nz = counts_nz{w}./T_nz_w; % probabilities of symbols
    %nonobsevs = nsyms_w - length(symbprobs);
    %symbprobs  = cat(1,symbprobs,zeros(nonobsevs,1));
    [w_ent(w),maxent(w),w_dkl(w)] = compute_seq_ent(length(symbprobs),symbprobs);
    [w_ent_nz(w),maxent_nz(w),w_dkl_nz(w)] = compute_seq_ent(length(symbprobs_nz),symbprobs_nz);
    
end
    
    
%% computing with all possible patterns  
nwords = maxlen;
all_w_ent = zeros(nwords,1);
all_w_dkl = all_w_ent;
all_w_ent_nz = all_w_ent;
all_w_dkl_nz = all_w_dkl;
all_maxent = all_w_dkl;
all_maxent_nz = all_maxent;

for w=1:nwords
    nsyms = length(unique(ens_seq_stim{1,4}));
    nsyms_nz = length(unique(ens_seq_nz_stim{1,4}));
    nsyms_w = nsyms^w;
    nsyms_w_nz = nsyms_nz^w;
    
    T = length(ens_seq_stim{1,4});
    T_w = T-w+1;
    T_nz = length(ens_seq_nz_stim{1,4});
    T_nz_w = T_nz -w +1;
    
    symbprobs = counts{w}./T_w; % probabilities of symbols
    symbprobs_nz = counts_nz{w}./T_nz_w; % probabilities of symbols
    
    nonobsevs = nsyms_w - length(symbprobs);
    nonobsevs_nz = nsyms_w_nz - length(symbprobs_nz);
    
%     symbprobs  = cat(1,symbprobs,zeros(nonobsevs,1));
%     symbprobs_nz  = cat(1,symbprobs_nz,zeros(nonobsevs_nz,1));
    [all_w_ent(w),all_maxent(w),all_w_dkl(w)] = compute_seq_ent(nsyms_w,symbprobs);
    [all_w_ent_nz(w),all_maxent_nz(w),all_w_dkl_nz(w)] = compute_seq_ent(nsyms_w_nz,symbprobs_nz);
    
end
    
  
%%
figure

subplot(121)
plot(w_ent,'bo-');hold on; 
plot(w_dkl,'ro-');hold on; 
plot(maxent,'go-');hold on; 

plot(w_ent_nz,'s-');hold on; 
plot(w_dkl_nz,'rs-');hold on; 
plot(maxent_nz,'gs-');hold on; 

title('Only observed words')
ylim([10^-6 maxlen])
%set(gca,'yscale','log','xscale','log')

subplot(122)
plot(all_w_ent,'bo-');hold on; 
plot(all_w_dkl,'ro-');hold on;
plot(all_maxent,'go-');hold on; 

plot(all_w_ent_nz,'bs-');hold on; 
plot(all_w_dkl_nz,'rs-');hold on; 
plot(all_maxent_nz,'gs-');hold on; 
title('All possible words')  
ylim([10^-0.5 10^2])
%set(gca,'yscale','log','xscale','log')
    
    
    
    
    
