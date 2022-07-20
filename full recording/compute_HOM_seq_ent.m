function [seq_ent,seq_maxent,dkl] = compute_HOM_seq_ent(nsymbs,symbprobs)
% Computes the entropy for a k-th order markov chain sequence
% seq is a sequence of integers (symbols) going from 0 to max(seq)
%
%nsymbs = length(symbprobs);
%minprob = eps;
%symbprobs(symbprobs==0) = minprob;
symbprobs=symbprobs(symbprobs>0);
seq_maxent = -(nsymbs*(1/nsymbs)*log2(1/nsymbs)); % entropy of uniform seq
seq_ent = -(sum(symbprobs.*log2(symbprobs))); % entropy of sequence
% dkl is just the difference between both entropies given thsat reference
% distribution is uniform
dkl = seq_maxent - seq_ent;

%prob_ratio = symbprobs./repmat(1/nsymbs,size(symbprobs)); % ratio of probabilities
%dkl = sum(symbprobs.*log2(prob_ratio));% Kullback-Leibler divergence