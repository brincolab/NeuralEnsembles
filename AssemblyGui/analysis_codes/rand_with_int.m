function [rand_nums] = rand_with_int(N,delta)
% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

N = 10; % number of integers required
delta = 2; % minimum difference required

a = randperm(100);
idx = 1;
b = a(idx);
%%
while(length(b) < N && idx < length(a))
    idx = idx+1;
    c = abs(b - a(idx));
    if any(c < delta)
        continue;
    end
    b = [b; a(idx)];
end

b