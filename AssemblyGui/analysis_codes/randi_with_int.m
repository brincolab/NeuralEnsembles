function [rand_nums] = randi_with_int(N,nmax,delta)
% N = 10; % number of integers required
% delta = 2; % minimum difference required

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

a = randperm(nmax);
idx = 1;
rand_nums = a(idx);
while(length(rand_nums) < N && idx < length(a))
    idx = idx+1;
    c = abs(rand_nums - a(idx));
    if any(c < delta)
        continue;
    end
    rand_nums = [rand_nums; a(idx)];
end
