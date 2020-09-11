function [rand_nums] = rand_with_int(N,delta)
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