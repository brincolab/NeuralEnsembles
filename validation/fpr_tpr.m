function [fpr, tpr] = fpr_tpr(data_real, data_est)
% Calculate the FPR and TPR values.
% INPUT:
%   data_real: <logical> D-by-T matrix with the real data, where D is the number of data to 
%     compare in a trial and T is the number of trials.
%   data_est: <logical> D-by-T matrix with the estimated data, where D is the number of data to 
%     compare in a trial and T is the number of trials.
% OUTPUT:
%   fpr: <float> 1-by-T vector with FPR (False Positive Rate) between data_real columns
%     and data_est columns.
%   tpr: <float> 1-by-T vector with TPR (True Positive Rate) between data_real columns
%     and data_est columns.

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020

% Get size of data_real and data_est
% dim_real = size(data_real);
% dim_est = size(data_est);
% 
% % Check if both size are equal
% if(~isequal(dim_real,dim_est))
%     error('Size of data_real and data_est must be equal');
% end
% 
% % Check if both data is logical
% if(~islogical(data_real) || ~islogical(data_est))
%     error('Inputs must be logical arrays');
% end

% Calculate TPR and FPR
fpr = sum(data_est & ~data_real)./sum(~data_real);
tpr = sum(data_est & data_real)./sum(data_real);

end