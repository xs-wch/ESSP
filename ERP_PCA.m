function [A,Z,nc_est] = ERP_PCA(simu_EEG,res_pct)
% ERP_ICA1 - extract ERP using PCA
%            trial averaged first
%            res_pct% energy is kept      
%                
% Syntax:  [A,Z,nc_est] = ERP_PCA(simu_EEG,res_pct)
%
% Inputs:
%          simu_EEG: EEG signal, chan*len*trial
%          res_pct: percent of the energy kept in PCA preprocessing 
%
% Outputs:
%          A: spatial pattern
%          Z: temporal wave 
%          nc_est:  number of ERP components 
%
% Example: 
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% 

% Author: Chaohua Wu
% Department of Biomedical Engineering, Tsinghua University
% email: xs.wuchaohua@outlook.com
%
% Jan 29 ; first version
%------------- BEGIN CODE --------------

[chan,len,trial] = size(simu_EEG);

%%%%%%% remove base line
EEG_mean = mean(simu_EEG,2);
EEG_mean = repmat(EEG_mean,1,len);
EEG_mean = reshape(EEG_mean,chan,len,trial);
EEG_demean = simu_EEG - EEG_mean;

%%%%%%%% PCA-preprocess
EEG_trialavr = mean(EEG_demean,3);

[U,S,V] = svd(EEG_trialavr);
pow_sq = cumsum(diag(S).^2);
I_kept = cumsum(pow_sq) <= res_pct/100*pow_sq(end);

nc_est = sum(I_kept);

A = U(:,1:nc_est);
Z = S(1:nc,1:nc)*V(:,1:nc_est)';
end