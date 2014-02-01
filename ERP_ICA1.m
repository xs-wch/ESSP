function [A,Z,nc_est] = ERP_ICA1(simu_EEG,res_pct)
% ERP_ICA1 - extract ERP using ICA
%            data concatated trial by trial
%            PCA-preprocessed, res_pct% power is kept      
%                
% Syntax:  [A,Z,nc_est] = ERP_ICA1(simu_EEG,res_pct)
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
% Other m-files required: eeglab toolbox
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
% EEG_mean = mean(simu_EEG,2);
% EEG_mean = repmat(squeeze(EEG_mean),len,1);
% EEG_mean = reshape(EEG_mean,chan,len,trial);

EEG_mean = mean(simu_EEG,2);
EEG_mean = repmat(EEG_mean,1,len);
EEG_demean = simu_EEG - EEG_mean;

%%%%%%%% PCA-preprocess
EEG_concat = EEG_demean(:,:);

S = eig(EEG_concat*EEG_concat');S = sort(S,'descend');
pow_sq = cumsum(S);
I_kept = pow_sq <= res_pct/100*pow_sq(end);

nc_est = sum(I_kept)+1;

[weights,sphere] = runica(EEG_concat,'pca',nc_est);
unmix = weights;
Z = zeros(chan,len);

Zt = unmix*EEG_concat;Zt = reshape(Zt,[nc_est,len,trial]);
Z(1:nc_est,:)  = mean(Zt,3);
A = zeros(chan,chan);
A(:,1:nc_est) = pinv(unmix);





