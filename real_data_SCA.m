function [A,Z,SF,nc_est] = real_data_SCA(simu_EEG,ERP_template)
% ERP_ICA1 - extract ERP using SCA
%            the template is the square template
%            
%            
%                
% Syntax:  [A,Z,nc_est] = ERP_SCA(simu_EEG1,simu_EEG2,sqtempl)
%
% Inputs:
%          simu_EEG 1/2: EEG signal, chan*len*trial
%          nc: number of ERP components used
%          sq_templ: the square wave as template
%          
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
% Jun 30 ; first trial average, then use sca to extract ERP. To do this to reduce computation load 
% Jun 18 ; first version
%------------- BEGIN CODE --------------



if nargin  < 1
	%disp('at least 2 parameters required');
    addpath /home/test2/wch_spf_ERP/extr_onesrc2
	[simu_EEG,SP,SF,Z]= sERP_data_gen(-10,45);
    load sq_template1.mat
end

[chan,len,trial] = size(simu_EEG);
%[chan2,len2,trial2] = size(simu_EEG2);
% if (chan1~=chan2)||(len1~=len2)
% 	disp('ERP data dimension does not match')
% 	pause
% end
nc = size(ERP_template,1);


% simu_EEG1_2D = simu_EEG1(:,:);
% simu_EEG2_2D = simu_EEG2(:,:);
%%%%%%% remove base line
% EEG_mean1 = mean(simu_EEG1_2D,2);
% EEG_mean1 = repmat(EEG_mean1,1,len1*trial1);
% simu_EEG1= simu_EEG1_2D - EEG_mean1;
% 
% EEG_mean2 = mean(simu_EEG2_2D,2);
% EEG_mean2 = repmat(EEG_mean2,1,len2*trial2);
% simu_EEG2 = simu_EEG2_2D - EEG_mean2;
for i = 1:trial
    simu_EEG(:,:,i) = detrend(squeeze(simu_EEG(:,:,i))',0)';
  %  simu_EEG2(:,:,i) = detrend(squeeze(simu_EEG2(:,:,i))','constant')';
end
%%%%%%%% make template for SCA
X = mean(simu_EEG,3);
[X_v,X_eig] = eig(X*X');
X_eig = diag(X_eig);

I_eig = cumsum(X_eig) > sum(X_eig)*10^-5;

X_v = X_v(:,I_eig);
X = X_v'*X;
% noise_to_keep_pos = 10^-6*randn(60,512);
% X = X+noise_to_keep_pos;

sqreptempl = ERP_template;
extr_method = 1;
transform_type = 1;
s_est = zeros(nc,len);
w = zeros(sum(I_eig),nc);
for i = 1:nc
	[s_est(i,:),w(:,i)] = extr_sparse(X, squeeze(sqreptempl(i,:)), extr_method,transform_type);
end
w = (w'*X_v')';
%%%%% normalize spatial filter
ampw = sqrt(diag(w'*w)).^-1;
w = w*diag(ampw);
s_est = diag(ampw)*s_est;
%%%%%%%%%%

cofmat = corrcoef(s_est,ERP_template);
if cofmat(1,2)>0
    neg_flag = -1;
else 
    neg_flag = 1;
end
%s_est = reshape(s_est,[nc,len,trial]);
Z = zeros(nc,len);
Z(1:nc,:) = neg_flag*s_est;
A = zeros(chan,nc);
A(:,1:nc) = neg_flag*pinv(w');
SF = neg_flag*w';
nc_est = nc;
end

