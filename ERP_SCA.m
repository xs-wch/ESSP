function [A,Z,nc_est] = ERP_SCA(simu_EEG,ERP_template)
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
% Jun 18 ; first version
%------------- BEGIN CODE --------------



if nargin  < 1
	%disp('at least 2 parameters required');
    addpath /home/wch/Documents/ERP_freq_max/extr_onesrc2
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
    simu_EEG(:,:,i) = detrend(squeeze(simu_EEG(:,:,i))','constant')';
  %  simu_EEG2(:,:,i) = detrend(squeeze(simu_EEG2(:,:,i))','constant')';
end
%%%%%%%% make template for SCA
X = simu_EEG(:,:);
sqreptempl = repmat(ERP_template,1,trial);
extr_method = 1;
transform_type = 2;
s_est = zeros(nc,len*trial);
w = zeros(chan,nc);
for i = 1:nc
	[s_est(i,:),w(:,i)] = extr_sparse(X, sqreptempl(i,:), extr_method,transform_type);
end

%%%%% normalize spatial filter
ampw = sqrt(diag(w'*w)).^-1;
w = w*diag(ampw);
s_est = diag(ampw)*s_est;
%%%%%%%%%%

s_est = reshape(s_est,[nc,len,trial]);
Z = zeros(chan,len);
Z(1:nc,:) = mean(s_est,3);
A = zeros(chan,chan);
A(:,1:nc) = pinv(w');
nc_est = nc;

close all

