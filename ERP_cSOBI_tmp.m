function [A,Z,nc_est] = ERP_cSOBI_tmp(simu_EEG,ERP_template)
% ERP_ICA1 - extract ERP using constrained SOBI
%            the template is the averaged ERP
%            select 4 channels randomly
%            
%                
% Syntax:  [A,Z,nc_est] = ERP_cSOBI(simu_EEG1,simu_EEG2,nc)
%
% Inputs:
%          simu_EEG 1/2: EEG signal, chan*len*trial
%          nc: number of ERP components used
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

if nargin < 1
   [simu_EEG,SP,SF,Z]= sERP_data_gen(-20,1);
   load sq_template1.mat
   addpath ~/Documents/STEP/functions
end
if nargin ==1
   load sq_template.mat
end
if nargin > 2
	disp('too many parameters')
    pause
end

[chan,len,trial] = size(simu_EEG);
nc = size(ERP_template,1);
%[chan2,len2,trial2] = size(simu_EEG2);
% if (chan1~=chan2)||(len1~=len2)
% 	disp('ERP data dimension does not match')
% 	pause
% end
% if nargin < 4
% 	dataref = zeros(nc,len1*trial1+len2*trial2);
% end 



%simu_EEG_2D = simu_EEG(:,:);
%%%%%%% remove base line
for i = 1:trial
    simu_EEG(:,:,i) = detrend(squeeze(simu_EEG(:,:,i))','constant')';
 %   simu_EEG2(:,:,i) = detrend(squeeze(simu_EEG2(:,:,i))','constant')';
end

%%%%%%%% make template for CSOBI

X = simu_EEG(:,:);
dataref = repmat(ERP_template,1,trial);

y = zeros(nc,trial*len);
W = zeros(chan,nc);
for i = 1:nc
	[y(i,:),W(:,i),CPU_time,step]=c_sobi6(X,dataref(i,:),'mse');
end

%%%%% normalize Spatial Filter
ampw = sqrt(diag(W'*W)).^-1;
y = diag(ampw)*y;
W = W*diag(ampw);
%%%%%%%%%%%%%%%%

y = reshape(y,[nc,len,trial]);
Z = zeros(chan,len);
Z(1:nc,:) = mean(y,3);
A = zeros(chan,chan);
A(:,1:nc) = pinv(W');

nc_est = nc;

close all

