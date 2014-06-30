function [A,Z,nc_est] = ERP_PLSOBI(simu_EEG,real_ERP,lambda)
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
   [simu_EEG,SP,SF,Z_real]= sERP_data_gen(-15,43);
   real_ERP = Z;
   lambda = 0.8;
   eeglab
end
[chan,len,trial] = size(simu_EEG);
%[chan2,len2,trial2] = size(simu_EEG2);
% if (chan1~=chan2)||(len1~=len2)
% 	disp('ERP data dimension does not match')
% 	pause
% end

if (lambda>1)||(lambda<0)
	disp('lambda should be between 0 and 1');
	pause
end

[nc,len_real] = size(real_ERP);
if (len~=len_real)
	disp('data ERP length and real ERP length do not match');
	pause
end

%%%%%%% remove base line

for i = 1:trial
    simu_EEG(:,:,i) = detrend(squeeze(simu_EEG(:,:,i))','constant')';
end

%%%%%%%% make template for PLSOBI
EEG_mean = mean(simu_EEG,3);
EEG_template = repmat(EEG_mean,1,trial);

simu_EEG_PL = (1-lambda)*simu_EEG(:,:) + lambda*EEG_template;

X = simu_EEG_PL;
[winv,act] = sobi(X);
%%%%% normalize spatial filter
SF = pinv(winv);
ampw = sqrt(diag(SF*SF')).^-1;
act = diag(ampw)*act;
SF = diag(ampw)*SF;
winv = pinv(SF);
%%%%%%%%%%%%

act = reshape(act,[chan,len,trial]);
act = mean(act,3);
%%%%%%% calculate correlation coefficient matrix
cofmat = zeros(chan, nc);

for i = 1:nc
	for j = 1:chan
		tempcor = corrcoef(act(j,:)',real_ERP(i,:)');
        cofmat(j,i) = tempcor(2,1);
	end
end

[Y,I] = max(abs(cofmat));
Z = zeros(chan,len);
Z(1:nc,:) = act(I,:);
A = zeros(chan,chan);
A(:,1:nc) = winv(:,I);
nc_est = sum(sum(abs(cofmat)>0.6));

close all

