function [A,C,Z,nc_est] = ERP_CCA_sparsev4_2(simu_EEG,frq_up,nc_up,fs,lambda)
%ERP_CCA_sparse - ERP extraction with sparse CCA
%
%
% Syntax:  [A,C] = ERP_CCA_sparse(simu_EEG,frq_up,nc_up,fs)
%
% Inputs:
%    simuEEG - EEG data, channel*length*trial
%    frq_up - frequency uplimit (defult 100Hz)
%    nc_up - number of ERP components uplimit (defult 10)
%    fs - sample rate (defult 1000)
%    lambda - parameter to control sparseness
%
% Outputs:
%    A - spatial pattern of the ERP components chan*nc_up
%    C - coefficient matrix nc_up*2K
%    nc_est - number of components estimated by this method
% Example: 
%    [A,C,nc_est] = ERP_CCA_sparse(simu_EEG,100,10,1000,4)
%
% Other m-files required: ERP_date_gen.m f_alpha_gaussian.m
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Chaohua Wu
%Department of Biomedical Engineering, Tsinghua University
% email: xs.wuchaohua@outlook.com
%
% Jan 29 2013; output number of components estimated and the ERP waveform 
% Jan 14 2013; update to 4.2, estimating C using Empirical Bayesian or Lasso embeded in MATLAB
% Dec 19 2013; update to version 4, use new form
% Dec 14 2013; update A with contrain ||a_i|| = 1 
% Nov 11 2013; update A without constrain but keep its norm to be 1
%
% Nov 8 2013; using least square with constrains; 
%             lambda is calculated with the method in "Blind Source Separation by Sparse Decomposition in a Signal
% 			  Dictionary "
% Nov 6 2013; update component by component; 
%             size of C and A are not modified, the elements less than 10^-6 are set to zeros
%             update C (l1 norm contrain optimization) using a new toolbox. Hope it works better!!
% Oct 30 2013; Last revision: 

%------------- BEGIN CODE --------------
% if nargin < 1
% 	error('EEG data matrix is required!!')
% end
% cd /home/chaohua/Documents/CCA_for_ERP/cvx
% cvx_setup
% cd ../

%addpath /home/chaohua/Documents/CCA_for_ERP/package/utilities

if nargin == 0
    clc
    clear all
    close all
	[simu_EEG,sp,twave]= ERP_data_gen(-15);
	frq_up = 100;
	nc_up = 10;
	fs = 1000;
    lambda = 4;
end

if nargin == 1
	frq_up = 100;
	nc_up = 10;
	fs = 1000;
    lambda = 4;
end

[chan,len,trial] = size(simu_EEG);
nc = nc_up; % number of components
K = floor(frq_up/(fs/len)); % number of frequency

concatEEG = simu_EEG(:,:);
base = zeros(2*K,len*trial);
for i = 1:K
	base(i,:) = sin(2*pi*i/len*[0:len*trial-1]);
	base(i+K,:) = cos(2*pi*i/len*[0:len*trial-1]);
end 

%%%% base norm = 1
base_norm = diag(base*base');
base = diag(1./sqrt(base_norm))*base;

%%%%% initialize with PCA

EEG_mean = mean(simu_EEG,3);
eta = concatEEG - repmat(EEG_mean,1,trial);
Phi = eta*eta'/(len*trial);
X_tilde = Phi^(-0.5)*concatEEG;

[A_tilde,S,V] = svd(X_tilde*base','econ');
A_tilde = A_tilde(:,1:nc);
A = Phi^0.5*A_tilde;
C = S(1:nc,1:nc)*V(:,1:nc)';

i = 1;dA = 0;

nc_flag = true(nc,1);
%matlabpool open
while ((i< 30)||(dA>0.001))&&(i< 50)

%for i = 1:250
%%%% update C
[C,nc_flag] = updateC(concatEEG,A,C,Phi,base,nc_flag,lambda);

%%%% update A
A_pre = A;
A = updateA(concatEEG,A,C,Phi,base);

%%%% updata Phi
eta = concatEEG - A*C*base;
Phi = eta*eta'/(len*trial);

i = i+1;

% L = 0.5*log(det(Phi))*trial*len+lambda*sum(sum(abs(C)))+0.5*trace(A'*A);
dA = norm(A_pre-A,'fro')./norm(A_pre,'fro');
% disp(['iter = ',num2str(i),'L = ',num2str(L),'\Delta A = ',num2str(dA)]);
end

%matlabpool close
I = sum(abs(A),1) > 10^-5;
%A = A(:,I);
%C = C(I,:);
nc_est  = sum(I);
% I = sum(abs(C),2) > 10^-8;
% C = C(I,:);

% norm(sp*twave-A*C*base(:,1:len),'fro')/norm(sp*twave,'fro')
[U, S, V] = svd(A*C*base(:,1:len),0);
A = zeros(chan,nc);
A(:,1:nc_est) = U(:,1:nc_est);
Z = zeros(nc,len);
Zt = S*V';
Z(1:nc_est,:) = Zt(1:nc_est,:);
% ai = amari(A_or,sp,3);
% wave = (A\A_or)^-1*C*base;
% wave = wave(:,1:len);
end

function [C,nc_flag] = updateC(concatEEG,A,C,Phi,base,nc_flag,lambda)

[chan, nc] = size(A); K = size(C,2);
%nc_flag = true(1,nc);
X_tilde = Phi^(-0.5)*concatEEG*base';
A_tilde = Phi^(-0.5)*A;

%W_tilde = pinv(A_tilde);
%C_value = zeros(nc,2*K);
%X_temp = W_tilde*X_tilde;

% options.RegX = lambda;
% options.RegType = 1;
% C = LS(A_tilde', X_tilde', options);

%%%%%%%%%%%%% lasso provided by MATLAB
C_temp = C; 
for ik = 1:K
    C_temp(:,ik)= lasso(A_tilde,X_tilde(:,ik),'Lambda',lambda);
end

C=C_temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% LS with EB
% C_temp = C; 
% parfor ik = 1:K
%     C_temp(:,ik)= EB_LS(X_tilde(:,ik),A_tilde,10^-5);
%  end
% 
% C=C_temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for inc = 1:nc
%   	if nc_flag(inc)
% 	 	if norm(C(inc,:),'fro') < 10^-6
%             nc_flag(inc) = false;
%             C(inc,:) = 0;
%         end
%  	end
% end

end


function A = updateA(concatEEG,A,C,Phi,base)

[chan,nc] = size(A);
X_tilde = Phi^(-0.5)*concatEEG*base';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.RegX = 1;
options.RegType = 2;
A_tilde = LS(C, X_tilde, options)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = Phi^0.5*A_tilde;

% nc_nz = sum(nc_flag);
% A_temp = zeros(chan,nc);
% [U,D] = eig(Phi);
% for ic = 1:nc
%     if nc_flag(ic) ~= 0
%         x_bar = concatEEG*base'-A*C+A(:,ic)*C(ic,:);
%         A_temp(:,ic) = U*(D+C(ic,:)*C(ic,:)'*eye(chan))^-1*U'*x_bar*C(ic,:)';
%     end
%     
% end
% 
% A = A_temp;


% [U1,D1] = eig(Phi);
% [V1,D2] = eig(C*C');
% 
% X_hat = U1'*concatEEG*base'*C'*V1;
% 
% cof_mat = repmat(diag(D1),1,nc)+repmat(diag(D2)',chan,1);
% 
% A_temp = X_hat./cof_mat;
% 
% A = U1*A_temp*V1';
end
%------------- END OF CODE --------------

