function [SP,SF,C,Z,nc_est,Phi] = ERP_CCA_scond(simu_EEG,frq_up,nc_up,fs)
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
%    lambda - parameter to control sparseness; a vectot containing all possible parameters
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
% Department of Biomedical Engineering, Tsinghua University
% email: xs.wuchaohua@outlook.com

% Jun 10 2014; setting parameter via cross validation
% May 19 2014; setting parameter in lasso using Generalized Cross Validation; function update C is modified
% Jan 29 2014; output number of components estimated and the ERP waveform 
% Jan 14 2014; update to 4.2, estimating C using Empirical Bayesian or Lasso embeded in MATLAB
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
startup

if nargin == 0
    clc
    clear all
    close all
	[simu_EEG,SP,SF,Z_or] = sERP_data_gen(-20,1);
	frq_up = 100;
	nc_up = 10;
	fs = 1000;
%    lambda = [10^-2, 10^-1, 1,4,8,16,40,10^2 ];
end

if nargin == 1
	frq_up = 100;
	nc_up = 10;
	fs = 1000;
%    lambda = [10^-2, 10^-1, 1,4,8,16,40,10^2 ];
end

[chan,len,trial] = size(simu_EEG);
%nc = nc_up; % number of components
K = floor(frq_up/(fs/len)); % number of frequency


base = zeros(2*K,len);
for i = 1:K
	base(i,:) = sin(2*pi*i/len*[0:len-1]);
	base(i+K,:) = cos(2*pi*i/len*[0:len-1]);
end 

% %%%% base norm = 1
% base_norm = diag(base*base');
% base = diag(1./sqrt(base_norm))*base;

%%%%% remove base line

for i = 1:trial
    simu_EEG(:,:,i) = detrend(squeeze(simu_EEG(:,:,i))','constant')';
 
end
%%%% initialize Phi, covariance matrix
concatEEG = simu_EEG(:,:);
EEG_mean = mean(simu_EEG,3);
eta = concatEEG - repmat(EEG_mean,1,trial);
Phi = eta*eta'/(len*trial);

AC_old = zeros(chan,2*K);
n_fold = 5;
delta_AC = 1;
i = 0;
while (delta_AC > 10^-6)&& (i<20)

i = i+1;
%%%% update A C

[AC,mu_cv,num_comp]= updateAC(simu_EEG,base,Phi,n_fold);
if i > 1
    delta_AC = norm(AC-AC_old,'fro')/norm(AC_old,'fro');
end
disp(['update rate: ',num2str(delta_AC),' \mu after cross-validation',num2str(mu_cv)]);
AC_old = AC;

%%%% update Phi
base_all = repmat(base,1,trial);
%base_all = base_all/norm(base_all(1,:),'fro');
Phi = (concatEEG-AC*base_all)*(concatEEG-AC*base_all)'/(len*trial);


end

nc_est = num_comp;
SF = zeros(nc_up,chan);
SP = zeros(chan,nc_up);
C = zeros(nc_up,2*K);
Z = zeros(nc_up,len);

[SF_inv,D,C_temp] = svds(Phi^-0.5*AC,nc_est);
SF(1:nc_est,:) = pinv(SF_inv);
SP(:,1:nc_est) = Phi^0.5*SF_inv;
C(1:nc_est,:) = D*C_temp';
Z(1:nc_est,:) = C(1:nc_est,:)*base;

end


function [AC,mu_cv,num_comp]= updateAC(simu_EEG,base,Phi,n_fold)
%%%%% solve the problem ||X-B*\Phi||_F^2+\mu*||B||_*
%%%%% using cross-validation to determine the right parameter \mu
Phi_nsqrt = Phi^-0.5;
Phi_sqrt = Phi^ 0.5;


[chan,len,trial] = size(simu_EEG);
K = size(base,1);
% A_temp = zeros(chan,nc_up,n_fold);
% C_temp = zeros(nc_up,K,n_fold);
mu_test = 10.^[-2:0.2:-0.4];
n_mu = length(mu_test);
error_cv = zeros(n_mu,n_fold);

trial_index = randperm(trial);
trial_perfold = floor(trial/n_fold);
if mod(trial,n_fold) ~= 0
    disp('warning: some data will be ignored in cross-validation');
end
trial_index = trial_index(1:n_fold*trial_perfold);
trial_index  = reshape(trial_index,n_fold,trial_perfold);

base_train = repmat(base,1,(n_fold-1)*trial_perfold);%% how to normalize?
norm_base = norm(base_train(1,:),'fro');
base_train = base_train/norm_base;
base_test = repmat(base,1,trial_perfold);
%base_test = base_test/norm_base;


for ifold = 1:n_fold
    fold_test = false(n_fold,1);
    fold_test(ifold) = true;
    fold_train = ~fold_test;
    train_index = trial_index(fold_train,:);
    train_index = train_index(:);
    test_index = trial_index(fold_test,:);
    test_index = test_index(:);

    X_tilde_train = simu_EEG(:,:,train_index);
    X_tilde_train = Phi_nsqrt*X_tilde_train(:,:);
    
    X_tilde_test = simu_EEG(:,:,test_index);
    X_tilde_test = Phi_nsqrt*X_tilde_test(:,:);
    %%%%%%%%%%%%%%%% re think about it
    %[U,S,V] = mexsvd(base_train',0); 
    G = sparse(norm_base*eye(min(size(base_train)))); 
    H = base_train*X_tilde_train'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bb = H(:); 
    Amap  = @(X) Amap_MLR(X,G);
    ATmap = @(y) ATmap_MLR(y,G);
    B = ATmap(bb); 

    
    options.tol = 1e-8; 
    mumax = svds(sparse(B),1,'L',options);
    par.tol     = 1e-4;
    par.verbose = 0;
    par.plotyes = 0;
    par.truncation       = 1; 
    par.truncation_gap   = 20; 
    par.maxiter  = 200; 
    problem_type = 'NNLS'; 

    for imu = 1: n_mu
        %M.U = V'*M.U; 

        %%
       % tstart = clock;

        mu_scaling = mu_test(imu);  
        mutarget   = mu_scaling*mumax;

        par.continuation_scaling = mu_scaling;  
     
        [AC_tilde,iter,time,sd,hist] = ...
        APGL(K,chan,problem_type,Amap,ATmap,bb,mutarget,0,par);

         %[U_ac,D_ac,V_ac] = svd(AC_tilde,'econ'); 
%         I_kpt = diag(D_ac) < 10^-4;
%         U_ac = U_ac(:,I_kpt);
%         D_ac = D_ac(I_kpt,I_kpt);
%         V_ac = V_ac(:,I_kpt);
%         AC_tilde = U_ac*D_ac*V_ac;

        error_cv(imu,ifold) = ...
        0.5*norm(X_tilde_test'-base_test'*AC_tilde,'fro')^2;% + mutarget*sum(sd);
    end

end

error_cv_mean = mean(error_cv,2);
error_cv_std = std(error_cv,0,2)/sqrt(n_fold);
[Y,I] = min(error_cv_mean);
I_d = error_cv_mean <= (Y+error_cv_std(I));

mu_cv = max(mu_test(I_d));


% error_cv_mean = mean(error_cv,2);
% error_cv_std = std(error_cv,0,2)/sqrt(n_fold);
% [Y,I] = min(error_cv_mean);
% mu_cv = mu_test(I);

%figure
%errorbar(log10(mu_test),error_cv_mean,error_cv_std)
% mu_cv = 10^-1;
%X_tilde_train = simu_EEG(:,:,train_index);
X_tilde = Phi_nsqrt*simu_EEG(:,:);
base_all = repmat(base,1,trial);
norm_baseall = norm(base_all(1,:),'fro');
base_all = base_all/norm_baseall;
    
G = sparse(norm_baseall*eye(min(size(base_all)))); 
H = base_all*X_tilde'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bb = H(:); 
Amap  = @(X) Amap_MLR(X,G);
ATmap = @(y) ATmap_MLR(y,G);
B = ATmap(bb); 

options.tol = 1e-8; 
mumax = svds(sparse(B),1,'L',options);
par.tol     = 1e-4;
par.verbose = 0;
par.plotyes = 0;
par.truncation       = 1; 
par.truncation_gap   = 20; 
par.maxiter  = 200; 
problem_type = 'NNLS'; 
mu_scaling = mu_cv;  
mutarget   = mu_scaling*mumax;
par.continuation_scaling = mu_scaling;  

[AC_tilde,iter,time,sd,hist] = ...
APGL(K,chan,problem_type,Amap,ATmap,bb,mutarget,0,par);


[U_ac,D_ac,V_ac] = svd(AC_tilde,'econ');
%trac = trace(D_ac);
I_kpt = diag(D_ac) > 10^-3;
num_comp = sum(I_kpt);
U_ac = U_ac(:,I_kpt);
D_ac = D_ac(I_kpt,I_kpt);
V_ac = V_ac(:,I_kpt);
AC_tilde = U_ac*D_ac*V_ac';

AC = Phi_sqrt*AC_tilde';


end
%------------- END OF CODE --------------

