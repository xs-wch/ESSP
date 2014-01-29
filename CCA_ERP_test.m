function CCA_ERP_test
% CCA_ERP_test - this function is to evaluate the performance of CCA for ERP
%                 
%                
% Syntax:  CCA_ERP_test
%
% Inputs:
%
%
% Outputs:
%    
%    
%
% Example: 
%
%
% Other m-files required: f_alpha_gaussian.m ERP_CCA_sparsev4_2.m ERP_data_gen.m
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

SNR = [-20:5:5]; num_SNR = length(SNR);
chan = 30;
len = 500;
trial = 50;
fs = 1000;
flag_fig = false;
MC_times = 100;
frq_up = 100;
nc_up = 10;
lambda = 4;

for iSNR = 1:num_SNR
	for iMC = 1:MC_times
		[simu_EEG,C_or,Z_or]= ERP_data_gen(SNR(iSNR),iMC,chan,len,trial,fs,flag_fig);
        
        [A,C,Z,nc_est] = ERP_CCA_sparsev4_2(simu_EEG,frq_up,nc_up,fs,lambda) %%% Z is waveform; C is coefficients; A is spatial pattern

        [A_ica,C_ica,Z_ica,nc_est_ica] = ERP_ICA();

        [A_pca,C_pca,Z_pca,nc_est_pca] = ERP_PCA();
         
        [amari_sp,amari_tw,mse_tw] = algo_eval(A,C,Z,C_or,Z_or);


	end
end

