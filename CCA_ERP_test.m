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
close all
clear all
clc
addpath /home/chaohua/Documents/CCA_for_ERP/package/utilities
eeglab

SNR = [-15:5:5]; num_SNR = length(SNR);
chan = 30;
len = 500;
trial = 100;
fs = 1000;
flag_fig = false;
MC_times = 100;
frq_up = 100;
nc_up = 10;
nc = 3;
lambda = 4;
K = 2*floor(frq_up/(fs/len)); 
res_pct = 90;

A = zeros(num_SNR,MC_times,chan,nc_up);
Z = zeros(num_SNR,MC_times,nc_up,len);
nc_est = zeros(num_SNR,MC_times);
C = zeros(num_SNR,MC_times,nc_up,K);

A_or = zeros(num_SNR,MC_times,chan,3);
Z_or = zeros(num_SNR,MC_times,3,len);
simu_EEG = zeros(num_SNR,MC_times,chan,len,trial);

A_ica1 = zeros(num_SNR,MC_times,chan,chan);
Z_ica1 = zeros(num_SNR,MC_times,chan,len);
nc_est_ica1 = zeros(num_SNR,MC_times);

A_ica2 = zeros(num_SNR,MC_times,chan,chan);
Z_ica2 = zeros(num_SNR,MC_times,chan,len);
nc_est_ica2 = zeros(num_SNR,MC_times);


A_pca = zeros(num_SNR,MC_times,chan,chan);
Z_pca = zeros(num_SNR,MC_times,chan,len);
nc_est_pca = zeros(num_SNR,MC_times);

amari_sp_cca = zeros(num_SNR,MC_times);
amari_tw_cca = zeros(num_SNR,MC_times);
rse_tw_cca = zeros(num_SNR,MC_times);


amari_sp_ica1 = zeros(num_SNR,MC_times);
amari_tw_ica1 = zeros(num_SNR,MC_times);
rse_tw_ica1 = zeros(num_SNR,MC_times);

amari_sp_ica2 = zeros(num_SNR,MC_times);
amari_tw_ica2 = zeros(num_SNR,MC_times);
rse_tw_ica2 = zeros(num_SNR,MC_times);

amari_sp_pca = zeros(num_SNR,MC_times);
amari_tw_pca = zeros(num_SNR,MC_times);
rse_tw_pca = zeros(num_SNR,MC_times);
matlabpool open
for iSNR = 1:num_SNR
	parfor iMC = 1:MC_times
		[simu_EEG(iSNR,iMC,:,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:)]= ERP_data_gen(SNR(iSNR),iMC,chan,len,trial,fs,flag_fig);
        
        [A(iSNR,iMC,:,:),C(iSNR,iMC,:,:),Z(iSNR,iMC,:,:),nc_est(iSNR,iMC)] = ERP_CCA_sparsev4_2(squeeze(simu_EEG(iSNR,iMC,:,:,:)),frq_up,nc_up,fs,lambda); %%% Z is waveform; C is coefficients; A is spatial pattern
        
        [amari_sp_cca(iSNR,iMC),amari_tw_cca(iSNR,iMC),rse_tw_cca(iSNR,iMC)] = algo_eval(A(iSNR,iMC,:,:),Z(iSNR,iMC,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:),nc_est(iSNR,iMC),nc);

        [A_ica1(iSNR,iMC,:,:),Z_ica1(iSNR,iMC,:,:),nc_est_ica1(iSNR,iMC)] = ERP_ICA1(squeeze(simu_EEG(iSNR,iMC,:,:,:)),res_pct);
        
        [amari_sp_ica1(iSNR,iMC),amari_tw_ica1(iSNR,iMC),rse_tw_ica1(iSNR,iMC)] = algo_eval(A_ica1(iSNR,iMC,:,:),Z_ica1(iSNR,iMC,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:),nc_est(iSNR,iMC),nc);

        [A_ica2(iSNR,iMC,:,:),Z_ica2(iSNR,iMC,:,:),nc_est_ica2(iSNR,iMC)] = ERP_ICA2(squeeze(simu_EEG(iSNR,iMC,:,:,:)),res_pct);
         
        [amari_sp_ica2(iSNR,iMC),amari_tw_ica2(iSNR,iMC),rse_tw_ica2(iSNR,iMC)] = algo_eval(A_ica2(iSNR,iMC,:,:),Z_ica2(iSNR,iMC,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:),nc_est(iSNR,iMC),nc);

        [A_pca(iSNR,iMC,:,:),Z_pca(iSNR,iMC,:,:),nc_est_pca(iSNR,iMC)] = ERP_PCA(squeeze(simu_EEG(iSNR,iMC,:,:,:)),res_pct);

        [amari_sp_pca(iSNR,iMC),amari_tw_pca(iSNR,iMC),rse_tw_pca(iSNR,iMC)] = algo_eval(A_pca(iSNR,iMC,:,:),Z_pca(iSNR,iMC,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:),nc_est(iSNR,iMC),nc);
        
        
       
    end
    save(['result_SNR_',num2str(iSNR),'.mat'],'simu_EEG','A_or','Z_or','A','C','Z','nc_est','amari_sp_cca','amari_tw_cca','rse_tw_cca',...
                                   'A_ica1','Z_ica1','nc_est_ica1','amari_sp_ica1','amari_tw_ica1','rse_tw_ica1',...
                                   'A_ica2','Z_ica2','nc_est_ica2','amari_sp_ica2','amari_tw_ica2','rse_tw_ica2',...
                                   'A_pca','Z_pca','nc_est_pca','amari_sp_pca','amari_tw_pca','rse_tw_pca');   
end

matlabpool close

end


%[amari_sp,amari_tw,mse_tw] = algo_eval(A,C,Z,C_or,Z_or);

