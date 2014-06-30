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
%clear all
clc
%%%%% add eeglab for run SOBI
% eeglab
% close all
%%%%% add C-SOBI function and template
addpath ~/Documents/STEP/functions
load sq_template1.mat
%%%%% add SCA function
% addpath /home/wch/Documents/ERP_freq_max/extr_onesrc2

SNR = [-20:5:5]; num_SNR = length(SNR);
chan = 30;
len = 500;
trial = 100;
fs = 1000;
flag_fig = false;
MC_times = 100;
frq_up = 100;
nc_up = 10;
nc = 3;
%lambda = 4;
K = 2*floor(frq_up/(fs/len)); 
%res_pct = 90;

A = zeros(num_SNR,MC_times,chan,nc_up);
Z = zeros(num_SNR,MC_times,nc_up,len);
nc_est = zeros(num_SNR,MC_times);
C = zeros(num_SNR,MC_times,nc_up,K);
Phi = zeros(num_SNR,MC_times,chan,chan);
SF = zeros(num_SNR,MC_times,nc_up,chan);

A_or = zeros(num_SNR,MC_times,chan,nc);
Z_or = zeros(num_SNR,MC_times,nc,len);
simu_EEG = zeros(num_SNR,MC_times,chan,len,trial);
SF_or = zeros(num_SNR,MC_times,nc,chan);

A_ica1 = zeros(num_SNR,MC_times,chan,chan);
Z_ica1 = zeros(num_SNR,MC_times,chan,len);
nc_est_ica1 = zeros(num_SNR,MC_times);

A_ica2 = zeros(num_SNR,MC_times,chan,chan);
Z_ica2 = zeros(num_SNR,MC_times,chan,len);
nc_est_ica2 = zeros(num_SNR,MC_times);

A_ica3 = zeros(num_SNR,MC_times,chan,chan);
Z_ica3 = zeros(num_SNR,MC_times,chan,len);
nc_est_ica3 = zeros(num_SNR,MC_times);

A_csobi = zeros(num_SNR,MC_times,chan,chan);
Z_csobi = zeros(num_SNR,MC_times,chan,len);
nc_est_csobi = zeros(num_SNR,MC_times);

A_sca = zeros(num_SNR,MC_times,chan,chan);
Z_sca = zeros(num_SNR,MC_times,chan,len);
nc_est_sca = zeros(num_SNR,MC_times);

amari_sp_cca = zeros(num_SNR,MC_times);
amari_tw_cca = zeros(num_SNR,MC_times);
snr_tw_cca = zeros(num_SNR,MC_times,nc);


amari_sp_ica1 = zeros(num_SNR,MC_times);
amari_tw_ica1 = zeros(num_SNR,MC_times);
snr_tw_ica1 = zeros(num_SNR,MC_times,nc);

amari_sp_ica2 = zeros(num_SNR,MC_times);
amari_tw_ica2 = zeros(num_SNR,MC_times);
snr_tw_ica2 = zeros(num_SNR,MC_times,nc);

amari_sp_ica3 = zeros(num_SNR,MC_times);
amari_tw_ica3 = zeros(num_SNR,MC_times);
snr_tw_ica3 = zeros(num_SNR,MC_times,nc);

amari_sp_csobi = zeros(num_SNR,MC_times);
amari_tw_csobi = zeros(num_SNR,MC_times);
snr_tw_csobi = zeros(num_SNR,MC_times,nc);

amari_sp_sca = zeros(num_SNR,MC_times);
amari_tw_sca = zeros(num_SNR,MC_times);
snr_tw_sca = zeros(num_SNR,MC_times,nc);

rng(44,'twister');
f_alpha_seed = randperm(num_SNR*MC_times*chan);
f_alpha_seed = reshape(f_alpha_seed,[num_SNR,MC_times,chan]);

ERP_template = ERP_template;
matlabpool open
for iSNR = 1:num_SNR
   % matlabpool open
    parfor iMC = 1:MC_times
      %  ERP_template = ERP_template;
	   [simu_EEG(iSNR,iMC,:,:,:),A_or(iSNR,iMC,:,:),SF_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:)]= sERP_data_gen(SNR(iSNR),iMC,chan,len,trial,fs,flag_fig,squeeze(f_alpha_seed(iSNR,iMC,:)));
        
       
       [A(iSNR,iMC,:,:),SF(iSNR,iMC,:,:),C(iSNR,iMC,:,:),Z(iSNR,iMC,:,:),nc_est(iSNR,iMC),Phi(iSNR,iMC,:,:)] = ERP_CCA_scond(squeeze(simu_EEG(iSNR,iMC,:,:,:)),frq_up,nc_up,fs); %%% Z is waveform; C is coefficients; A is spatial pattern
        
       [amari_sp_cca(iSNR,iMC),amari_tw_cca(iSNR,iMC),snr_tw_cca(iSNR,iMC,:)] = algo_eval(A(iSNR,iMC,:,:),Z(iSNR,iMC,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:),nc_est(iSNR,iMC),nc);

%        [A_ica1(iSNR,iMC,:,:),Z_ica1(iSNR,iMC,:,:),nc_est_ica1(iSNR,iMC)] = ERP_PLSOBI(squeeze(simu_EEG(iSNR,iMC,:,:,:)),squeeze(Z_or(iSNR,iMC,:,:)),0);
        
%        [amari_sp_ica1(iSNR,iMC),amari_tw_ica1(iSNR,iMC),snr_tw_ica1(iSNR,iMC,:)] = algo_eval(A_ica1(iSNR,iMC,:,:),Z_ica1(iSNR,iMC,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:),nc_est_ica1(iSNR,iMC),nc);

%        [A_ica2(iSNR,iMC,:,:),Z_ica2(iSNR,iMC,:,:),nc_est_ica2(iSNR,iMC)] = ERP_PLSOBI(squeeze(simu_EEG(iSNR,iMC,:,:,:)),squeeze(Z_or(iSNR,iMC,:,:)),0.8);
         
%        [amari_sp_ica2(iSNR,iMC),amari_tw_ica2(iSNR,iMC),snr_tw_ica2(iSNR,iMC,:)] = algo_eval(A_ica2(iSNR,iMC,:,:),Z_ica2(iSNR,iMC,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:),nc_est_ica2(iSNR,iMC),nc);
        
%        [A_ica3(iSNR,iMC,:,:),Z_ica3(iSNR,iMC,:,:),nc_est_ica3(iSNR,iMC)] = ERP_PLSOBI(squeeze(simu_EEG(iSNR,iMC,:,:,:)),squeeze(Z_or(iSNR,iMC,:,:)),1);
        
%        [amari_sp_ica3(iSNR,iMC),amari_tw_ica3(iSNR,iMC),snr_tw_ica3(iSNR,iMC,:)] = algo_eval(A_ica3(iSNR,iMC,:,:),Z_ica3(iSNR,iMC,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:),nc_est_ica3(iSNR,iMC),nc);
        

%        [A_csobi(iSNR,iMC,:,:),Z_csobi(iSNR,iMC,:,:),nc_est_csobi(iSNR,iMC)] = ERP_cSOBI_tmp(squeeze(simu_EEG(iSNR,iMC,:,:,:)),ERP_template);

%        [amari_sp_csobi(iSNR,iMC),amari_tw_csobi(iSNR,iMC),snr_tw_csobi(iSNR,iMC,:)] = algo_eval(A_csobi(iSNR,iMC,:,:),Z_csobi(iSNR,iMC,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:),nc_est_csobi(iSNR,iMC),nc);
        
%        [A_sca(iSNR,iMC,:,:),Z_sca(iSNR,iMC,:,:),nc_est_sca(iSNR,iMC)] = ERP_SCA(squeeze(simu_EEG(iSNR,iMC,:,:,:)),ERP_template);
        
%        [amari_sp_sca(iSNR,iMC),amari_tw_sca(iSNR,iMC),snr_tw_sca(iSNR,iMC,:)] = algo_eval(A_sca(iSNR,iMC,:,:),Z_sca(iSNR,iMC,:,:),A_or(iSNR,iMC,:,:),Z_or(iSNR,iMC,:,:),nc_est_sca(iSNR,iMC),nc);
       
    end
  %  matlabpool close
    save(['result_SNR_',num2str(SNR(iSNR)),'cca.mat'],'simu_EEG','A_or','Z_or','A','C','Z','nc_est','Phi',...
                                   'amari_sp_cca','amari_tw_cca','snr_tw_cca');
end
matlabpool close  
    %,...                                  % 'A_ica1','Z_ica1','nc_est_ica1','amari_sp_ica1','amari_tw_ica1','snr_tw_ica1',...
                                   %'A_ica2','Z_ica2','nc_est_ica2','amari_sp_ica2','amari_tw_ica2','snr_tw_ica2',...
                                   %'A_ica3','Z_ica3','nc_est_ica3','amari_sp_ica3','amari_tw_ica3','snr_tw_ica3');%,...
                                   %'A_sca','Z_sca','nc_est_sca','amari_sp_sca','amari_tw_sca','snr_tw_sca');   




end


%[amari_sp,amari_tw,mse_tw] = algo_eval(A,C,Z,C_or,Z_or);

