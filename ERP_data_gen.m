function [simu_EEG1,C,Z1]= ERP_data_gen(SNR,iMC,chan,len,trial,fs,flag_fig)
% ERP_data_gen - generate artifical EEG data with given SNR channel length trial and sample-rate
%                ERP: three gamma function 
%                spontanous EEG: 1/f noise with Gaussian distribution
% Syntax:  concatEEG = ERP_data_gen(SNR,chan,len,trial,fs)
%
% Inputs:
%    SNR - Signal to Noise Ratio of ERP
%    chan - channel number of EEG data
%    len - length of EEG in each trial
%    trial - number of ERP trials
%    fs - sample rate
%    flag_fig - true: plot the figures; false: do not plot any figure
%
% Outputs:
%    concatEEG - concatenated EEG data 
%    
%
% Example: 
%    concatEEG = ERP_data_gen(SNR,chan,len,trial,fs)
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: f_alpha_gaussian.m
% Subfunctions: none
% MAT-files required: none
%
% 

% Author: Chaohua Wu
%Department of Biomedical Engineering, Tsinghua University
% email: xs.wuchaohua@outlook.com
%
% Nov 4 2013; output spatial pattern C & components Z1
% Oct 29 2013; Last revision: 

%------------- BEGIN CODE --------------
if nargin == 0
SNR = 0; 
chan = 30;
len = 500;
trial = 100;
fs = 1000;
flag_fig = false;
end

if nargin == 1
chan = 30;
len = 500;
trial = 60;
fs = 1000;
flag_fig = false;
end

SNR_white = 20;
%iMC = iMC*randi(19);
nc = 3;
%%%% signal generation %%%%
Z1 =10*[gampdf([0:len-1],15,5);gampdf([0:len-1],45,6);gampdf([0:len-1],70,5)];
Z1 = detrend(Z1','constant')';

rng((564+SNR)*iMC,'twister'); 
C = randn(chan,nc); 

[C_temp,S_temp,V_temp] = svd(C*Z1,'econ');
C = C_temp(:,1:nc);
Z1 = S_temp(i:nc,1:nc)*V_temp(:,1:nc)';
if flag_fig
	figure;
	subplot(1,2,1)
	plot(Z1(1,:));
	xlabel('time(ms)');
	ylabel('amplitude(\muV)');
	subplot(1,2,2)
	stem([0:len-1]./len*fs,abs(fft(Z1(1,:)))./len,'filled');
    xlabel('frequency (Hz)');
    ylabel('amplitude (\muV)')
    temprepERP =  repmat(Z1,[1 trial]);

    figure;
	subplot(1,2,1)
	plot(temprepERP(1,:));
	xlabel('time(ms)');
	ylabel('amplitude(\muV)');
	subplot(1,2,2)
	stem([0:len*trial-1]./(len*trial)*fs,abs(fft(temprepERP(1,:)))./(len*trial),'filled');
    xlabel('frequency (Hz)');
    ylabel('amplitude (\muV)')
end





ERP_simu = repmat(C*Z1,[1 1 trial]);
  
%%%% spontaneous EEG
nspEEG = 30;
rng((421+SNR)*iMC,'twister'); 
spEEGmix1 =  randn(chan,nspEEG);
spEEG1 = zeros(chan,len,trial);
spEEGsource = zeros(nspEEG,len*trial); 

%for it= 1:trial
for insp = 1:nspEEG
    % rng(insp*trial+iMC,'twister');
    temp = f_alpha_gaussian(len*trial,abs(insp*(5*iMC+SNR)));
    spEEGsource(insp,:) = detrend(temp,'constant')';
end

spEEGsource = reshape(spEEGmix1*spEEGsource,[nspEEG,len,trial]);
spEEG1 = spEEGsource(:,:,randperm(trial));
%spEEG1(:,:,it) = spEEGmix1*spEEGsource(:,:,it);
%end
%temp = f_alpha_gaussian(len*trial,abs())





%%%% mix %%%%%

tempERP = C*Z1;
%tempERP = detrend(tempERP','constant')';
power_ERP1 = trial.*trace(tempERP*tempERP');
%power_spEEG1 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it= 1:trial
    power_spEEG1 = power_spEEG1 +trace(spEEG1(:,:,it)*spEEG1(:,:,it)');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
power_spEEG1 = trace(spEEG1(:,:)*spEEG1(:,:)');

ampN1 = sqrt(power_ERP1/(10^(SNR/10)*power_spEEG1));
spEEG1 = ampN1.*spEEG1;

simu_EEG1 = zeros(chan,len,trial);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it = 1:trial
    simu_EEG1(:,:,it) = ERP_simu(:,:,it)+spEEG1(:,:,it);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simu_EEG1 = ERP_simu(:,:)+spEEG1(:,:);
simu_EEG1 = reshape(simu_EEG1,[chan,len,trial]);


%%%% add white noise %%%%
rng((4517+SNR)*iMC,'twister'); 
white_noise1 = randn(chan,len,trial);

power_whn1 = trace(white_noise1(:,:)*white_noise1(:,:)');

power_EEG1 = trace(simu_EEG1(:,:)*simu_EEG1(:,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it= 1:trial
    power_EEG1 = power_EEG1 +trace(simu_EEG1(:,:,it)*simu_EEG1(:,:,it)');

    power_whn1 = power_whn1 +trace(white_noise1(:,:,it)*white_noise1(:,:,it)');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

power_EEG1
ampN1 = sqrt(power_EEG1/(10^(SNR_white/10)*power_whn1));
white_noise1 = ampN1.*white_noise1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it = 1:trial
    simu_EEG1(:,:,it) = simu_EEG1(:,:,it)+white_noise1(:,:,it);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simu_EEG1 = simu_EEG1(:,:)+white_noise1(:,:);
simu_EEG1 = reshape(simu_EEG1,[chan,len,trial]);

concatEEG = simu_EEG1(:,:);

if flag_fig
	figure
	subplot(1,2,1)
    plot(concatEEG(1,:));
	xlabel('time(ms)');
	ylabel('amplitude(\muV)');

    subplot(1,2,2)
	stem([0:len*trial-1]./(len*trial)*fs,abs(fft(concatEEG(1,:)))./(len*trial),'filled');
    xlabel('frequency (Hz)');
    ylabel('amplitude (\muV)')
 
end




%------------- END OF CODE --------------

