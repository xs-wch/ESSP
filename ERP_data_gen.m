function [simu_EEG1,simu_EEG2,SP,SF,dZ]= ERP_data_gen(SNR,iMC,chan,len,trial,fs,flag_fig);
% ERP_data_gen - generate artifical EEG data with given SNR channel length trial and sample-rate
%                ERP: three gamma function 
%                spontanous EEG: 1/f noise with Gaussian distribution
% Syntax:  concatEEG = ERP_data_gen(SNR,chan,len,trial,fs)
%
% Inputs:
%    SNR - Signal to Noise Ratio of ERP
%    iMC - index for Monto Carlo 
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
% Department of Biomedical Engineering, Tsinghua University
% email: xs.wuchaohua@outlook.com
%
% Nov 4 2013; output spatial pattern C & components Z1
% Oct 29 2013; Last revision: 

%------------- BEGIN CODE --------------
if nargin == 0
SNR = -10; 
chan = 30;
len = 500;
trial = 100;
fs = 1000;
flag_fig = true;
iMC = 5;
end

if nargin == 1
chan = 30;
len = 500;
trial = 100;
fs = 1000;
flag_fig = false;
iMC = 1;
end


if nargin == 2
chan = 30;
len = 500;
trial = 100;
fs = 1000;
flag_fig = false;
end

SNR_white = 20;
nc = 4;
dnc = 3;
%%%% ERP signal generation %%%%
Z1 =10*[gampdf([0:len-1],15,5);gampdf([0:len-1],27,6);gampdf([0:len-1],33,8);gampdf([0:len-1],38,9)];

dZ =[gampdf([0:len-1],16,5);gampdf([0:len-1],30,6);gampdf([0:len-1],40,8)];
mask =ones(1,len);
fbound = 60;
mask_index = floor(fbound/(fs/len));
mask(mask_index+1:len-mask_index-1) = 0;
for i = 1:dnc
    tempfft = fft(dZ(i,:));
    tempfft = tempfft.*mask;
    dZ(i,:) = real(ifft(tempfft));
end
dZ = detrend(dZ','constant')';
[Ut,St,Vt] = svd(dZ,'econ');
dZ = 7*diag([-1,-1,1])*St*Vt';
dZ = [dZ;zeros(1,len)];
if flag_fig 
    figure
    for i = 1:nc
        subplot(2,2,i)
        plot(dZ(i,:));
    end
end

Z2 = Z1+dZ;
if flag_fig
    figure
    for i = 1:nc
        subplot(2,2,i)
        plot(Z2(i,:),'r--');hold on;
        plot(Z1(i,:));
    end
end

% ERP_template.Z1 = zeros(nc,len);
% ERP_template.Z2 = zeros(nc,len);
% for i = 1:nc
%     ERP_template.Z1(i,:) = Z1(i,:)>0.05;
%     ERP_template.Z2(i,:) = Z2(i,:)>0.05;
% % Z1_template(2,:) = Z1(2,:)>0.05;Z2_template(2,:) = Z2(2,:)>0.05;
% % Z1_template(3,:) = Z1(3,:)>0.05;Z2_template(3,:) = Z2(3,:)>0.05;
% % Z1_template(4,:) = Z1(4,:)>0.05;Z2_template(4,:) = Z2(4,:)>0.05;
% end
% if flag_fig
%     figure
%     for i = 1:nc
%         subplot(2,2,i)
%         plot(ERP_template.Z2(i,:),'r--');hold on;
%         plot(ERP_template.Z1(i,:));
%     end
% end
% 
% save('sq_template.mat','ERP_template');

%%%%%%%%%%%%%%%% spatial filter
rng((564+SNR)*iMC,'twister'); 
SF = randn(chan,nc); 
SF = orth(SF);

% [C,S_temp,V_temp] = svds(C*Z1,nc);
% C = C_temp(:,1:nc);
% Z1 = S_temp*V_temp';
% if flag_fig
% 	figure;
% 	subplot(1,2,1)
% 	plot(Z1(1,:));
% 	xlabel('time(ms)');
% 	ylabel('amplitude(\muV)');
% 	subplot(1,2,2)
% 	stem([0:len-1]./len*fs,abs(fft(Z1(1,:)))./len,'filled');
%     xlabel('frequency (Hz)');
%     ylabel('amplitude (\muV)')
%     temprepERP =  repmat(Z1,[1 trial]);
% 
%     figure;
% 	subplot(1,2,1)
% 	plot(temprepERP(1,:));
% 	xlabel('time(ms)');
% 	ylabel('amplitude(\muV)');
% 	subplot(1,2,2)
% 	stem([0:len*trial-1]./(len*trial)*fs,abs(fft(temprepERP(1,:)))./(len*trial),'filled');
%     xlabel('frequency (Hz)');
%     ylabel('amplitude (\muV)')
% end

% ERP_simu = repmat(C*Z1,[1 1 trial]);
  
%%%% spontaneous EEG
nspEEG = 30;
rng((421+SNR)*iMC,'twister'); 
spEEGmix =  randn(chan,nspEEG);
%spEEG1 = zeros(chan,len,trial);
spEEGsource1 = zeros(nspEEG,len*trial); 
spEEGsource2 = zeros(nspEEG,len*trial); 
%for it= 1:trial
for insp = 1:nspEEG
    % rng(insp*trial+iMC,'twister');
    temp = f_alpha_gaussian(len*trial,abs(10*insp*(5*iMC+SNR)));
    spEEGsource1(insp,:) = detrend(temp','constant')';
    %spEEGsource1(insp,:) = temp;
    temp = f_alpha_gaussian(len*trial,abs(70*insp*(345*iMC+5*SNR)));
    spEEGsource2(insp,:) = detrend(temp','constant')';
    %spEEGsource2(insp,:) = temp;
end
spEEG_all = spEEGmix*[spEEGsource1,spEEGsource2];
spEEG_cov = spEEG_all*spEEG_all'/(len*trial*2);

%%%%%%%%%% spatial pattern
C = spEEG_cov^0.5*SF;
SP =C;
SF = pinv(SF);
%%%%%%%%%% mixed ERP
ERP_simu1 = repmat(C*Z1,[1 1 trial]);
ERP_simu2 = repmat(C*Z2,[1 1 trial]);
%%%%%%%%%% mixed spontaneous EEG
% spEEGsource1 = reshape(spEEGmix*spEEGsource1,[nspEEG,len,trial]);
% spEEG1 = spEEGsource1(:,:,randperm(trial));
% 
% spEEGsource2 = reshape(spEEGmix*spEEGsource2,[nspEEG,len,trial]);
% spEEG2 = spEEGsource1(:,:,randperm(trial));
spEEG1 = spEEGmix*spEEGsource1;
spEEG2 = spEEGmix*spEEGsource2;
%%%% mix %%%%%
tempdERP = C*dZ;
power_dERP = trial.*trace(tempdERP*tempdERP');
tempERP1 = C*Z1;tempERP2 = C*Z2;
power_ERP = trial.*(trace(tempERP1*tempERP1')+trace(tempERP2*tempERP2'));
power_spEEG1 = trace(spEEG1*spEEG1');
power_spEEG2 = trace(spEEG2*spEEG2');
PR = 10^(SNR/10);% power ratio
ampN = sqrt(((1/PR+1)*power_dERP-power_ERP)/(power_spEEG1+power_spEEG2));
% ampN = 0;
spEEG1 = ampN.*spEEG1;
spEEG2 = ampN.*spEEG2;

simu_EEG1 = ERP_simu1(:,:)+spEEG1(:,:);
simu_EEG1 = reshape(simu_EEG1,[chan,len,trial]);
simu_EEG2 = ERP_simu2(:,:)+spEEG2(:,:);
simu_EEG2 = reshape(simu_EEG2,[chan,len,trial]);

%%%% add white noise %%%%
rng((4517+SNR)*iMC,'twister'); 
white_noise1 = randn(chan,len,trial);

rng((456517+SNR)*iMC,'twister'); 
white_noise2 = randn(chan,len,trial);
power_whn1 = trace(white_noise1(:,:)*white_noise1(:,:)');
power_whn2 = trace(white_noise2(:,:)*white_noise2(:,:)');
power_EEG1 = trace(simu_EEG1(:,:)*simu_EEG1(:,:)');
power_EEG2 = trace(simu_EEG2(:,:)*simu_EEG2(:,:)');
ampN1 = sqrt((power_EEG1+power_EEG2)/(10^(SNR_white/10))/(power_whn1+power_whn2));
white_noise1 = ampN1.*white_noise1;
white_noise2 = ampN1.*white_noise2;
simu_EEG1 = simu_EEG1(:,:)+white_noise1(:,:);
simu_EEG1 = reshape(simu_EEG1,[chan,len,trial]);
simu_EEG2 = simu_EEG2(:,:)+white_noise2(:,:);
simu_EEG2 = reshape(simu_EEG2,[chan,len,trial]);
if flag_fig 
    figure
    plot(simu_EEG1(1,:,1));
    hold on;
    plot(simu_EEG2(1,:,1),'r--');
end

% concatEEG = simu_EEG1(:,:);
% 
% if flag_fig
% 	figure
% 	subplot(1,2,1)
%     plot(concatEEG(1,:));
% 	xlabel('time(ms)');
% 	ylabel('amplitude(\muV)');
% 
%     subplot(1,2,2)
% 	stem([0:len*trial-1]./(len*trial)*fs,abs(fft(concatEEG(1,:)))./(len*trial),'filled');
%     xlabel('frequency (Hz)');
%     ylabel('amplitude (\muV)')
%  
% end




%------------- END OF CODE --------------

