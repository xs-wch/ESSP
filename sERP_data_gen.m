function [simu_EEG,SP,SF,Z]= sERP_data_gen(SNR,iMC,chan,len,trial,fs,flag_fig,f_alpha_seed)
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
trial = 50;
fs = 1000;
flag_fig = true;
iMC = 5;
f_alpha_seed = [1:30];
end

if nargin == 1
chan = 30;
len = 500;
trial = 100;
fs = 1000;
flag_fig = false;
iMC = 1;
f_alpha_seed = [1:30];
end


if nargin == 2
chan = 30;
len = 500;
trial = 100;
fs = 1000;
flag_fig = false;
f_alpha_seed = [1:30];
end

SNR_white = 20;
nc = 3;
%dnc = 3;
%%%% ERP signal generation %%%%
Z =10*[gampdf([0:len-1],15,5);gampdf([0:len-1],27,6);gampdf([0:len-1],38,8)];
% if flag_fig 
%     figure
%     for i = 1:nc
%         subplot(2,2,i)
%         plot(Z(i,:));
%     end
% end


mask =ones(1,len);
fbound = 60;
mask_index = floor(fbound/(fs/len));
mask(mask_index+1:len-mask_index-1) = 0;
for i = 1:nc
    tempfft = fft(Z(i,:));
    tempfft = tempfft.*mask;
    Z(i,:) = real(ifft(tempfft));
end
Z = detrend(Z','constant')';
[Ut,St,Vt] = svd(Z,'econ');
Z = 10*diag([-1,1,1])*St*Vt';
if flag_fig 
    figure
    for i = 1:nc
        subplot(2,2,i)
        plot(Z(i,:));
    end
end


%  ERP_template = zeros(nc,len);
% % ERP_template.Z2 = zeros(nc,len);
%  for i = 1:nc
%      ERP_template(i,:) = Z(i,:)>0.1;
% %     ERP_template.Z2(i,:) = Z2(i,:)>0.05;
% % % Z1_template(2,:) = Z1(2,:)>0.05;Z2_template(2,:) = Z2(2,:)>0.05;
% % % Z1_template(3,:) = Z1(3,:)>0.05;Z2_template(3,:) = Z2(3,:)>0.05;
% % % Z1_template(4,:) = Z1(4,:)>0.05;Z2_template(4,:) = Z2(4,:)>0.05;
%  end
% if flag_fig
%     figure
%     for i = 1:nc
%         subplot(2,2,i)
%         %plot(ERP_template(i,:),'r--');hold on;
%         plot(ERP_template(i,:));
%     end
% end
% 
% save('sq_template.mat','ERP_template');

%%%%%%%%%%%%%%%% spatial filter
rng((564+SNR)*iMC,'twister'); 
SF = randn(chan,nc); 
SF = orth(SF);
  
%%%% spontaneous EEG
nspEEG = chan;
rng((421+SNR)*iMC,'twister'); 
spEEGmix =  randn(chan,nspEEG);

spEEGsource = zeros(nspEEG,len*trial); 

for insp = 1:nspEEG
    
    temp = f_alpha_gaussian(len*trial,f_alpha_seed(insp));
    spEEGsource(insp,:) = detrend(temp','constant')';

end
%%%%%%%%%% mixed spontaneous EEG
spEEG = spEEGmix*spEEGsource;
spEEG_cov = spEEG*spEEG'/(len*trial);

%%%%%%%%%% spatial pattern
C = spEEG_cov^0.5*SF;
if ~isreal(C)
    disp('stop')
    pause
end
SP = C;
SF = pinv(SF);
%%%%%%%%%% mixed ERP
ERP_simu = repmat(C*Z,[1 1 trial]);

%%%% mix %%%%%
power_ERP = trace(ERP_simu(:,:)*ERP_simu(:,:)');
power_spEEG = trace(spEEG*spEEG');

PR = 10^(SNR/10);% power ratio
ampN = sqrt((power_ERP)/PR/power_spEEG);
spEEG = ampN.*spEEG;

simu_EEG = ERP_simu(:,:)+spEEG(:,:);
simu_EEG = reshape(simu_EEG,[chan,len,trial]);

%%%% add white noise %%%%
rng((4517+SNR)*iMC,'twister'); 
white_noise = randn(chan,len,trial);

power_whn = trace(white_noise(:,:)*white_noise(:,:)');
power_EEG = trace(simu_EEG(:,:)*simu_EEG(:,:)');

ampN1 = sqrt((power_EEG)/(10^(SNR_white/10))/(power_whn));
white_noise = ampN1.*white_noise;

simu_EEG = simu_EEG(:,:)+white_noise(:,:);
simu_EEG = reshape(simu_EEG,[chan,len,trial]);

if flag_fig 
    figure
    plot(simu_EEG(5,:,1));

end
% if ~isreal(simu_EEG)
%     disp('stop here')
%     pause
% end
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

