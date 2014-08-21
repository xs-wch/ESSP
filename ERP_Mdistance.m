%%%%% function to evaluate N170, M-distance?

function [Pv,amplitude_diff,Y,tpr,fpr] = ERP_Mdistance(EEG1,EEG2,SF,win,figure_flag)
% EEG1: chan*len*trial 
% EEG2: chan*len*trial
% SF: spatial filter 1*chan
% win: ERP analysis time window

[chan,len,trial] = size(EEG1);
if length(SF) ~= chan
	error('Wrong spaltial filter');
end
SF1 = SF/norm(SF);
SF2 = SF1;

EEG1_comp = zeros(len,trial);
%EEG1_feature = zeros(trial,1);
EEG2_comp = zeros(len,trial);
%EEG2_feature = zeros(trial,1);
for i = 1:trial
	EEG1_comp(:,i) = SF1*squeeze(EEG1(:,:,i));
	EEG2_comp(:,i) = SF2*squeeze(EEG2(:,:,i));
end

EEG1mean = mean(EEG1_comp,2);
EEG2mean = mean(EEG2_comp,2);
[Y1,I1] = min(EEG1mean(win));
[Y2,I2] = min(EEG2mean(win));

amplitude_diff = Y1-Y2;
EEG1_feature = EEG1_comp(win(1)+I1-1,:);
EEG2_feature = EEG2_comp(win(1)+I2-1,:);

if figure_flag
    figure
    plot(mean(EEG1_comp,2),'r');
    hold on
    plot(mean(EEG2_comp,2),'b');
    legend('obj','face')
end
[H,Pv] = ttest2(EEG1_feature,EEG2_feature);

feature_min = min([EEG1_feature,EEG2_feature]);
feature_max = max([EEG1_feature,EEG2_feature]);
Y = ([EEG1_feature,EEG2_feature]-feature_min)/(feature_max-feature_min);
T = [ones(1,trial),zeros(1,trial)];
[tpr,fpr,th] = roc(T,Y);

if figure
    plotroc(T,Y);
end
end




