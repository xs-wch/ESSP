%%%% this file to apply JE_wch.m on real EEG data
function real_data_SSF
clear all
close all

addpath /home/wch/Documents/ERP_freq_max/CCA_ERP_singlecond
addpath /home/wch/Documents/ERP_freq_max/CCA_ERP_singlecond/PROPACKmod
addpath /home/wch/Documents/ERP_freq_max/CCA_ERP_singlecond/solver 

load ERP_LYF_20111123p1_filterd;
srate = mycnt.header.rate;
channum = length(mycnt.electloc);
eloc = readlocs('neuros_lead68.loc');
chandelete = {'CB1','CB2','HEO','VEO','EKG','EMG','M1','M2'};%,'FP1','FP2','FPZ','AF3','AF4'};

% choose the channel of interest
chanuseful = {'T7','T8','TP7','P8','FT7','FT8','CZ','OZ'};
%[SP,SF,C,Z,nc_est,Phi] = ERP_CCA_scond_bkup(simu_EEG,frq_up,nc_up,fs)

%%
% select the useful channel as setting
[chanindex,mycnt,eloc] = chanselect(mycnt, eloc, chandelete, chanuseful);
channum = channum-length(chandelete);

ana_begin = 0.01*srate; % begining
ana_end = 0.522*srate;    % end
[target_data,face_data,face_invert_data,object_data] = dividedata(mycnt,ana_begin,ana_end);
[chan,len,trial] = size(face_data);

for i = 1:trial
    face_data(:,:,i) = detrend(squeeze(face_data(:,:,i))','constant')';
    object_data(:,:,i) = detrend(squeeze(object_data(:,:,i))','constant')';
end
wavetemplate = mean(face_data(chanindex(4),:,:),3);
ERP_template = [wavetemplate(1:250),wavetemplate(250)*ones(1,512-250)];
%[face_data, latency_column] = woody_filtering(face_data,chanindex(4));

%face_data = face_data(:,100:601,:);

[SP_face,SF_face,C_face,Z_face,nc_est_face,Phi_face] = ERP_CCA_scond_bkup(face_data,0,60,10,1000,false);

cofmat = zeros(1,nc_est_face);
for i = 1:nc_est_face
    tempcorcof = corrcoef(Z_face(i,:),ERP_template);
    cofmat(i) = tempcorcof(1,2);
end
[Y,I] = max(abs(cofmat));
%I = 1;
if cofmat(I)<0
    SF_face = -SF_face(I,:);
else
    SF_face = SF_face(I,:);
end


% close all
% h = figure
% for i = 1:nc_est
%     subplot(nc_est,2,2*i-1)
%     plot(Z(i,:));
%     subplot(nc_est,2,2*i)
%     topoplot(SP(:,i),eloc);
% end
% hgexport(h,'face_figure')



[SP_object,SF_object,C_object,Z_object,nc_est_object,Phi_object] = ERP_CCA_scond_bkup(object_data,0,60,10,1000,false);
% close all
% h = figure
% for i = 1:nc_est
%     subplot(nc_est,2,2*i-1)
%     plot(Z(i,:));
%     subplot(nc_est,2,2*i)
%     topoplot(SP(:,i),eloc);
% end
% hgexport(h,'object_figure')

cofmat = zeros(1,nc_est_object);
for i = 1:nc_est_object
    tempcorcof = corrcoef(Z_object(i,:),ERP_template);
    cofmat(i) = tempcorcof(1,2);
end
[Y,I] = max(abs(cofmat));
%I = 1;
if cofmat(I)<0
    SF_object= -SF_object(I,:);
else
    SF_object= SF_object(I,:);
end

[Pv,Y_face,T_face,tpr_face,fpr_face]= ERP_Mdistance(object_data,face_data,SF_object,SF_face,round(0.170*srate),round(0.040*srate));

auc = myauc(fpr_face,tpr_face);
end

% function [face_data, latency_column] = woody_filtering(face_data,chan_index)
% 
% [chan,len,trial] = size(face_data);
% latency_column = zeros(trial,1);
% face_average = squeeze(mean(face_data(chan_index,:,:),3))';
% data_shift = squeeze(face_data(chan_index,:,:));
% delta_faceavr = 1;
% 
% while delta_faceavr>10^-4
%     for j = 1:trial
%         tempdata = data_shift(:,j);
%         xtemp = xcorr(tempdata, face_average);
%         [Y,I] = max(xtemp(len-15:len+15));
%         deltaI = I-15-1;
%         latency_column(j) = latency_column(j)+deltaI; 
%         data_shift(:,j) = circshift(tempdata, deltaI);
%     end
%     face_average_new = mean(data_shift,2);
%     delta_faceavr = norm(face_average-face_average_new,'fro')/norm(face_average,'fro');
%     face_average = face_average_new;
%     
%     
% end
% 
% for j = 1:trial
%     face_data_temp = squeeze(face_data(:,:,j));
%     face_data(:,:,j) = circshift(face_data_temp',latency_coloumn(j))';
% 
% end
% 
% 
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% randomly chosen 20 trial, repeat 50 times
% rt = 50;
% rn = 10;
% Ncomp = zeros(1,rn);
% ZR = zeros(1,rn);
% Tao_repeat = zeros(1,rn);
% A_repeat = zeros(2,rn);
% Z_repeat = cell(rn,1);
% C_repeat = cell(rn,1);
% 
% matlabpool open
% parfor irn = 1:rn
% rng(irn,'twister');
% Ichose = randperm(trial);
% face_data_part = face_data(:,:,Ichose(1:rt));
% face_invert_data_part = face_invert_data(:,:,Ichose(rt+1:2*rt));
% [C,Z1,Z2,Tao,A1,A2,errflag]=ERP_GA_VB4(face_data_part,face_invert_data_part,10^-8);
% 
% A1d = diag(A1);
% A2d = diag(A2);
% [Y,I] = sort(diag(A1),'descend');
% I = I(Y>0);
% C = C(:,I);
% Z1 = Z1(I,:);
% Z2 = Z2(I,:);
% Tao = Tao(I);
% A1 = diag(A1d(I));
% A2 = diag(A2d(I));
% 
% nc = length(I);
% Ncomp(irn) = nc;
% tempZR = zeros(1,nc);
% tempz = template.JE;
% for ic = 1:nc
%     tc = corrcoef(Z1(ic,:)',tempz');
%     tempZR(ic) = tc(2,1);
% end
% [Y,I_maxzr] = max(abs(tempZR));
% ZR(irn) = Y;
% Tao_repeat(irn) = Tao(I_maxzr);
% A_repeat(:,irn) = [A1(I_maxzr,I_maxzr); A2(I_maxzr,I_maxzr)];
% Z_repeat{irn} = Z1;
% C_repeat{irn} = C;
% 
% end
% matlabpool close
% save('repeat_result_301.mat','Ncomp','Tao_repeat','ZR','A_repeat','Z_repeat','C_repeat',...
%     'Ncomp_pcaica_repeat','ZR_pcaica_repeat','Tao_pcaica_repeat','A_pcaica_repeat','Z_pcaica_repeat','C_pcaica_repeat');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


