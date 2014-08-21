clear all
close all

addpath /home/wch/Documents/ERP_freq_max/extr_onesrc2
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

% select the useful channel as setting
[chanindex,mycnt,eloc] = chanselect(mycnt, eloc, chandelete, chanuseful);
channum = channum-length(chandelete);

ana_begin = 0.01*srate; % begining
ana_end = 0.266*srate;    % end
[target_data,face_data,face_invert_data,object_data] = dividedata(mycnt,ana_begin,ana_end);
[chan,len,trial] = size(face_data);

for i = 1:trial
    face_data(:,:,i) = detrend(squeeze(face_data(:,:,i))','constant')';
    object_data(:,:,i) = detrend(squeeze(object_data(:,:,i))','constant')';
end

figure_flag = true;
%%%% make a template 
load tempz10.mat
ERP_template = tempz;

Pv = zeros(7,1);
peak_diff = zeros(7,1);
TPR = zeros(7,2*trial);
FPR = zeros(7,2*trial);
auc = zeros(7,1);
nc_est = ones(7,1);
%%%%%%%%%%%%%% SSF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

[Pv(1),peak_diff(1),Y_face,TPR(1,:),FPR(1,:)]= ERP_Mdistance(object_data,face_data,SF_face,[round(0.140*srate):round(0.190*srate)],figure_flag);

% auc(1) = myauc(fpr_face,tpr_face);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLSOBI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = [0,0.5,0.8,1];
for i = 1:4
    [A_plsobi_face,Z_plsobi_face,SF_plsobi_face,nc_est(i+1)] = real_data_PLSOBI(face_data,ERP_template,lambda(i));

    %%%%% just for N170 first
    [Pv(i+1),peak_diff(i+1),Y_plsobi,TPR(i+1,:),FPR(i+1,:)]= ERP_Mdistance(object_data,face_data,SF_plsobi_face,[round(0.140*srate):round(0.190*srate)],figure_flag);
end
%%%%%%%%%%%%%%%%%%%%% SCA1: use synthetic N170 as template %%%%%%%%%%%%%%%% 
[A_sca_face,Z_sca_face,SF_sca_face,nc_est_sca_face] = real_data_SCA(face_data,ERP_template);
cofmat = corrcoef(Z_sca_face, ERP_template);
if cofmat(2,1) < 0
    SF_sca_face = -SF_sca_face;
end
[Pv(6),peak_diff(6),Y_sca1,TPR(6,:),FPR(6,:)]= ERP_Mdistance(object_data,face_data,SF_sca_face,[round(0.140*srate):round(0.190*srate)],figure_flag);

%%%%%%%%%%%%%%%%%%%% SCA2: use square wave as tempalte %%%%%%%%%%%%%%%%%%%%
ERP_template = zeros(1,256);
ERP_template(150:190) = 1;
[A_sca_face,Z_sca_face,SF_sca_face,nc_est_sca_face] = real_data_SCA(face_data,ERP_template);
cofmat = corrcoef(Z_sca_face, ERP_template);
if cofmat(2,1) < 0
    SF_sca_face = -SF_sca_face;
end
[Pv(7),peak_diff(7),Y_sca2,TPR(7,:),FPR(7,:)]= ERP_Mdistance(object_data,face_data,SF_sca_face,[round(0.140*srate):round(0.190*srate)],figure_flag);


%%%%%%%%%%%%%%%%%%% plot figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(FPR(1,:),TPR(1,:),'k');
hold on;
plot(FPR(4,:),TPR(4,:),'r');
hold on;
plot(FPR(6,:),TPR(6,:),'b');
hold on;
plot(FPR(7,:),TPR(7,:),'g');

for i = 1:7
    auc(i) = myauc(FPR(i,:),TPR(i,:));
end
%%%%% use different type of template
% synthetic data
% square wave

%%%%% data length is 512ms







