function [amari_sp,amari_tw,rse_tw] = algo_eval(A,Z,A_real,Z_real,nc_est,nc_real)
% algo_eval - this function is to evaluate performance of different  
%             amari index of spatial pattern
%             amari index of the temporal wave
%             relative square error of the temporal wave      
%                
% Syntax:  [amari_sp,amari_tw,mse_tw] = algo_eval(A,Z,A_or,Z_or)
%
% Inputs:
%          A: estimated spatial pattern 
%          Z: estimated temporal wave
%          A_real: real spatial pattern
%          Z_real: real temporal wave
%          nc_est: number of components estimated
% Outputs:
%          amari_sp: amari index of spatial pattern
%          amari_tw: amari index of temporal wave 
%          rse_tw: relative square error of the temporal wave, definition: norm(Z-Z_real,'fro')/norm(Z_real,'fro')  
%        
% Example: 
%
%
% Other m-files required: amari.m
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

%%%%%%%%%%% component match %%%%%%%%%%%%%%%%%%%%%
A = squeeze(A);A = A(:,1:nc_est);
Z = squeeze(Z);Z = Z(1:nc_est,:);
A_real = squeeze(A_real);
Z_real = squeeze(Z_real);

pow_A = diag(A'*A);
A_norm = A*diag(pow_A.^-0.5);

pow_A_real = diag(A_real'*A_real);
A_real_norm = A_real*diag(pow_A_real.^-0.5);

pow_Z = diag(Z*Z');
Z_norm = diag(pow_Z.^-0.5)*Z;

pow_Z_real = diag(Z_real*Z_real');
Z_real_norm = diag(pow_Z_real.^-0.5)*Z_real;

[Y,I] = max(abs(Z_norm*Z_real_norm'));

Z_norm = Z_norm(I,:);
A_norm = A_norm(:,I);

amari_sp = amari(A_real_norm,A_norm,nc_real);

amari_tw = amari(Z_real_norm',Z_norm',nc_real);

%%%%%%%%%% mse_tw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_chose = Z(I,:);
rse_tw = norm(Z_chose-Z_real,'fro')/norm(Z_real,'fro');