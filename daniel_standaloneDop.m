% ============ Standalone: RMS 및 DOP 분석 ============
clear all;
close all;
clc;
%load('EXP_20250929_group1-2_data.mat');
load('EXP_20251020_group1-2_data.mat');
c = 299792458; % 빛의 속도 (m/s)

%% time setting
ini_epoch = 5*60;              % 초기 미지정수를 위한 정지구간
n_epoch = size(true.GPSTime,1); % 총 데이터 길이
maskangle = 15;                 
nt = 32;
elapsedTime = true.GPSTime - true.GPSTime(1); 

%% Visible satellite selection
vis_sat = zeros(1,32);
for ii = 1:32
    temp_el1 = find(ref.El(:,ii) < maskangle| user.El(:,ii) < maskangle); 
    temp_el2 = find(ref.pr(:,ii)<10 |user.pr(:,ii)<10);
    if isempty(temp_el1) && isempty(temp_el2)
        vis_sat(ii)=1;
    end
end
SV_vis = find(vis_sat==1);
n_vis = length(SV_vis);

%% ECEF -> ENU 변환 행렬 (Rtran) 계산
lat = atan2(ref_xyz(3), sqrt(ref_xyz(1)^2 + ref_xyz(2)^2));
lon = atan2(ref_xyz(2), ref_xyz(1));
Rtran = [-sin(lon),              cos(lon),               0;
         -sin(lat)*cos(lon),     -sin(lat)*sin(lon),     cos(lat);
         cos(lat)*cos(lon),      cos(lat)*sin(lon),      sin(lat)];

%% Standalone 측위 + DOP 계산
 user.pos_standalone = zeros(n_epoch,3);
 user.pos_standalone_enu = zeros(n_epoch,3);
 DOP_values_spp = zeros(n_epoch, 5); % [GDOP, PDOP, HDOP, VDOP, TDOP]
 
 fprintf('Standalone 측위 및 DOP 계산 중...\n');
 
 for ti = 1: n_epoch
     dx = 100;
     if ti == 1
         R_user = ref_xyz; 
     else
         R_user = user.pos_standalone(ti-1, :); 
     end
     B_ur = 0;              
     x_old = [R_user.';B_ur];
     
     iter = 0;
     max_iter = 10;
     
     while(dx > 10^-4 && iter < max_iter)
         H_standalone = zeros(n_vis,4);
         z_standalone =  zeros(n_vis,1);
         for jj = 1: n_vis
             R_su = user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)) - R_user;
             e_hat = R_su/norm(R_su);
             H_standalone(jj,1:3) = e_hat;
             H_standalone(jj,4) = -1;
             z_standalone(jj,1) = e_hat*user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)).'-(user.pr(ti,SV_vis(jj))+user.b(ti,SV_vis(jj))*c);
         end
         
         x = pinv(H_standalone)*z_standalone;
         R_user = x(1:3).';
         B_ur = x(4);
         dx = norm(x-x_old);
         x_old = x;
         iter = iter + 1;
     end
     
     user.pos_standalone(ti,:) = R_user;
     user.pos_standalone_enu(ti,:) = (Rtran*(user.pos_standalone(ti,:)-ref_xyz).').';
     
     % --- DOP 계산 ---
     Q_ecef = inv(H_standalone' * H_standalone);
     R_block = eye(4);
     R_block(1:3, 1:3) = Rtran;
     Q_enu = R_block * Q_ecef * R_block';
     
     sigma_ee_sq = Q_enu(1,1);
     sigma_nn_sq = Q_enu(2,2);
     sigma_uu_sq = Q_enu(3,3);
     sigma_tt_sq = Q_enu(4,4);
     
     HDOP = sqrt(sigma_ee_sq + sigma_nn_sq);
     VDOP = sqrt(sigma_uu_sq);
     PDOP = sqrt(sigma_ee_sq + sigma_nn_sq + sigma_uu_sq);
     TDOP = sqrt(sigma_tt_sq);
     GDOP = sqrt(sigma_ee_sq + sigma_nn_sq + sigma_uu_sq + sigma_tt_sq);
     
     DOP_values_spp(ti, :) = [GDOP, PDOP, HDOP, VDOP, TDOP];
 end
 
 fprintf('계산 완료.\n');

%% --- 1. RMS 오차 수치 분석 (Standalone) ---
analysis_range = (ini_epoch + 1):n_epoch;

errors_spp_enu = user.pos_standalone_enu(analysis_range,:) - true.enu(analysis_range,:);
rms_spp_e = sqrt(mean(errors_spp_enu(:,1).^2));
rms_spp_n = sqrt(mean(errors_spp_enu(:,2).^2));
rms_spp_u = sqrt(mean(errors_spp_enu(:,3).^2));
rms_spp_3d = sqrt(mean(sum(errors_spp_enu.^2, 2)));

fprintf('\n--- [Standalone] RMS 오차 분석 결과 (%.0f초 이후) ---\n', ini_epoch);
fprintf('   East RMS   : %.4f (m)\n', rms_spp_e);
fprintf('   North RMS  : %.4f (m)\n', rms_spp_n);
fprintf('   Up RMS     : %.4f (m)\n', rms_spp_u);
fprintf('   3D RMS     : %.4f (m)\n', rms_spp_3d);
fprintf('--------------------------------------------------\n');

%% --- 2. DOP 그래프 출력 (Standalone) ---
figure('Name', 'DOP (Standalone)');
plot(elapsedTime, DOP_values_spp);
grid on;
title('DOP values over Time (Standalone)');
xlabel('Elapsed Time (s)');
ylabel('DOP (unitless)');
legend('GDOP', 'PDOP', 'HDOP', 'VDOP', 'TDOP', 'Location', 'best');
xlim([0, elapsedTime(end)]);
% Y축 스케일 자동 조정을 위해 ylim 제거 (가독성 향상)