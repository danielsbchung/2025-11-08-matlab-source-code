% ============ CDGPS: RMS 및 DOP 분석 ============
clear all;
close all;
clc;
%load('EXP_20250929_group1-2_data.mat');
load('EXP_20251020_group1-2_data.mat');
c = 299792458; % 빛의 속도 (m/s)
lambda_L1 = c / 1575.42e6; % L1 파장

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
    if isempty(temp_el1)
        if isempty(temp_el2)
            vis_sat(ii)=1;
        end
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

%% 1. Standalone 측위 (CDGPS 초기값 계산을 위해 필요)
 user.pos_standalone = zeros(n_epoch,3);
 
 fprintf('1. Standalone 측위 계산 중 (CDGPS 초기값용)...\n');
 
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
 end
 
%% 2. CDGPS 측위 + DOP 계산
% 2.1 미지정수 결정
fprintf('2. CDGPS 미지정수 계산 중...\n');
avg_el = mean(user.El(1:ini_epoch, SV_vis));
[~, ref_sat_idx_in_vis] = max(avg_el);
ref_sat_prn = SV_vis(ref_sat_idx_in_vis);
other_sats_vis_indices = 1:n_vis;
other_sats_vis_indices(ref_sat_idx_in_vis) = [];
n_other_sats = n_vis - 1;

N_estimates = zeros(ini_epoch, n_other_sats);
for ti = 1:ini_epoch
    ref_sat_pos = user.svpos(ti, 3*ref_sat_prn-2:3*ref_sat_prn);
    true_pos_ti = true.xyz(ti,:);
    for jj = 1:n_other_sats
        current_sat_idx = other_sats_vis_indices(jj);
        current_sat_prn = SV_vis(current_sat_idx);
        current_sat_pos = user.svpos(ti, 3*current_sat_prn-2:3*current_sat_prn);
        
        phi_u_diff = user.cp(ti, current_sat_prn) - user.cp(ti, ref_sat_prn);
        phi_r_diff = ref.cp(ti, current_sat_prn) - ref.cp(ti, ref_sat_prn);
        phi_dd_m = phi_u_diff - phi_r_diff;
        
        d_u_j = norm(current_sat_pos - true_pos_ti);
        d_u_k = norm(ref_sat_pos - true_pos_ti);
        d_r_j = norm(current_sat_pos - ref_xyz);
        d_r_k = norm(ref_sat_pos - ref_xyz);
        d_dd = (d_u_j - d_u_k) - (d_r_j - d_r_k);
        
        N_estimates(ti, jj) = (phi_dd_m - d_dd) / lambda_L1;
    end
end
N_fixed = round(mean(N_estimates, 1));
fprintf('   미지정수 N 고정 완료. (기준 위성: PRN %d)\n', ref_sat_prn);

% 2.2 CDGPS 위치 해 및 DOP 계산
user.pos_cdgps      = zeros(n_epoch,3);
user.pos_cdgps_enu  = zeros(n_epoch,3);
DOP_cdgps = zeros(n_epoch, 3); % [PDOP, HDOP, VDOP]

fprintf('   CDGPS 위치 및 DOP 계산 중...\n');
for ti = 1:n_epoch
    R_user = user.pos_standalone(ti,:); % Standalone 결과를 초기값으로 사용
    max_iter = 10;
    tol = 1e-4;
    iter = 0;
    delta_r_norm = Inf;
    
    H_cdgps = zeros(n_other_sats, 3); % H 행렬 선언
    
    while(delta_r_norm > tol && iter < max_iter)
        iter = iter + 1;
        z_cdgps = zeros(n_other_sats, 1);
        ref_sat_pos = user.svpos(ti, 3*ref_sat_prn-2:3*ref_sat_prn);
        
        for jj = 1:n_other_sats
            current_sat_idx = other_sats_vis_indices(jj);
            current_sat_prn = SV_vis(current_sat_idx);
            current_sat_pos = user.svpos(ti, 3*current_sat_prn-2:3*current_sat_prn);
            
            rho_j0 = current_sat_pos - R_user;
            rho_k0 = ref_sat_pos - R_user;
            rj0 = norm(rho_j0);
            rk0 = norm(rho_k0);
            e_j = rho_j0 / rj0;
            e_k = rho_k0 / rk0;
            H_cdgps(jj, :) = (e_j - e_k);
            
            phi_u_diff = user.cp(ti, current_sat_prn) - user.cp(ti, ref_sat_prn);
            phi_r_diff = ref.cp(ti, current_sat_prn) - ref.cp(ti, ref_sat_prn);
            phi_dd_m = phi_u_diff - phi_r_diff;
            
            d_r_j = norm(current_sat_pos - ref_xyz);
            d_r_k = norm(ref_sat_pos - ref_xyz);
            N_fixed_m = N_fixed(jj) * lambda_L1;
            d_dd_pred = (rj0 - rk0) - (d_r_j - d_r_k);
            
            z_cdgps(jj) = ( d_dd_pred ) - (phi_dd_m - N_fixed_m);
        end
        
        delta_r = H_cdgps \ z_cdgps;
        R_user = R_user + delta_r.';
        delta_r_norm = norm(delta_r);
    end
    
    user.pos_cdgps(ti,:)      = R_user;
    user.pos_cdgps_enu(ti,:)  = (Rtran * (R_user - ref_xyz).').';
    
    % --- DOP 계산 (CDGPS) ---
    Q_ecef_cdgps = inv(H_cdgps' * H_cdgps); % (3 x 3)
    Q_enu_cdgps = Rtran * Q_ecef_cdgps * Rtran';
    
    sigma_ee_sq = Q_enu_cdgps(1,1);
    sigma_nn_sq = Q_enu_cdgps(2,2);
    sigma_uu_sq = Q_enu_cdgps(3,3);
    
    HDOP = sqrt(sigma_ee_sq + sigma_nn_sq);
    VDOP = sqrt(sigma_uu_sq);
    PDOP = sqrt(sigma_ee_sq + sigma_nn_sq + sigma_uu_sq);
    
    DOP_cdgps(ti, :) = [PDOP, HDOP, VDOP];
end
fprintf('모든 계산 완료.\n');

%% --- 1. RMS 오차 수치 분석 (CDGPS) ---
analysis_range = (ini_epoch + 1):n_epoch;

errors_cdgps_enu = user.pos_cdgps_enu(analysis_range,:) - true.enu(analysis_range,:);
rms_cdgps_e = sqrt(mean(errors_cdgps_enu(:,1).^2));
rms_cdgps_n = sqrt(mean(errors_cdgps_enu(:,2).^2));
rms_cdgps_u = sqrt(mean(errors_cdgps_enu(:,3).^2));
rms_cdgps_3d = sqrt(mean(sum(errors_cdgps_enu.^2, 2)));

fprintf('\n--- [CDGPS] RMS 오차 분석 결과 (%.0f초 이후) ---\n', ini_epoch);
fprintf('   East RMS   : %.4f (m)\n', rms_cdgps_e);
fprintf('   North RMS  : %.4f (m)\n', rms_cdgps_n);
fprintf('   Up RMS     : %.4f (m)\n', rms_cdgps_u);
fprintf('   3D RMS     : %.4f (m)\n', rms_cdgps_3d);
fprintf('--------------------------------------------------\n');

%% --- 2. DOP 그래프 출력 (CDGPS) ---
figure('Name', 'DOP (CDGPS)');
plot(elapsedTime, DOP_cdgps);
grid on;
title('DOP values over Time (CDGPS)');
xlabel('Elapsed Time (s)');
ylabel('DOP (unitless)');
legend('PDOP', 'HDOP', 'VDOP', 'Location', 'best');
xlim([0, elapsedTime(end)]);
% Y축 스케일 자동 조정을 위해 ylim 제거 (가독성 향상)