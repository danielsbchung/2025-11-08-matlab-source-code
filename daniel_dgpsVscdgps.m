% ============== 실험 2: DGPS vs CDGPS 비교 ==============
clc
clear
close all

% --- 0. 데이터 불러오기 및 기본 설정 ---
load('EXP_20251020_group1-2_data.mat'); % 이동(글자) 데이터
c = 299792458; % 빛의 속도 (m/s)
lambda_L1 = c / 1575.42e6; % L1 파장

% --- 0.1. 기본 설정 ---
n_epoch     = size(true.GPSTime,1);     % 총 데이터 길이
maskangle   = 15;                       % 차폐각 (deg)
nt          = 32;                       % GPS 위성 개수
ini_epoch   = 5*60;                     % 미지정수 계산용 초기 정지시간 (sec)
elapsedTime = true.GPSTime - true.GPSTime(1);

% --- 0.2. 가시 위성 선별 ---
vis_sat = zeros(1,nt);
for ii = 1:nt
    temp_el1 = find(ref.El(:,ii) < maskangle | user.El(:,ii) < maskangle);
    temp_el2 = find(ref.pr(:,ii) < 10 | user.pr(:,ii) < 10);
    if isempty(temp_el1) && isempty(temp_el2)
        vis_sat(ii) = 1;
    end
end
SV_vis  = find(vis_sat == 1);
n_vis   = length(SV_vis);

% --- 0.3. ECEF -> ENU 변환 행렬 계산 ---
lat = atan2(ref_xyz(3), sqrt(ref_xyz(1)^2 + ref_xyz(2)^2));
lon = atan2(ref_xyz(2), ref_xyz(1));
Rtran = [-sin(lon),              cos(lon),               0;
         -sin(lat)*cos(lon),     -sin(lat)*sin(lon),     cos(lat);
         cos(lat)*cos(lon),      cos(lat)*sin(lon),      sin(lat)];
         
%% 1. Standalone 측위 (CDGPS 초기값 계산을 위해 필요)
fprintf('1. Standalone 측위 계산 중 (CDGPS 초기값용)...\n');
user.pos_standalone        = zeros(n_epoch,3);
user.pos_standalone_enu    = zeros(n_epoch,3);
for ti = 1:n_epoch
    dx = 100;
    if ti == 1
        R_user = ref_xyz; 
    else
        R_user = user.pos_standalone(ti-1, :); 
    end
    B_user = 0;
    x_old = [R_user.'; B_user];
    iter = 0;
    max_iter = 10;
    while(dx > 10^-4 && iter < max_iter)
        H_standalone = zeros(n_vis,4);
        z_standalone =  zeros(n_vis,1);
        for jj = 1: n_vis
            sv_pos = user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj));
            R_su = sv_pos - R_user;
            e_hat = R_su/norm(R_su);
            H_standalone(jj,1:3)   = e_hat;
            H_standalone(jj,4)     = -1;
            pr_measured = user.pr(ti,SV_vis(jj)) + user.b(ti,SV_vis(jj))*c;
            z_standalone(jj,1) = e_hat*sv_pos.' - pr_measured;
        end
        x = pinv(H_standalone)*z_standalone;
        R_user = x(1:3).';
        B_user = x(4);
        dx = norm(x - x_old);
        x_old = x;
        iter = iter + 1;
    end
    user.pos_standalone(ti,:)     = R_user;
end

%% 2. DGPS 측위 (Warm Start 적용)
fprintf('2. DGPS 측위 계산 중...\n');
user.pos_dgps     = zeros(n_epoch,3);
user.pos_dgps_enu = zeros(n_epoch,3);
for ti = 1:n_epoch
    dx = 100;
    if ti == 1
        R_user = ref_xyz; 
    else
        R_user = user.pos_dgps(ti-1, :); 
    end
    B_user = 0;
    x_old = [R_user.'; B_user];
    iter = 0;
    max_iter = 10;
    while(dx > 1e-4 && iter < max_iter)
        H_dgps = zeros(n_vis,4);
        z_dgps = zeros(n_vis,1);
        for jj = 1: n_vis
            sv_pos = user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj));
            R_su = sv_pos - R_user;
            e_hat = R_su/norm(R_su);
            H_dgps(jj,1:3) = e_hat;
            H_dgps(jj,4)   = -1;
            dr_ref = norm(ref.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)) - ref_xyz);
            pr_corrected_dgps = user.pr(ti,SV_vis(jj)) - ref.pr(ti,SV_vis(jj)) + dr_ref;
            z_dgps(jj,1) = e_hat*sv_pos.' - pr_corrected_dgps;
        end
        x_new = pinv(H_dgps)*z_dgps;
        R_user = x_new(1:3).';
        B_user = x_new(4);
        dx = norm(x_new - x_old);
        x_old = x_new;
        iter = iter + 1;
    end
    user.pos_dgps(ti,:) = R_user;
    user.pos_dgps_enu(ti,:) = (Rtran*(R_user - ref_xyz).').';
end

%% 3. CDGPS 측위
fprintf('3. CDGPS: 미지정수 계산 중...\n');
% 3.1 기준 위성 선정
avg_el = mean(user.El(1:ini_epoch, SV_vis));
[~, ref_sat_idx_in_vis] = max(avg_el);
ref_sat_prn = SV_vis(ref_sat_idx_in_vis);
fprintf('   기준 위성으로 PRN %d가 선정되었습니다.\n', ref_sat_prn);
other_sats_vis_indices = 1:n_vis;
other_sats_vis_indices(ref_sat_idx_in_vis) = [];
n_other_sats = n_vis - 1;

% 3.2 미지정수(N) 계산 (True 값 사용)
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
fprintf('   미지정수 N 고정 완료.\n');

% 3.3 CDGPS 위치 계산
user.pos_cdgps      = zeros(n_epoch,3);
user.pos_cdgps_enu  = zeros(n_epoch,3);
fprintf('   CDGPS 위치 계산 중...\n');
for ti = 1:n_epoch
    R_user = user.pos_standalone(ti,:); % Standalone 해를 초기값으로 사용
    max_iter = 10;
    tol = 1e-4;
    iter = 0;
    delta_r_norm = Inf;
    while(delta_r_norm > tol && iter < max_iter)
        iter = iter + 1;
        H_cdgps = zeros(n_other_sats, 3);
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
end
fprintf('모든 계산 완료. 그래프 출력 중...\n');

%% 4. 그래프 출력 (DGPS vs CDGPS)
% --- Figure 1: North-East 궤적 ---
figure('Name', '실험 2: DGPS vs CDGPS (North-East)');
hold on;
plot(true.enu(:,1), true.enu(:,2), '-r', 'LineWidth', 2);
plot(user.pos_dgps_enu(:,1), user.pos_dgps_enu(:,2), '.b');
plot(user.pos_cdgps_enu(:,1), user.pos_cdgps_enu(:,2), '.g');
grid on; axis equal;
title('True vs DGPS vs CDGPS (North-East)');
legend('True', 'DGPS', 'CDGPS', 'Location', 'best');
xlabel('East (m)'); ylabel('North (m)');
hold off;

% --- Figure 2: Up vs Time 궤적 ---
figure('Name', '실험 2: DGPS vs CDGPS (Up vs Time)');
hold on;
plot(elapsedTime, true.enu(:,3), '-r', 'LineWidth', 2);
plot(elapsedTime, user.pos_dgps_enu(:,3), '.b');
plot(elapsedTime, user.pos_cdgps_enu(:,3), '.g');
grid on;
title('True vs DGPS vs CDGPS (Up Trajectory)');
legend('True', 'DGPS', 'CDGPS', 'Location', 'best');
xlabel('Elapsed Time (s)'); ylabel('Up (m)');
maxTime = ceil(elapsedTime(end)/100)*100;
xticks(0:100:maxTime);
xlim([0, elapsedTime(end)]);
hold off;