clc
clear
close all

% --- 0. 데이터 불러오기 및 기본 설정 ---
load('EXP_20250929_group1-2_data.mat'); % 정지 데이터
%load('EXP_20251020_group1-2_data.mat'); % 이동(글자) 데이터
c = 299792458; % 빛의 속도 (m/s)
lambda_L1 = c / 1575.42e6; % L1 파장

% --- 0.1. 기본 설정 ---
n_epoch     = size(true.GPSTime,1);     % 총 데이터 길이
maskangle   = 15;                       % 차폐각 (deg)
nt          = 32;                       % GPS 위성 개수
ini_epoch   = 5*60;                     % 미지정수 계산용 초기 정지시간 (sec)
elapsedTime = true.GPSTime - true.GPSTime(1); % UpTime 그래프용 시간

% --- 0.2. 가시 위성 선별 ---
fprintf('가시 위성 선별 중...\n');
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

%% ========================================================================
%   1. Standalone 측위
% =========================================================================
fprintf('1. Standalone 측위 계산 중...\n');
user.pos_standalone = zeros(n_epoch,3);
user.pos_standalone_enu = zeros(n_epoch,3);
error.standalone_enu = zeros(n_epoch,3); % ENU 오차 저장용

for ti = 1: n_epoch
    dx = 100;
    if ti == 1
        R_user = ref_xyz; % 첫 번째 epoch는 기준국 위치로 초기화
    else
        R_user = user.pos_standalone(ti-1, :); % 이전 epoch 결과 사용
    end
    B_ur = 0;              % 수신기 시계오차
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
    error.standalone_enu(ti,:) = user.pos_standalone_enu(ti,:)  - true.enu(ti,:);
end

%% ========================================================================
%   2. DGPS 측위
% =========================================================================
fprintf('2. DGPS 측위 계산 중...\n');
user.pos_dgps = zeros(n_epoch,3); 
user.pos_dgps_enu = zeros(n_epoch,3);
error.dgps_enu = zeros(n_epoch,3); % ENU 오차 저장용
dr = zeros(n_epoch, n_vis);

for ti = 1: n_epoch
    dx = 100;
    if ti == 1
        R_user = ref_xyz; % 첫 번째 epoch는 기준국 위치로 초기화
    else
        R_user = user.pos_dgps(ti-1, :); % 이전 epoch 결과 사용
    end
    B_user = 0;              % 수신기 시계오차
    x_old = [R_user.';B_user];
    
    iter = 0;
    max_iter = 10;
    
    while(dx > 10^-4 && iter < max_iter)
        H_standalone = zeros(n_vis,4); % (H 매트릭스 변수명 재활용)
        z_standalone = zeros(n_vis,1); % (z 벡터 변수명 재활용)
        for jj = 1: n_vis
            D_su = user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)) - R_user;
            e_hat = D_su/norm(D_su);
            
            H_standalone(jj,1:3) = e_hat;
            H_standalone(jj,4) = -1;
            
            dr(ti, jj) = norm(ref.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)) - ref_xyz);
            z_standalone(jj,1) = e_hat*user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)).'-(user.pr(ti,SV_vis(jj))-ref.pr(ti,SV_vis(jj))+dr(ti, jj));
        end
        
        x = pinv(H_standalone)*z_standalone;
        R_user = x(1:3).';
        B_user = x(4);
        
        dx = norm(x-x_old);
        x_old = x;
        iter = iter + 1;
    end
    
    user.pos_dgps(ti,:) = R_user;
    user.pos_dgps_enu(ti,:) = (Rtran*(user.pos_dgps(ti,:)-ref_xyz).').';
    error.dgps_enu(ti,:) = user.pos_dgps_enu(ti,:)  - true.enu(ti,:);
end

%% ========================================================================
%   3. CDGPS 측위
% =========================================================================
%   Part 1: 미지정수(Integer Ambiguity) 결정
% =========================================================================
fprintf('3. CDGPS: 초기 %d초 데이터로 미지정수 계산 중...\n', ini_epoch);
avg_el = mean(user.El(1:ini_epoch, SV_vis));
[~, ref_sat_idx_in_vis] = max(avg_el);
ref_sat_prn = SV_vis(ref_sat_idx_in_vis); % 기준 위성의 PRN 번호
fprintf('   CDGPS 기준 위성: PRN %d\n', ref_sat_prn);

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
        
        % 측정된 반송파 위상의 이중 차분 (meters)
        phi_u_diff = user.cp(ti, current_sat_prn) - user.cp(ti, ref_sat_prn);
        phi_r_diff = ref.cp(ti, current_sat_prn) - ref.cp(ti, ref_sat_prn);
        phi_dd_m = phi_u_diff - phi_r_diff; % (meters)
        
        % 참값을 이용한 기하학적 거리의 이중 차분 (meters)
        d_u_j = norm(current_sat_pos - true_pos_ti);
        d_u_k = norm(ref_sat_pos - true_pos_ti);
        d_r_j = norm(current_sat_pos - ref_xyz);
        d_r_k = norm(ref_sat_pos - ref_xyz);
        d_dd = (d_u_j - d_u_k) - (d_r_j - d_r_k);
        
        % 미지정수 추정치 계산 (cycles)
        N_estimates(ti, jj) = (phi_dd_m - d_dd) / lambda_L1;
    end
end
N_fixed = round(mean(N_estimates, 1));
fprintf('   고정된 정수 미지정수 값 (N_j - N_k) [cycles]:\n');
disp(N_fixed);

% =========================================================================
%   Part 2: CDGPS 위치 해 계산
% =========================================================================
user.pos_cdgps      = zeros(n_epoch,3);
user.pos_cdgps_enu  = zeros(n_epoch,3);
error.cdgps_enu     = zeros(n_epoch,3); % ENU 오차 저장용

fprintf('   CDGPS 위치 계산 중...\n');
for ti = 1:n_epoch
    
    % 현재 epoch의 초기 위치 추정값 (Standalone 결과 사용)
    R_user = user.pos_standalone(ti,:); % ECEF (row vector)
    
    max_iter = 50; 
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
    
    R_user_enu = (Rtran * (R_user - ref_xyz).').';
    
    user.pos_cdgps(ti,:)      = R_user;
    user.pos_cdgps_enu(ti,:)  = R_user_enu;
    error.cdgps_enu(ti,:)     = R_user_enu - true.enu(ti,:);
end

%% ========================================================================
%   5. 통합 오차 그래프 (참값 기준)
% =========================================================================
fprintf('모든 계산 완료. 오차 그래프 출력 중...\n');

% --- 그래프 1: East-North 오차 (참값 = 0,0) ---
figure('Name', '측위 방식별 East-North 오차 (참값 기준)');
hold on;
plot(error.standalone_enu(:,1), error.standalone_enu(:,2), '.k', 'MarkerSize', 4); % Standalone (검은색 점)
plot(error.dgps_enu(:,1), error.dgps_enu(:,2), '.b', 'MarkerSize', 4); % DGPS (파란색 점)
plot(error.cdgps_enu(:,1), error.cdgps_enu(:,2), '.g', 'MarkerSize', 8); % CDGPS (초록색 점)
plot(0, 0, '+r', 'MarkerSize', 12, 'LineWidth', 2); % 참값 위치 (빨간색 십자가)
hold off;
grid on;
axis equal;
title('East-North 오차 (참값 = 원점)');
xlabel('East Error (m)');
ylabel('North Error (m)');
legend('Standalone', 'DGPS', 'CDGPS', 'True Position (0,0)');
set(gca, 'FontSize', 12);

% --- 그래프 2: Up 오차 (시간별) ---
figure('Name', '측위 방식별 Up 오차 (시간별)');
hold on;
plot(elapsedTime, error.standalone_enu(:,3), '-k', 'LineWidth', 1); % Standalone (검은색 실선)
plot(elapsedTime, error.dgps_enu(:,3), '-b', 'LineWidth', 1); % DGPS (파란색 실선)
plot(elapsedTime, error.cdgps_enu(:,3), '-g', 'LineWidth', 1.5); % CDGPS (초록색 굵은 실선)
yline(0, 'r--', 'LineWidth', 1.5, 'Label', ' '); % 참값 (빨간색 점선)
hold off;
grid on;
title('Up 오차 (시간별)');
xlabel('Elapsed Time (s)');
ylabel('Up Error (m)');
legend('Standalone', 'DGPS', 'CDGPS');
xlim([0, elapsedTime(end)]);
maxTime = ceil(elapsedTime(end)/100)*100;
xticks(0:100:maxTime);
set(gca, 'FontSize', 12);

fprintf('그래프 출력을 완료했습니다.\n');