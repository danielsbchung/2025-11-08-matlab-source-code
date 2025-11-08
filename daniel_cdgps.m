% CDGPS 측위 (두 번째 코드 기반) 및 그래프 3종 출력
clc
clear
close all

% --- 0. 데이터 불러오기 및 기본 설정 ---
%load('EXP_20250929_group1-2_data.mat'); % 정지 데이터
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
    user.pos_standalone_enu(ti,:) = (Rtran*(R_user - ref_xyz).').';
end

%% 2. CDGPS 측위 (두 번째 코드로 구현)
% =========================================================================
%   Part 1: 미지정수(Integer Ambiguity) 결정
% =========================================================================
fprintf('2. CDGPS: 초기 %d초 데이터를 사용하여 미지정수를 계산합니다...\n', ini_epoch);

% 1.1 기준 위성 선정 (초기 데이터 구간에서 평균 앙각이 가장 높은 위성)
avg_el = mean(user.El(1:ini_epoch, SV_vis));
[~, ref_sat_idx_in_vis] = max(avg_el);
ref_sat_prn = SV_vis(ref_sat_idx_in_vis); % 기준 위성의 PRN 번호
fprintf('   기준 위성으로 PRN %d가 선정되었습니다.\n', ref_sat_prn);

% 기준 위성을 제외한 나머지 위성 목록 생성
other_sats_vis_indices = 1:n_vis;
other_sats_vis_indices(ref_sat_idx_in_vis) = [];
n_other_sats = n_vis - 1;

% 1.2 초기 정지 구간(ini_epoch) 동안의 미지정수 추정
N_estimates = zeros(ini_epoch, n_other_sats);

for ti = 1:ini_epoch
    ref_sat_pos = user.svpos(ti, 3*ref_sat_prn-2:3*ref_sat_prn); % 기준 위성 위치
    true_pos_ti = true.xyz(ti,:); % user true ECEF (meters)
    
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
% 1.3 추정된 값들을 평균내고 반올림하여 최종 정수 미지정수(N)로 고정 (cycles)
N_fixed = round(mean(N_estimates, 1));
fprintf('   고정된 정수 미지정수 값 (N_j - N_k) [cycles]:\n');
disp(N_fixed);


% =========================================================================
%   Part 2: CDGPS 위치 해 계산
% =========================================================================
% 결과 저장 변수 초기화
user.pos_cdgps      = zeros(n_epoch,3);
user.pos_cdgps_enu  = zeros(n_epoch,3);
error.cdgps         = zeros(n_epoch,3);
error.cdgps_enu     = zeros(n_epoch,3);

fprintf('   CDGPS 위치 계산 중...\n');
for ti = 1:n_epoch
    
    % 현재 epoch의 초기 위치 추정값 (Standalone 결과 사용)
    R_user = user.pos_standalone(ti,:); % ECEF (row vector)
    
    max_iter = 50; % 무한 루프 방지를 위한 최대 반복 횟수
    tol = 1e-4; % 0.1 mm
    iter = 0;
    delta_r_norm = Inf;

    % 반복 계산을 통해 정확한 위치 해 탐색
    while(delta_r_norm > tol && iter < max_iter)
        iter = iter + 1;
        H_cdgps = zeros(n_other_sats, 3);
        z_cdgps = zeros(n_other_sats, 1);
        
        ref_sat_pos = user.svpos(ti, 3*ref_sat_prn-2:3*ref_sat_prn);
        
        for jj = 1:n_other_sats
            current_sat_idx = other_sats_vis_indices(jj);
            current_sat_prn = SV_vis(current_sat_idx);
            current_sat_pos = user.svpos(ti, 3*current_sat_prn-2:3*current_sat_prn);
            
            % 현재 추정 위치(R_user) 기반 e_hat 벡터 계산
            rho_j0 = current_sat_pos - R_user;
            rho_k0 = ref_sat_pos - R_user;
            rj0 = norm(rho_j0);
            rk0 = norm(rho_k0);
            e_j = rho_j0 / rj0;
            e_k = rho_k0 / rk0;
            
            % H 매트릭스 구성
            H_cdgps(jj, :) = (e_j - e_k);
            
            % 측정값(obs): 이중 차분된 반송파 위상 (meters)
            phi_u_diff = user.cp(ti, current_sat_prn) - user.cp(ti, ref_sat_prn);
            phi_r_diff = ref.cp(ti, current_sat_prn) - ref.cp(ti, ref_sat_prn);
            phi_dd_m = phi_u_diff - phi_r_diff; % (meters)
            
            % 예측값(pred) 계산을 위한 항들
            d_r_j = norm(current_sat_pos - ref_xyz);
            d_r_k = norm(ref_sat_pos - ref_xyz);
            N_fixed_m = N_fixed(jj) * lambda_L1; % (cycles) * (m/cycle) = (m)
            
            % 예측된 기하학적 거리 이중차분 @ R_user
            d_dd_pred = (rj0 - rk0) - (d_r_j - d_r_k);
            
            % z 벡터 구성 (예측값 - 실제 측정값)
            z_cdgps(jj) = ( d_dd_pred ) - (phi_dd_m - N_fixed_m);
        end
        
        % 위치 보정값(delta_r) 계산
        delta_r = H_cdgps \ z_cdgps; % mldivide (robust)
        
        % 사용자 위치 업데이트
        R_user = R_user + delta_r.'; % delta_r은 column vector이므로 transpose
        
        % 수렴 조건 확인
        delta_r_norm = norm(delta_r);
    end
    
    % 최종 계산된 위치를 ECEF -> ENU로 변환
    R_user_enu = (Rtran * (R_user - ref_xyz).').';
    
    % 결과 저장
    user.pos_cdgps(ti,:)      = R_user;
    user.pos_cdgps_enu(ti,:)  = R_user_enu;
    error.cdgps(ti,:)         = R_user - true.xyz(ti,:);
    error.cdgps_enu(ti,:)     = R_user_enu - true.enu(ti,:);
end

fprintf('모든 계산 완료. 그래프 출력 중...\n');

%% --- 4.3. 참값 vs CDGPS (3개 창) ---
% ECEF 그래프 가독성 개선을 위한 상대 위치 계산
mean_pos = mean(true.xyz, 1);
true_relative_xyz = true.xyz - mean_pos;
cdgps_relative_xyz = user.pos_cdgps - mean_pos;

figure('Name', 'CDGPS: 3D ECEF');
hold on;
plot3(true_relative_xyz(:,1), true_relative_xyz(:,2), true_relative_xyz(:,3), '-r', 'LineWidth', 1.5);
plot3(cdgps_relative_xyz(:,1), cdgps_relative_xyz(:,2), cdgps_relative_xyz(:,3), '.g');
grid on; axis equal; view(30, 20);
title('True vs CDGPS (3D ECEF, Relative)');
legend('True', 'CDGPS');
xlabel('Relative X (m)'); ylabel('Relative Y (m)'); zlabel('Relative Z (m)');
hold off;

figure('Name', 'CDGPS: North-East');
hold on;
plot(true.enu(:,1), true.enu(:,2), '-r', 'LineWidth', 1.5);
plot(user.pos_cdgps_enu(:,1), user.pos_cdgps_enu(:,2), '.g');
grid on; axis equal;
title('True vs CDGPS (North-East)');
legend('True', 'CDGPS');
xlabel('East (m)'); ylabel('North (m)');
hold off;

figure('Name', 'CDGPS: Up vs Time');
hold on;
plot(elapsedTime, true.enu(:,3), '-r', 'LineWidth', 1.5);
plot(elapsedTime, user.pos_cdgps_enu(:,3), '.g');
grid on;
title('True vs CDGPS (Up Trajectory)');
legend('True', 'CDGPS');
xlabel('Elapsed Time (s)'); ylabel('Up (m)');
maxTime = ceil(elapsedTime(end)/100)*100;
xticks(0:100:maxTime);
xlim([0, elapsedTime(end)]);
hold off;