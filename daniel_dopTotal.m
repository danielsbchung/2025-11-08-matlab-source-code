clc
clear
close all

% --- 0. 데이터 불러오기 및 기본 설정 ---
%load('EXP_20250929_group1-2_data.mat'); % 정지 데이터
load('EXP_20251020_group1-2_data.mat'); % 이동(글자) 데이터
c = 299792458; % 빛의 속도 (m/s)
lambda_L1 = c / 1575.42e6; % L1 파장

% --- 0.1. 기본 설정 ---
ini_epoch   = 5*60;                     % 초기 미지정수 계산용 초기 정지시간 (sec)
n_epoch     = size(true.GPSTime,1);     % 총 데이터 길이
maskangle   = 15;                       % 차폐각 (deg)
nt          = 32;                       % GPS 위성 개수
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
%   1. Standalone 측위 + DOP 계산
% =========================================================================
fprintf('1. Standalone 측위 및 DOP 계산 중...\n');
user.pos_standalone = zeros(n_epoch,3);
DOP_values_spp = zeros(n_epoch, 5); % [GDOP, PDOP, HDOP, VDOP, TDOP]

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
    H_standalone = zeros(n_vis,4); % H 행렬 선언
    
    while(dx > 10^-4 && iter < max_iter)
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
    
    % --- DOP 계산 (Standalone) ---
    Q_ecef = inv(H_standalone' * H_standalone);
    R_block = eye(4);
    R_block(1:3, 1:3) = Rtran;
    Q_enu = R_block * Q_ecef * R_block';
    
    sigma_ee_sq = Q_enu(1,1);
    sigma_nn_sq = Q_enu(2,2);
    sigma_uu_sq = Q_enu(3,3);
    sigma_tt_sq = Q_enu(4,4);
    
    GDOP = sqrt(sigma_ee_sq + sigma_nn_sq + sigma_uu_sq + sigma_tt_sq);
    PDOP = sqrt(sigma_ee_sq + sigma_nn_sq + sigma_uu_sq);
    HDOP = sqrt(sigma_ee_sq + sigma_nn_sq);
    VDOP = sqrt(sigma_uu_sq);
    TDOP = sqrt(sigma_tt_sq);
    
    DOP_values_spp(ti, :) = [GDOP, PDOP, HDOP, VDOP, TDOP];
end

%% ========================================================================
%   2. DGPS 측위 + DOP 계산
% =========================================================================
fprintf('2. DGPS 측위 및 DOP 계산 중...\n');
user.pos_dgps = zeros(n_epoch,3);
DOP_values_dgps = zeros(n_epoch, 5); % [GDOP, PDOP, HDOP, VDOP, TDOP]
dr = zeros(n_epoch, n_vis);

for ti = 1: n_epoch
    dx = 100;
    if ti == 1
        R_user = ref_xyz; 
    else
        R_user = user.pos_dgps(ti-1, :); 
    end
    B_user = 0;             
    x_old = [R_user.';B_user];
    
    iter = 0;
    max_iter = 10;
    H_standalone = zeros(n_vis,4); % H 행렬 이름은 SPP와 동일 (4x4)
    
    while(dx > 10^-4 && iter < max_iter)
        z_standalone = zeros(n_vis,1);
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
     
    % --- DOP 계산 (DGPS) ---
    % H 행렬은 Standalone과 동일하므로 DOP 값도 동일함
    Q_ecef = inv(H_standalone' * H_standalone);
    R_block = eye(4);
    R_block(1:3, 1:3) = Rtran;
    Q_enu = R_block * Q_ecef * R_block';
     
    sigma_ee_sq = Q_enu(1,1);
    sigma_nn_sq = Q_enu(2,2);
    sigma_uu_sq = Q_enu(3,3);
    sigma_tt_sq = Q_enu(4,4);
     
    GDOP = sqrt(sigma_ee_sq + sigma_nn_sq + sigma_uu_sq + sigma_tt_sq);
    PDOP = sqrt(sigma_ee_sq + sigma_nn_sq + sigma_uu_sq);
    HDOP = sqrt(sigma_ee_sq + sigma_nn_sq);
    VDOP = sqrt(sigma_uu_sq);
    TDOP = sqrt(sigma_tt_sq);
     
    DOP_values_dgps(ti, :) = [GDOP, PDOP, HDOP, VDOP, TDOP];
end
 
%% ========================================================================
%   3. CDGPS 측위 + DOP 계산
% =========================================================================
fprintf('3. CDGPS DOP 계산 중 (PDOP/HDOP/VDOP만 계산됨)...\n');
% 3.1 미지정수 결정 (기준 위성 선정을 위해 필요)
avg_el = mean(user.El(1:ini_epoch, SV_vis));
[~, ref_sat_idx_in_vis] = max(avg_el);
ref_sat_prn = SV_vis(ref_sat_idx_in_vis);
other_sats_vis_indices = 1:n_vis;
other_sats_vis_indices(ref_sat_idx_in_vis) = [];
n_other_sats = n_vis - 1;

% 3.2 CDGPS DOP 계산
DOP_cdgps = zeros(n_epoch, 3); % [PDOP, HDOP, VDOP]

for ti = 1:n_epoch
    % R_user 근사값으로 ref_xyz 사용 (DOP 경향성 파악에는 무리 없음)
    R_user = ref_xyz; 
    
    H_cdgps = zeros(n_other_sats, 3); % H 행렬 선언
    ref_sat_pos = user.svpos(ti, 3*ref_sat_prn-2:3*ref_sat_prn);
    
    for jj = 1:n_other_sats
        current_sat_idx = other_sats_vis_indices(jj);
        current_sat_prn = SV_vis(current_sat_idx);
        current_sat_pos = user.svpos(ti, 3*current_sat_prn-2:3*current_sat_prn);
        
        rho_j0 = current_sat_pos - R_user;
        rho_k0 = ref_sat_pos - R_user;
        e_j = rho_j0 / norm(rho_j0);
        e_k = rho_k0 / norm(rho_k0);
        H_cdgps(jj, :) = (e_j - e_k);
    end
    
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
fprintf('모든 계산 완료. 그래프 출력 중...\n');

%% ========================================================================
%   4. 통합 DOP 그래프 출력 (요청 사항 반영)
% =========================================================================

% --- 그래프 1: Standalone vs DGPS (GDOP) ---
figure('Name', 'GDOP 비교 (Standalone vs DGPS)');
hold on;
plot(elapsedTime, DOP_values_spp(:,1), '-r', 'LineWidth', 4);
plot(elapsedTime, DOP_values_dgps(:,1), '--b', 'LineWidth', 4);
hold off;
grid on;
title('GDOP (Standalone vs DGPS)');
xlabel('Elapsed Time (s)');
ylabel('GDOP (unitless)');
% (참고: 두 값은 H행렬이 동일하므로 완전히 겹쳐서 검은색 선만 보일 것입니다)
legend('Standalone', 'DGPS');
xlim([0, elapsedTime(end)]);

% --- 그래프 2: CDGPS (PDOP, HDOP, VDOP) ---
figure('Name', 'DOP (CDGPS)');
plot(elapsedTime, DOP_cdgps);
grid on;
title('DOP values over Time (CDGPS)');
xlabel('Elapsed Time (s)');
ylabel('DOP (unitless)');
legend('PDOP', 'HDOP', 'VDOP', 'Location', 'best');
xlim([0, elapsedTime(end)]);

fprintf('그래프 출력을 완료했습니다.\n');