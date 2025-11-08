% ============== 실험 1: Standalone vs DGPS 비교 ==============
clc
clear
close all

% --- 0. 데이터 불러오기 및 기본 설정 ---
load('EXP_20250929_group1-2_data.mat'); % 이동(글자) 데이터
c = 299792458; % 빛의 속도 (m/s)

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
         
%% 1. Standalone 측위 (Warm Start 적용)
fprintf('1. Standalone 측위 계산 중...\n');
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

fprintf('계산 완료. 그래프 출력 중...\n');

%% 3. 그래프 출력 (Standalone vs DGPS)
% --- Figure 1: North-East 궤적 ---
figure('Name', '실험 1: Standalone vs DGPS (North-East)');
hold on;
plot(true.enu(:,1), true.enu(:,2), '-r', 'LineWidth', 2);
plot(user.pos_standalone_enu(:,1), user.pos_standalone_enu(:,2), '.k');
plot(user.pos_dgps_enu(:,1), user.pos_dgps_enu(:,2), '.b');
grid on; axis equal;
title('True vs Standalone vs DGPS (North-East)');
legend('True', 'Standalone', 'DGPS', 'Location', 'best');
xlabel('East (m)'); ylabel('North (m)');
hold off;

% --- Figure 2: Up vs Time 궤적 ---
figure('Name', '실험 1: Standalone vs DGPS (Up vs Time)');
hold on;
plot(elapsedTime, true.enu(:,3), '-r', 'LineWidth', 2);
plot(elapsedTime, user.pos_standalone_enu(:,3), '.k');
plot(elapsedTime, user.pos_dgps_enu(:,3), '.b');
grid on;
title('True vs Standalone vs DGPS (Up Trajectory)');
legend('True', 'Standalone', 'DGPS', 'Location', 'best');
xlabel('Elapsed Time (s)'); ylabel('Up (m)');
maxTime = ceil(elapsedTime(end)/100)*100;
xticks(0:100:maxTime);
xlim([0, elapsedTime(end)]);
hold off;