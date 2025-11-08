clear all;
close all;
clc;

% --- 0. 데이터 불러오기 ---
% Puzzle 데이터를 로드합니다. (이 파일에는 'true' 또는 'puzzle' 변수가 없습니다)
load(['EXP_20251020_group1-2_data.mat'])
c = 299792458; % 빛의 속도 (m/s)

% --- 'user' 구조체 확인 ---
if ~exist('user', 'var')
    error("'.mat' 파일에 'user' 구조체가 없습니다. .mat 파일 내용을 확인해주세요.");
end

%% time setting
ini_epoch = 5*60;              % 초기 미지정수를 위한 정지구간
n_epoch = size(user.GPSTime,1); % 총 데이터 길이 (user.GPSTime 사용)
maskangle = 15;                 
nt = 32;
elapsedTime = user.GPSTime - user.GPSTime(1); % UpTime 그래프용 시간 (user.GPSTime 사용)

%% Visible satellite selection
PRN = [1:32];
vis_sat = zeros(1,32);
for ii = 1:32
    % mask angle보다 elevation angle이 높고, 측정치가 있는 위성 찾기
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
         
%% 1. Standalone 측위 (CDGPS 계산 없음)
fprintf('1. Standalone 측위 계산 중...\n');
user.pos_standalone        = zeros(n_epoch,3);
% user.pos_standalone_enu    = zeros(n_epoch,3); % ENU 계산 불필요

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
        
        % [수정] 가시 위성이 4개 미만일 경우 pinv 오류 방지
        if size(H_standalone, 1) < 4
             % 이전 위치를 그대로 사용하거나 초기값 사용
             if ti > 1
                 R_user = user.pos_standalone(ti-1, :);
             else
                 R_user = ref_xyz;
             end
             break; % while 루프 탈출
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

fprintf('모든 계산 완료. 지구 중심 ECEF 그래프 출력 중...\n');

%% --- 4. [수정] Standalone 궤적 (지구 중심 Absolute ECEF 기준) ---
% Puzzle 데이터에는 참값이 없으므로, 계산된 궤적만 표시합니다.

% --- 그래프 1: 3D ECEF (Absolute) ---
figure('Name', 'Standalone: 3D (Absolute ECEF)');
hold on;
% [수정] true_data 플롯 제거
plot3(user.pos_standalone(:,1), user.pos_standalone(:,2), user.pos_standalone(:,3), '.k');
grid on; axis equal; view(30, 20);
title('Standalone 궤적 (3D Absolute ECEF)');
legend('Standalone');
xlabel('Absolute X (m)'); ylabel('Absolute Y (m)'); zlabel('Absolute Z (m)');
fprintf('참고: 3D ECEF 그래프는 궤적의 절대 위치(수백만 m)를 표시합니다.\n');
hold off;

% --- 그래프 2: X-Y ECEF (Absolute) ---
figure('Name', 'Standalone: X-Y (Absolute ECEF)');
hold on;
% [수정] true_data 플롯 제거
plot(user.pos_standalone(:,1), user.pos_standalone(:,2), '.k');
grid on; axis equal;
title('Standalone 궤적 (X-Y Absolute ECEF)');
legend('Standalone');
xlabel('Absolute X (m)'); ylabel('Absolute Y (m)');
fprintf('참고: X-Y ECEF 그래프는 궤적의 절대 위치(수백만 m)를 표시합니다.\n');
hold off;

% --- 그래프 3: Z ECEF (Absolute) vs Time ---
figure('Name', 'Standalone: Z (Absolute ECEF) vs Time');
hold on;
% [수정] true_data 플롯 제거
plot(elapsedTime, user.pos_standalone(:,3), '.k');
grid on;
title('Standalone 궤적 (Z Trajectory, Absolute ECEF)');
legend('Standalone');
xlabel('Elapsed Time (s)'); ylabel('Absolute Z (m)');
maxTime = ceil(elapsedTime(end)/100)*100;
xticks(0:100:maxTime);
xlim([0, elapsedTime(end)]);
hold off;