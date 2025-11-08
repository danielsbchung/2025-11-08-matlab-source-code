clear all;
close all;
clc;
%load(['EXP_20250929_group1-2_data.mat'])
load(['EXP_20251020_group1-2_data.mat']) % 1조 데이터 활용하는 경우
c = 299792458; % 빛의 속도 (m/s)

%% time setting
ini_epoch = 5*60;              % 초기 미지정수를 위한 정지구간
n_epoch = size(true.GPSTime,1); % 총 데이터 길이
maskangle = 15;                 
nt = 32;
elapsedTime = true.GPSTime - true.GPSTime(1); % UpTime 그래프용 시간

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

%% ECEF -> ENU 변환 행렬 (Rtran) 계산 (그래프 출력을 위해 추가)
lat = atan2(ref_xyz(3), sqrt(ref_xyz(1)^2 + ref_xyz(2)^2));
lon = atan2(ref_xyz(2), ref_xyz(1));
Rtran = [-sin(lon),              cos(lon),               0;
         -sin(lat)*cos(lon),     -sin(lat)*sin(lon),     cos(lat);
         cos(lat)*cos(lon),      cos(lat)*sin(lon),      sin(lat)];

%% DGPS (변수명은 standalone, 로직은 DGPS)
 user.pos_dgps = zeros(n_epoch,3); % (변수명을 dgps로 명확히 하는 것을 권장)
 user.pos_dgps_enu = zeros(n_epoch,3);
 R_ur = zeros(n_epoch,3);
 deltar = zeros(n_epoch, n_vis);
 dr = zeros(n_epoch, n_vis);
 
 for ti = 1: n_epoch
     
     dx = 100;
     
     % --- 수정된 부분: Warm Start 적용 ---
     if ti == 1
         R_user = ref_xyz; % 첫 번째 epoch는 기준국 위치로 초기화
     else
         R_user = user.pos_dgps(ti-1, :); % 이전 epoch 결과 사용
     end
     % ------------------------------------
     
     % R_user = zeros(3,1).'; % [기존 코드] 초기 위치 가정 (Cold Start)
     B_user = 0;              % 수신기 시계오차
     x_old = [R_user.';B_user];
     
     iter = 0;
     max_iter = 10;
     
     while(dx > 10^-4 && iter < max_iter)
         H_standalone = zeros(n_vis,4);
         z_standalone = zeros(n_vis,1);
         for jj = 1: n_vis
             % Step 2. e_hat 계산
             D_su = user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)) - R_user;
             e_hat = D_su/norm(D_su);
             
             % Step 3. 사용자 위치 계산 하기
             % Step 3.1.1. H matrix 계산
             H_standalone(jj,1:3) = e_hat;
             H_standalone(jj,4) = -1;
             
             % Step 3.1.2. reference 위치 계산시 차이나는 거리 계산
             dr(ti, jj) = norm(ref.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)) - ref_xyz);
             % Step 3.2. measurement vector 계산 (DGPS 로직)
             z_standalone(jj,1) = e_hat*user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)).'-(user.pr(ti,SV_vis(jj))-ref.pr(ti,SV_vis(jj))+dr(ti, jj));
         end
         
         % Step 3.3. user position 계산
         x = pinv(H_standalone)*z_standalone;
         R_user = x(1:3).';
         B_user = x(4);
         
         % Step 4. 사용자 위치 비교 후 update
         dx = norm(x-x_old);
         x_old = x;
         iter = iter + 1;
     end
     
     user.pos_dgps(ti,:) = R_user;
     user.pos_dgps_enu(ti,:) = (Rtran*(user.pos_dgps(ti,:)-ref_xyz).').';
     % error.standalone(ti,:) = user.pos_dgps(ti,:)  - true.xyz(ti,:);
     
 end
 
%% --- 4.2. 참값 vs DGPS (3개 창) ---
% ECEF 그래프 가독성 개선을 위한 상대 위치 계산
mean_pos = mean(true.xyz, 1);
true_relative_xyz = true.xyz - mean_pos;
dgps_relative_xyz = user.pos_dgps - mean_pos; % (user.pos_standalone 대신)

figure('Name', 'DGPS: 3D ECEF');
hold on;
plot3(true_relative_xyz(:,1), true_relative_xyz(:,2), true_relative_xyz(:,3), '-r', 'LineWidth', 1.5);
plot3(dgps_relative_xyz(:,1), dgps_relative_xyz(:,2), dgps_relative_xyz(:,3), '.b');
grid on; axis equal; view(30, 20);
title('True vs DGPS (3D ECEF, Relative)');
legend('True', 'DGPS');
xlabel('Relative X (m)'); ylabel('Relative Y (m)'); zlabel('Relative Z (m)');
hold off;

figure('Name', 'DGPS: North-East');
hold on;
plot(true.enu(:,1), true.enu(:,2), '-r', 'LineWidth', 1.5);
plot(user.pos_dgps_enu(:,1), user.pos_dgps_enu(:,2), '.b');
grid on; axis equal;
title('True vs DGPS (North-East)');
legend('True', 'DGPS');
xlabel('East (m)'); ylabel('North (m)');
hold off;

figure('Name', 'DGPS: Up vs Time');
hold on;
plot(elapsedTime, true.enu(:,3), '-r', 'LineWidth', 1.5);
plot(elapsedTime, user.pos_dgps_enu(:,3), '.b');
grid on;
title('True vs DGPS (Up Trajectory)');
legend('True', 'DGPS');
xlabel('Elapsed Time (s)'); ylabel('Up (m)');
maxTime = ceil(elapsedTime(end)/100)*100;
xticks(0:100:maxTime);
xlim([0, elapsedTime(end)]);
hold off;