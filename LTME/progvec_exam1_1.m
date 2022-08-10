% 조류 연속벡터도를 작성하는 예 
% m-file: progvec_exam1.m
% by cho jae gab & peter

clc;clear all;close all;
% Data 불러오기
[yy,mm,dd,hr,min,temp,speed,n_dir,u_vel,v_vel]=textread('suyeong9805.dat','');
ti = 3600;          % time interval 1hour = 3600sec
n = length(yy);     % data의 size

% 시작점 위치 (0,0)
nu = zeros(n,1); nv =zeros(n,1);
for i=2:1:n
    nu(i) = nu(i-1) + u_vel(i-1)/100000*ti;   % nu는 u_vel의 다음 점 정보
    nv(i) = nv(i-1) + v_vel(i-1)/100000*ti;   % nv는 v_vel의 다음 점 정보
end
% 시작점 기준위치 (0,0)
plot(0,0,'ro');hold on
plot(nu,nv,'r');hold on;     % pgrogressive vector 도시
plot([-50 50],[0 0],':b');   % plot([x1 x2],[y1 y2],'b') 
plot([0 0],[-50 50],':b');
xlim([-50,200]);ylim([-50,200]);
axis equal
legend('start point','Depth 5m','Location','SouthEast');
title('Progressive Vector Diagram');   % 제목
xlabel('x-direction(km)');ylabel('y-direction(km)');         % y축 이름