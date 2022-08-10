% 인천 경기만 부근에서 관측된 조류관측 조류 data의 Scatterplot 예
% m-file: scattermap_exam.m
% input file: 00IC-1.dat, 02IC-2.dat 

clc;clear all;close all
% data load 00ic-1.dat
[yy,mm,dd,hr,min,spd,n_dir,temp]=textread('00IC-1.dat','','headerlines',15);
dir = 450 - n_dir;               % 각도변환
u1 = spd.*cos(dir.*pi./180);     % speed와 정북기준각에서 북방성분 계산
v1 = spd.*sin(dir.*pi./180);     % speed와 정북기준각에서 동방성분 계산
u1 = 126.4097+u1/3000            % 원점을 위도 126.4097로 변환된 u1값
v1 = 37.2781+v1/3000             % 원점을 경도 37.2781로 변환  v1값
% data load 02ic-2.dat
[yy,mm,dd,hr,min,spd,n_dir,temp,cond]=textread('02IC-2.dat','','headerlines',15);
dir = 450 - n_dir;                   % 각도변환
u2 = spd.*cos(dir.*pi./180);     % speed와 정북기준각에서 북방성분 계산
v2 = spd.*sin(dir.*pi./180);      % speed와 정북기준각에서 동방성분 계산
u2 = 126.5222+u2/3000           % 원점을 위도 126.4097로 변환된 u1값
v2 = 37.3736+v2/3000             % 원점을 경도 37.2781로 변환  v1값
% 사용자 해안선 자료를 읽고 그림 그리고 부분
m_proj('mercator','lat',[37.2 37.7],'long',[126.2 126.8]);
m_usercoast('icn_map.mat');       
m_grid('linestyle','none','box','fancy','tickdir','in');  % 위도 및 경도선 그리기
hold on
% scatterplot 표현
m_plot(u1,v1,'r.','MarkerSize',3)   % scatter plot 
hold on
m_plot(u2,v2,'r.','MarkerSize',3)   % scatter plot 
hold on
% 기본 정점원점 그리기
m_plot(126.4097,37.2781,'b+','MarkerSize',5)
m_plot(126.5222,37.3736,'b+','MarkerSize',5)