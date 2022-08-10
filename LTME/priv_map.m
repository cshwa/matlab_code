% Private map making and display method
% m-file: priv_map.m
%
clc;clear all;close all
sta_x = [126.4097,126.4933,126.5222];  % 4개의 point표시
sta_y = [37.2781,37.3183,37.3736];  
m_proj('mercator','lat',[37.2 37.7],'long',[126.2 126.8]);
[x y] = m_ll2xy(sta_x,sta_y);              % 정점 좌표를 그림좌표로 변환
plot(x,y,'ro');hold on;
% 사용자 해안선 자료를 읽고 그림 그리고 부분
m_usercoast('icn_map.mat');

%m_grid              % 위도 및 경도선 그리기
m_grid('linestyle','none','box','fancy','tickdir','in');
print -dpsc icn_map.ps    % 칼라 Post Scrpt file 로 저장