% 임의의 해류자료 표현 
% m-file: curmap_exam1.m
% 
clc;clear all;close all;
m_proj('mercator','lon',[128 137],'lat',[33 39.5]);  % 도법, 범위지정
m_gshhs_i('color','k');            % 준고해상도 지도해안선 검은색
m_gshhs_i('patch',[.8 .8 .8]);       % 지도색을 회색으로 지정
m_grid('box','fancy','tickdir','in','linewidth',1);
m_text(130.5,38.5,'EAST SEA','vertical','top');   %지명을 지도상에 표시한다.
hold on
%curtest.dat파일을 불러서 동방 및 북방성분으로 변경
data = load('curtest.dat')
long = data(:,1);lat = data(:,2);speed = data(:,3);
direction=450-data(:,4);
v=speed.*sin(direction.*pi./180); % 유향유속을 가지고 북방성분 계산
u=speed.*cos(direction.*pi./180); % 유향유속을 가지고 동방성분 계산
m_quiver(long,lat,u,v);           % 해류 Vector 표현
m_plot(long,lat,'ro')
set(gcf,'Color','w');              % 배경색 White 
