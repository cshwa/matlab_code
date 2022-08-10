% M_Map을 이용한 등수심선 도시
% m?file: contourmap_exam1.m
% input put data: Depth_w229.dat (해도 W229 수심data)
% 

clc;clear all;close all
% South korea map display
m_proj('mercator','lon',[127 129.5],'lat',[34 35.5])
m_gshhs_h('color','k')
m_gshhs_h('patch',[.84 .84 .75]);  %지도 색 지정
m_grid('box','fancy','tickdir','in');
hold on
% contour line making
data = load('Depth_w229.dat');
lon = data(:,1);
lat = data(:,2);
depth = data(:,3);
x1 = min(lon);x2 = max(lon);
y1 = min(lat);y2 = max(lat);
Xp = [x1:.01:x2];
Yp = [y1:.01:y2];
[Xi,Yi] = meshgrid(Xp,Yp);
Zi = griddata(lon,lat,depth, Xi, Yi);
[cs,h] = m_contour(Xi,Yi,Zi,[200:10:0]);
clabel(cs,h,'fontsize',6);
set(gcf,'color','w');
hold off