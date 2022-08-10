% Co-Tidal line plot using m_ contour
% in-put file: nao1.dat
% m-file: cotidal_map.m
% 주의: 360도와 0도 사이가 조밀하게 나타나는
% 부분을 개선할 필요 있음

clc;clear all;close all
% korea map display
m_proj('mercator','lon',[123 131],'lat',[30 38])
m_gshhs_h('color','k')
m_gshhs_h('patch',[.84 .84 .75]);  %지도색 지정
%m_grid('linestyle','none','tickdir','in','linewidth',1);
m_grid('box','fancy','tickdir','in');
hold on
% data load
data = load('nao1.dat');
lon = data(:,1);
lat = data(:,2);
q_phase = data(:,4);
q_amp = data(:,3);

[xi,yi] = meshgrid([123:.1:131],[30:.1:38]);
Zi = griddata(lon,lat,q_phase, xi, yi);
[cs,h] = m_contour(xi,yi,Zi,[0:30:360])
clabel(cs,h,'fontsize',6);
hold on;
%m_contour(lon,lat,Zi)
hold off