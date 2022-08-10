% m-file: verdist_exam2.m
%
clc;clear all;close all;
[dep,lon,tem] = textread('temp_TEE.dat','');
dep = -1*dep;                 % 수심 data 변환
x1 = min(lon);x2 = max(lon);
y1 = min(dep);y2 = max(dep); 
Xp = [x1:.0015:x2];
Yp = [y1:.07:y2];
[Xi,Yi] = meshgrid(Xp, Yp);
Zi = griddata(lon,dep,tem, Xi, Yi);
hold on
contourf(Xi,Yi,Zi,30)
shading flat
colorbar('vert');
set(gcf,'Color','w')
title('Temperature vertical distribution','FontSize',13);
ylabel('depth(m)');xlabel('longitude') 
print -dpsc 'Temp_verdist2.ps'