% Temperature horizontal distribution method
% verdist_exam1.m
% input data: temp_sample.dat
% 5개 ctd data를 이용해서 작성된 수온자료
%
clc;clear all;close all;
data = load('temp_sample.dat');
Dist = data(:,1);     % Distance(x) position.
Depth = data(:,2);    % Depth data
Temp = data(:,3);     % Tempeterture data
X1 = round(min(Dist));X2 = round(max(Dist));
Y1 = round(min(Depth));Y2 = round(max(Depth));
Xp = [X1:2:X2];Yp = [Y1:2:Y2];
[Xi,Yi] = meshgrid(Xp, Yp);
Zi = griddata(Dist, Depth, Temp, Xi, Yi);
hold on
% [Xi,Yi] = contour(Zi);
colormap('jet');
pcolor(Xi,Yi,Zi);
hold on;
shading interp;
% clabel(Xi,Yi,'fontsize',10,'color','k');
colorbar('vert');
set(gcf,'Color','w')
title('Temperature vertical distribution','FontSize',13);
ylabel('depth(m)');xlabel('distance') 
print -dpsc 'Temp_verdist1.ps'