% 이중축 표현을 통해서 한눈에 서로값들을 시각적 비교
% TSD_2axes.m
clc;clear all;close all
data = load('CTD_TED3.dat')
depth = data(:,2);temp = data(:,3);salinity = data(:,7);DO= data(:,9);
% 수온과 염분을 2중축상에 표현
subplot(1,2,1)
hl1 = line(temp,-1*depth,'Color','b');
ax1 = gca;
set(ax1,'XColor','b','YColor','k')
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','left',...
           'Color','none',...
           'XColor','k','YColor','k');
hl2 = line(salinity,-1*depth,'Color','r');
grid on
legend([hl1,hl2],'Temp.','Sal.','Location','SouthWest') 
title('TED3 Temp & Salinity')

% 용존산소와 염분을 2중 축상에 표현
subplot(1,2,2)
hl1 = line(DO,-1*depth,'Color','g'); % sound velocity는 green color
ax1 = gca;
set(ax1,'XColor','g','YColor','k')
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','left',...
           'Color','none',...
           'XColor','k','YColor','k');
hl2 = line(salinity,-1*depth,'Color','r');
grid on
legend([hl1,hl2],'DO.','Sal.','Location','SouthWest') 
title('TED3 DO. & Salinity')