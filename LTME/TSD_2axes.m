% ������ ǥ���� ���ؼ� �Ѵ��� ���ΰ����� �ð��� ��
% TSD_2axes.m
clc;clear all;close all
data = load('CTD_TED3.dat')
depth = data(:,2);temp = data(:,3);salinity = data(:,7);DO= data(:,9);
% ���°� ������ 2����� ǥ��
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

% ������ҿ� ������ 2�� ��� ǥ��
subplot(1,2,2)
hl1 = line(DO,-1*depth,'Color','g'); % sound velocity�� green color
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