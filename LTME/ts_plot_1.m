% plot.m�� �̿��� ���� �� ���� ����������
% m-file: ts_plot.m

clc;clear all;close all
% Data load
data = load('CTD_TED6.dat');
depth = data(:,2);
temp = data(:,3);
salinity = data(:,7);
oxygen = data(:,9);

% ���ɿ� ���� ���º��� �׸���
subplot(1,2,1)
LinHd(1)=plot(temp, -1*depth);  
set(LinHd, 'color', 'b');
grid on
% label �����
xlabel('Temp. (degree)');
ylabel('Depth (m)');
title('TED6 ���� ����������');
subplot(1,2,2)
% ���ɿ� ���� ���к�����
LinHd= plot(salinity, -1*depth);  
set(LinHd,'color','r'); % line���� red��
set(gcf,'color','w'); % ������ white��
grid on   % �׸��� 
% Label �����
xlabel('Salinity (PSU)');
ylabel('Depth (m)');
title('TED6 ���� ����������');