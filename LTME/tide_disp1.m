% Data load Jeju tide data (2002)
% Data= jeju2002.dat
% m-file = tide_disp1.m
%
clc;clear all;close all;
[YY,MM,DD,HR,MIN,TH]= textread('jeju2002.dat','%n%n%n%n%n%n');
Dn= datenum(YY,MM,DD,HR,MIN,0);   % Date serial number
plot(Dn,TH, 'b.','MarkerSize',2);        % �Ķ��� ������ ǥ�� 
%plot(Dn,TH)
datetick('x',2)     % 'mm/dd/yy'
xlabel('date')
ylabel('Tide Height(cm)')
title('2002�� �������������� ��������')