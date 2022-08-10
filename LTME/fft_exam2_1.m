% m-file: fft_exam2.m
% data�� interval�� 1hour
%
clc;clear all;close all;
[YY,MM,DD,HR,TH]= textread('jeju0510.dat','');
Dn= datenum(YY,MM,DD,HR,0,0);   % Date serial number
%  �ܼ��� �����ڷ� ǥ��
subplot(2,1,1)
plot(Dn,TH,'b');        % �Ķ��� ������ ǥ��
datetick('x',6);        % 'mm/dd'
xx=TH-mean(TH);
% FFT�� ǥ��
subplot(2,1,2)
fftplot(xx)
xlabel('cph')