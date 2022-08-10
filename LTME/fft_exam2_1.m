% m-file: fft_exam2.m
% data의 interval은 1hour
%
clc;clear all;close all;
[YY,MM,DD,HR,TH]= textread('jeju0510.dat','');
Dn= datenum(YY,MM,DD,HR,0,0);   % Date serial number
%  단순히 조석자료 표현
subplot(2,1,1)
plot(Dn,TH,'b');        % 파란색 점으로 표현
datetick('x',6);        % 'mm/dd'
xx=TH-mean(TH);
% FFT를 표현
subplot(2,1,2)
fftplot(xx)
xlabel('cph')