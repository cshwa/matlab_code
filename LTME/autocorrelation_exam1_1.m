% autocorrelation  활용
% autocorrel_exam1.m
% input: busantemp_1996.dat
clear all;close all;clc;
[yy,mm,dd,temp]=textread('busantemp_1996.dat','','headerlines',1);
n=length(yy); 
subplot(2,1,1)
% 전체자료를 적용해서 구하기, Lags는 'Series길이보다 작아야 한다.(length - 1)
% Bounds=2는 95 percent confidence interval
[tempACF,Lags,Bounds] = autocorr(temp,n-1,2); 
autocorr(temp,n-1,2);
title('부산 조위관측소 수온자료의 autocorrelation 결과');
ylim([-1, 1]);

subplot(2,1,2)
plot(Lags,tempACF,'b');
title('autocorrelation 결과 Display')
xlabel('Lag');ylabel('Temperature ACF')