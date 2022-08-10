% autocorrelation  Ȱ��
% autocorrel_exam1.m
% input: busantemp_1996.dat
clear all;close all;clc;
[yy,mm,dd,temp]=textread('busantemp_1996.dat','','headerlines',1);
n=length(yy); 
subplot(2,1,1)
% ��ü�ڷḦ �����ؼ� ���ϱ�, Lags�� 'Series���̺��� �۾ƾ� �Ѵ�.(length - 1)
% Bounds=2�� 95 percent confidence interval
[tempACF,Lags,Bounds] = autocorr(temp,n-1,2); 
autocorr(temp,n-1,2);
title('�λ� ���������� �����ڷ��� autocorrelation ���');
ylim([-1, 1]);

subplot(2,1,2)
plot(Lags,tempACF,'b');
title('autocorrelation ��� Display')
xlabel('Lag');ylabel('Temperature ACF')