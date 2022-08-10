% m-file: fft_exam1.m
% 임의의 함수를 만들고 그 결과를 FFT실행하기
%
clc;clear all;close all

subplot(2,2,[1:2])
[year,spotnum]=textread('sunspot_tong.dat','');
plot(year,spotnum);title('Sunspot Data')
n=length(year);
tt=1:n;
fu=fft(spotnum); 
fu(1)=[];

subplot(2,2,3)
% Periodogram 작성
power = abs(fu(1:floor(n/2))).^2;
nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist;
plot(freq,power)
xlabel('cycles/year');title('Periodogram');

subplot(2,2,4)
% 최대값 표현
period=1./freq;
plot(period,power);
axis([0 40 0 2e+7]);
ylabel('Power');xlabel('period (years/cycle)');
hold on;
index=find(power==max(power));
mainPeriodStr=num2str(period(index));
plot(period(index),power(index),'r.', 'MarkerSize',25);
text(period(index)+2,power(index),['Period = ',mainPeriodStr]);
hold off
