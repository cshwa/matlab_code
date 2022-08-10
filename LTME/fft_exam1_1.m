% m-file: fft_exam1.m
% ������ �Լ��� ����� �� ����� FFT�����ϱ�
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
% Periodogram �ۼ�
power = abs(fu(1:floor(n/2))).^2;
nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist;
plot(freq,power)
xlabel('cycles/year');title('Periodogram');

subplot(2,2,4)
% �ִ밪 ǥ��
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
