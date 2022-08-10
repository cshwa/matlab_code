% plot tide_exam1

% m2 2003년 결과 A:2.8627, G:130.76
% S2 2003년 결과 A:1.1480, G:187.62
% O1 2003년 결과 A:0.2902, G:263.58
% k1 2003년 결과 A:0.3968, G:302.88 
%
clc;clear all;close all;
t=0:0.5:24*31;             % 31일 표현
m2=286.27*cos(28.984*pi/180*t-130.76*pi/180);
% plot(t,m2);hold on
o1=29.02*cos(13.943*pi/180*t-263.58*pi/180);
k1=29.02*cos(15.041*pi/180*t-302.88*pi/180);

xx=m2+k1+o1;
xx1=m2+2*k1+2*o1;

subplot(2,1,1)
%plot(t,xx,'r');hold on
plot(t,xx1,'b')
legend('M2+K1+O1')
set(gca,'xLim',[0 744]);
set(gca,'XTick',[0:24:744]);
ylim([-400,400]);

subplot(2,1,2)
plot(t,xx1,'r')
legend('M2+2*K1+2*O1')
set(gca,'xLim',[0 744]),
set(gca,'XTick',[0:24:744])
ylim([-400 400])