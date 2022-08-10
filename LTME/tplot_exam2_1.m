% plot tide_exam1
% Z0는 평균해면값  4.6197


% O1 A:0.2902, G:263.58
% k1 A:0.3968, G:302.88 
% N2 A:0.5335, G:111.10
% K2 A:0.3169, G:183.49
clc;clear all;close all;
t=0:0.5:48         % 2일 간 표현
msl=4.6197*100     % 평균해면(Mean Sea Level)
% m2 표현 w:28.9841, A:2.8627, G:130.76
am2=2.8627*100;     % cm단위로 변환
gm2=130.76*pi/180   % phase를 Radian으로 변환 
m2=am2*cos(28.9841*pi/180*t-gm2);
% S2 표현, w:30.000, A:1.1480, G:187.62 
as2=1.1480*100
gs2=187.62*pi/180    
s2=as2*cos(
% o1
o1=29.02*cos(13.943*pi/180*t-263.58*pi/180);
o1a=29.02*cos(14.492*pi/180*t-263.58*pi/180);
xx=m2+o1;
xx1=m2+o1a;

subplot(2,1,1)
%plot(t,xx,'r');hold on
plot(t,xx1,'b')
legend('m2+O1','change O1')
set(gca,'xLim',[0 744]);
set(gca,'XTick',[0:24:744]);
ylim([-400,400]);

subplot(2,1,2)
xx1=m2+3*o1
plot(t,xx1,'r')
legend('m2+O1*3')
set(gca,'xLim',[0 744]),
set(gca,'XTick',[0:24:744])
ylim([-400 400])