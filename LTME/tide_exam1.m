% m-file: tide_exam1.m
% Z0는 평균해면값 4.6197
% 
clc;clear all;close all;
t=0:0.5:24         % 2일 간 표현
MSL=4.6197*100     % 평균해면(Mean Sea Level)
% m2는 w:28.9841, A:2.8627, G:130.76
wm2=28.9841*pi/180
am2=2.8627*100;     % cm단위로 변환
gm2=130.76*pi/180   % phase를 Radian으로 변환 
m2=am2*cos(wm2*t-gm2);
% S2는 w:30.000, A:1.1480, G:187.62 
ws2=30.000*pi/180;          % s2의 각속도
as2=1.1480*100;
gs2=187.62*pi/180;    
s2=as2*cos(ws2*t-gs2);
% k1은 w:15.0410686, k1 A:0.3968, G:302.88 
wk1=15.0411*pi/180;ak1=0.3968*100;gk1=302.88*pi/180
k1=ak1*cos(wk1*t-gk1);
% O1은 w:13.9430, A:0.2902, G:263.58
wo1=13.9430*pi/180;ao1=0.2902*100;go1=263.58*pi/180
o1=ao1*cos(wo1*t-go1);
% 4대 분조의 합
TH=m2+s2+k1+o1;
% 그래프로 표현
plot(t,m2,'r',t,s2,'g',t,k1,'b',t,o1,'m');hold on
plot(t,TH,'k-','LineWidth',2);hold on
plot([0,25],[0,0],'k:');   % plot([x1,x2],[y1,y2]) x축표현
plot([12,12],[-600,600],'k:');
legend('M2','S2','O1','K1','Total height')
set(gca,'xLim',[0 24]);set(gca,'XTick',[0:1:24]);
ylim([-600,600]);
xlabel('time(h)');ylabel('elevation') 