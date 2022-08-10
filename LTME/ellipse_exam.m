% 4대 분조값 조류타원도
% m-file: ellipse_exam.m
% 
clc;clear all;close all
figure(1)
% M2 분조의 조류타원도 작성
Sma = 24.141;Smi = 1.174;D = 33;
ellipse(Sma,Smi,D*pi/180,0,0,'r');
hold on
%S2 분조의 조류타원도 작성
Sma = 10.104;Smi = 0.003;D=30.9;
ellipse(Sma,Smi,D*pi/180,0,0,'b'); hold on
plot(0,0,'r.');  % y축
xlim([-30 30]);ylim([-30 30]);axis equal
title('M_2 & S_2 Tidal current ellipse')
legend('M_2','S_2','Location','SouthEast');legend('boxoff')

% K1분조의 조류타원도 작성
figure(2)                    
Sma = 9.518;Smi = 0.504;D=39.8;
ellipse(Sma,Smi,D*pi/180,0,0,'m');      % 장축,단축,각도,중심좌표(0,0)
hold on
% O1 분조의 조류타원도 작성
Sma = 8.834;Smi = 0.890;D=46.9;
ellipse(Sma,Smi,D*pi/180,0,0,'g');     % 장축,단축,각도,중심좌표(0,0)
plot(0,0,'r.')
xlim([-30 30]);ylim([-30 30]);axis equal
title('k_1 & O_1 Tidal current ellipse');
legend('K_1','O_1','Location','SouthEast');legend('boxoff')