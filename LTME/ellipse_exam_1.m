% 4�� ������ ����Ÿ����
% m-file: ellipse_exam.m
% 
clc;clear all;close all
figure(1)
% M2 ������ ����Ÿ���� �ۼ�
Sma = 24.141;Smi = 1.174;D = 33;
ellipse(Sma,Smi,D*pi/180,0,0,'r');
hold on
%S2 ������ ����Ÿ���� �ۼ�
Sma = 10.104;Smi = 0.003;D=30.9;
ellipse(Sma,Smi,D*pi/180,0,0,'b'); hold on
plot(0,0,'r.');  % y��
xlim([-30 30]);ylim([-30 30]);axis equal
title('M_2 & S_2 Tidal current ellipse')
legend('M_2','S_2','Location','SouthEast');legend('boxoff')

% K1������ ����Ÿ���� �ۼ�
figure(2)                    
Sma = 9.518;Smi = 0.504;D=39.8;
ellipse(Sma,Smi,D*pi/180,0,0,'m');      % ����,����,����,�߽���ǥ(0,0)
hold on
% O1 ������ ����Ÿ���� �ۼ�
Sma = 8.834;Smi = 0.890;D=46.9;
ellipse(Sma,Smi,D*pi/180,0,0,'g');     % ����,����,����,�߽���ǥ(0,0)
plot(0,0,'r.')
xlim([-30 30]);ylim([-30 30]);axis equal
title('k_1 & O_1 Tidal current ellipse');
legend('K_1','O_1','Location','SouthEast');legend('boxoff')