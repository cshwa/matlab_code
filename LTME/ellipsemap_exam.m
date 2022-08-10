% ����Ÿ���� �׸���
% m-file: ellipsemap_exam.m
% using m_ellipse.m
% ������ ������ȭ���� ���
% NAME      SPEED       MAJOR   MINOR   INC    G      G+     G-
% 13 M2    0.08051140     0.730   0.088   68.9   61.5  352.6  130.4  00ic-1.dat
%  6 M2    0.08051140     0.655  -0.009   41.2   50.4    9.1   91.6  02ic-2.dat

close all;clear all;clc
% �����ۼ�
m_proj('mercator','lat',[37.0 37.6],'long',[126.15 126.9]);
m_gshhs_f('color','k');hold on    % �ʰ��ػ� ����
m_gshhs_f('patch',[.8 .8 .8]);    % ������ ����
m_grid('linestyle','none','box','fancy','tickdir','in');  % ���� �� �浵�� �׸���
hold on
% m_ellipse�� Ȱ���� tidal ellipse �ۼ�
scale = 15;                                    %  ���� Scale����
maj1=0.730/scale;min1=0.088/scale;inc1=68.9;    %  ����
inc1=inc1*pi/180;
maj2=0.655/scale;min2=-0.009/scale;inc2=41.2;   % ����
inc2=inc2*pi/180;
m_ellipse(maj1,min1,inc1,126.4097,37.2781,'r'); % ����,����,����,�߽���ǥ
hold on
m_ellipse(maj2,min2,inc2,126.5222,37.3736,'r');  % ����,����,����,�߽���ǥ
set(gcf,'Color','w')