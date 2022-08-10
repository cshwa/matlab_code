% ��õ ��⸸ �αٿ��� ������ �������� ���� data�� Scatterplot ��
% m-file: scattermap_exam.m
% input file: 00IC-1.dat, 02IC-2.dat 

clc;clear all;close all
% data load 00ic-1.dat
[yy,mm,dd,hr,min,spd,n_dir,temp]=textread('00IC-1.dat','','headerlines',15);
dir = 450 - n_dir;               % ������ȯ
u1 = spd.*cos(dir.*pi./180);     % speed�� ���ϱ��ذ����� �Ϲ漺�� ���
v1 = spd.*sin(dir.*pi./180);     % speed�� ���ϱ��ذ����� ���漺�� ���
u1 = 126.4097+u1/3000            % ������ ���� 126.4097�� ��ȯ�� u1��
v1 = 37.2781+v1/3000             % ������ �浵 37.2781�� ��ȯ  v1��
% data load 02ic-2.dat
[yy,mm,dd,hr,min,spd,n_dir,temp,cond]=textread('02IC-2.dat','','headerlines',15);
dir = 450 - n_dir;                   % ������ȯ
u2 = spd.*cos(dir.*pi./180);     % speed�� ���ϱ��ذ����� �Ϲ漺�� ���
v2 = spd.*sin(dir.*pi./180);      % speed�� ���ϱ��ذ����� ���漺�� ���
u2 = 126.5222+u2/3000           % ������ ���� 126.4097�� ��ȯ�� u1��
v2 = 37.3736+v2/3000             % ������ �浵 37.2781�� ��ȯ  v1��
% ����� �ؾȼ� �ڷḦ �а� �׸� �׸��� �κ�
m_proj('mercator','lat',[37.2 37.7],'long',[126.2 126.8]);
m_usercoast('icn_map.mat');       
m_grid('linestyle','none','box','fancy','tickdir','in');  % ���� �� �浵�� �׸���
hold on
% scatterplot ǥ��
m_plot(u1,v1,'r.','MarkerSize',3)   % scatter plot 
hold on
m_plot(u2,v2,'r.','MarkerSize',3)   % scatter plot 
hold on
% �⺻ �������� �׸���
m_plot(126.4097,37.2781,'b+','MarkerSize',5)
m_plot(126.5222,37.3736,'b+','MarkerSize',5)