% Private map making and display method
% m-file: priv_map.m
%
clc;clear all;close all
sta_x = [126.4097,126.4933,126.5222];  % 4���� pointǥ��
sta_y = [37.2781,37.3183,37.3736];  
m_proj('mercator','lat',[37.2 37.7],'long',[126.2 126.8]);
[x y] = m_ll2xy(sta_x,sta_y);              % ���� ��ǥ�� �׸���ǥ�� ��ȯ
plot(x,y,'ro');hold on;
% ����� �ؾȼ� �ڷḦ �а� �׸� �׸��� �κ�
m_usercoast('icn_map.mat');

%m_grid              % ���� �� �浵�� �׸���
m_grid('linestyle','none','box','fancy','tickdir','in');
print -dpsc icn_map.ps    % Į�� Post Scrpt file �� ����