% ������ �ط��ڷ� ǥ�� 
% m-file: curmap_exam1.m
% 
clc;clear all;close all;
m_proj('mercator','lon',[128 137],'lat',[33 39.5]);  % ����, ��������
m_gshhs_i('color','k');            % �ذ��ػ� �����ؾȼ� ������
m_gshhs_i('patch',[.8 .8 .8]);       % �������� ȸ������ ����
m_grid('box','fancy','tickdir','in','linewidth',1);
m_text(130.5,38.5,'EAST SEA','vertical','top');   %������ ������ ǥ���Ѵ�.
hold on
%curtest.dat������ �ҷ��� ���� �� �Ϲ漺������ ����
data = load('curtest.dat')
long = data(:,1);lat = data(:,2);speed = data(:,3);
direction=450-data(:,4);
v=speed.*sin(direction.*pi./180); % ���������� ������ �Ϲ漺�� ���
u=speed.*cos(direction.*pi./180); % ���������� ������ ���漺�� ���
m_quiver(long,lat,u,v);           % �ط� Vector ǥ��
m_plot(long,lat,'ro')
set(gcf,'Color','w');              % ���� White 
