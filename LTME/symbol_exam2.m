% ������ grid����,symbol �� ���浵 ��ġ ����
% symbol_exam2.m
%
clc;clear all;close all;
m_proj('UTM','lon',[124 130],'lat',[31.5 35]);
m_gshhs_i('color','k'); %�ؾȼ� ����
m_gshhs_i('patch',[.7 .7 .7]); %�����κ� �帰�������
% ����ǥ�� �� ����
m_grid('box','on','tickdir','in','xaxisloc','bottom','yaxisloc','right','clip','on');
m_line(128,33,'marker','square','MarkerEdgeColor','r',...
'MarkerFaceColor','g','MarkerSize',9,'color','r');
m_line(128,33.5,'marker','p','MarkerEdgeColor','k',...
'MarkerFaceColor','r','MarkerSize',10,'color','r');
m_line(128.7,34.15,'marker','o','MarkerEdgeColor','y',...
'MarkerFaceColor','m','MarkerSize',12,'color','r');
set(gcf,'Color','w');              % ���� White 