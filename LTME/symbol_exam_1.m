% ������ �ɺ�ǥ��
% symbol_exam.m
%
clc;clear all;close all;
m_proj('mercator','lon',[128 134],'lat',[33 39]);  % �浵128~134, ���� 33~39�����ǵ��� �ۼ�
m_gshhs_h('patch',[1 1 .7]);                       % ���� ������ ����
m_grid('linestyle','none','tickdir','in','linewidth',1)  % ���浵�� ǥ��
m_text(130.5,38.5,'EAST SEA','vertical','top');          % ������ �������� ǥ��
m_text(131.8,34.55,'JAPAN','vertical','top');m_text(128.8,37.5,'Donghae','vertical','top');
m_text(129.2,36.5,'Hupo','vertical','top');m_text(129.06,36.1,'Pohang','vertical','top');
m_text(129.2,35.6,'Ulsan','vertical','top');m_text(129,35.25,'Busan','vertical','top');
m_text(131.85,37.43,'Dokdo','vertical','top')
m_line(131,36,'marker','square','markersize',8,'color','r');     % �����׸� �ɺ�
m_line(130.8,37,'marker','diamond','markersize',8,'color','g');  % ��� ������ �ɺ�