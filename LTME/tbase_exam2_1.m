% tbase�� �̿��� �츮���� �α� ���ɵ�
% tbase_exam2.m
clc;clear all;close all;
m_proj('mercator','lon',[120 140],'lat',[30 50]);
m_tbase('contourf');
% ���浵 ǥ�ù��
m_grid('linestyle','none','tickdir','in','linewidth',2);
m_grid('box','fancy','tickdir','in');