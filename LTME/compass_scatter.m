% Scatterplot�� Compass grid�� ǥ���ϱ�
% m-file: compass_scatter.m
% input: suyeong9805.dat

clc;clear all;close all
% data load
[yy,mm,dd,hr,min,temp,speed,n_dir,u_comp,v_comp]=textread('suyeong9805.dat','');
% ��� compas �׸���
polar(1,100);hold on  % compass�� ���ɿ��� ������ 100
% scatterplot ǥ��
plot(u_comp,v_comp,'b.','MarkerSize',5)   % marksize�� 5
set(gcf,'Color','w');       % ���� White