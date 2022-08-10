% Scatterplot을 Compass grid상에 표현하기
% m-file: compass_scatter.m
% input: suyeong9805.dat

clc;clear all;close all
% data load
[yy,mm,dd,hr,min,temp,speed,n_dir,u_comp,v_comp]=textread('suyeong9805.dat','');
% 배경 compas 그리기
polar(1,100);hold on  % compass의 동심원의 반지름 100
% scatterplot 표현
plot(u_comp,v_comp,'b.','MarkerSize',5)   % marksize를 5
set(gcf,'Color','w');       % 배경색 White