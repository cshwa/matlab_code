% using SEAWATER Tool box and Function making simple T?S plotting program 
% 2004-4-22, by Jeff Book (U.S. Naval Resarch Laboratory) and PETER
% INPUT : CTD_TED1.dat, CTD_TED2.dat, CTD_TED3.dat
% m?file: TS_diagram2.m
% It need function file ('ts_plot_general.m'), and seawater toolbox
% 
clc;clear all;close all
clc;clear all;close all 
%load CTD_TED1.dat
data1=load('CTD_TED1.dat');
w(:,2) = data1(:,3);      % Temp data
w(:,3) = data1(:,7);      % Salinity data
n = size(w,1);
% CTD_TED2.dat DATA load
%load CTD_TED2.dat
data2=load('CTD_TED3.dat');
w(n+1:length(data2)+n,2) = data2(:,3);   % TED1의 Temp Data에 이어서 
w(n+1:length(data2)+n,3) = data2(:,7);   % TED2의 Salinity Data에 이어서
n = size(w,1);
% CTD_TED3.dat DATA load
% data3=load('CTD_TED3.dat');
% w(n+1:length(data3)+n,2) = data3(:,3); % TED1&TED2의 Temp Data에 이어서
% w(n+1:length(data3)+n,3) = data3(:,7); % TED1&TED2의 Salinity Data에 이어서
ylim;             % Temp Variation display
xlim;             % Salinity Variation display
% Call function 'ts_plot_general.m'
[c,h,x,y,dens_grid] = ts_plot_general(w,'r');  
clabel(c)                                       % Density 표시 
title('T-S Diagram')
set(gcf,'Color','w')               % 배경색 흰색