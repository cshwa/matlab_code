% EARDO 3D MAP using Matlab 
% input: scotra.dat
% M-file: map3d_exam.m
% 

clc;clear all;close all;
[long,lat,alt]=textread('scotra.dat','%n%n%n');
longmin = (min(long));longmax = (max(long));
latmin =  (min(lat));latmax = (max(lat));
Xp = [longmin:.0001:longmax]; 
Yp = [latmin:.0001:latmax]; 
[Xi,Yi] = meshgrid(Xp,Yp); 
Zi = griddata(long,lat,-alt,Xi,Yi); 
meshc(Xi,Yi,Zi) 
shading interp % or colormap jet or shading flat 
view(-45,45)
colorscale([2 0],[-60 0],100,'vert','position',[0.92 0.3 0.01 3]);
colorbar('vert')
set(gca,'YaxisLocation','right');
title('3D 이어도 해저수심도')