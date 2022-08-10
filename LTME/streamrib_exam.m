% streamribbon 작성
% m-file: streamrib_exam.m
%
clc;clear all;close all;

% data load, 최대 및 최소값 지정
load wind
xmin = min(x(:));xmax = max(x(:));
ymax = max(y(:));zmin = min(z(:));

% Add Slice Planes for Visual Context
wind_speed = sqrt(u.^2 + v.^2 + w.^2);
hsurfaces = slice(x,y,z,wind_speed,[xmin,100,xmax],ymax,zmin);
set(hsurfaces,'FaceColor','interp','EdgeColor','none');hold on;

% Add Contour Lines to the Slice Planes
hcont = contourslice(x,y,z,wind_speed,[xmin,100,xmax],ymax,zmin);
set(hcont,'EdgeColor',[.5,.5,.5],'LineWidth',.5)

% Define the Starting Points for the Stream Lines
[sx,sy,sz] = meshgrid(80,20:10:50,0:5:15);
hlines = streamline(x,y,z,u,v,w,sx,sy,sz);
set(hlines,'LineWidth',1.5,'Color','r')

% Define the view
view(3)
daspect([2,2,1]);axis tight
