% 수심 차이 구하기
clc; clear all; 
% close all;
%--- 80년대 수심, griddata -----
past = load('depth_80s_degree.dat');
lon1 = past(:,1);
lat1 = past(:,2);
mi_lon = min(lon1);
ma_lon = max(lon1);
mi_lat = min(lat1);
ma_lat = max(lat1);
Xi=[min(lon1):0.001:max(lon1)]; 
Yi=[min(lat1):0.001:max(lat1)];
[X,Y]=meshgrid(Xi,Yi);
z = -past(:,3)+2.71;
dep_past=griddata(lon1,lat1,z,X,Y);

%--- 최근 수심 -----------------
pre = load('depth_degree_silver.dat');
lon2 = pre(:,1);
lat2 = pre(:,2);
z2 = pre(:,3);
dep_pre=griddata(lon2,lat2,z2,X,Y);

%%
figure('position',[100 100 600 500])
pcolor(X,Y,dep_past);shading flat;
caxis([0 40]);
colorbar;

figure('position',[700 100 600 500])
pcolor(X,Y,dep_pre);shading flat;
caxis([0 40]);
colorbar;

%%
figure('position',[700 100 600 500])
pcolor(X,Y,dep_pre-dep_past);shading flat;
% caxis([0 10]);
colorbar;
%%
lon = reshape(X,1,[]);
lat = reshape(Y,1,[]);
ss = dep_pre-dep_past;
depth = reshape(ss,1,[]);


temp(:,1) = lon;
temp(:,2) = lat;
temp(:,3) = depth;

save diff_dep2.dat temp -ascii


