close all; clear;clc;
load sumjin_roms_v11_2dm.mat;
Lonr=lon_rho;Latr=lat_rho;
% h=griddata(X,Y,z,lon_rho,lat_rho); %memory overflow
load dep_interp_1970.mat

mask=ncread('grid_sumjin_v1970.nc','mask_rho');

% h(h>0)=NaN;
correction_MSL=2.713;
figure;pcolor(lon_rho,lat_rho,h.*(mask'./mask')+correction_MSL);
shading flat;colorbar;
grid on;
ylim([34.75 35]);
xlim([127.55 127.9]);
caxis([2 35])
set(gca,'fontsize', 12, 'fontw','bold');
whitebg([210/255 180/255 140/255])

h_ori=ncread('grid_sumjin_v11_bf_s.nc','h');
mask_ori=ncread('grid_sumjin_v11_bf_s.nc','mask_rho');
figure;pcolor(lon_rho,lat_rho,h_ori'.*(mask_ori'./mask_ori'));
shading flat;colorbar;
grid on;
ylim([34.75 35]);
xlim([127.55 127.9]);
caxis([2 35])
set(gca,'fontsize', 12, 'fontw','bold');
whitebg([210/255 180/255 140/255])

diff = h_ori'-(h+correction_MSL);

diff(diff<0)=NaN
figure;
pcolor(lon_rho,lat_rho,diff.*(mask'./mask'));
shading flat;colorbar;
grid on;
ylim([34.75 35]);
xlim([127.55 127.9]);
set(gca,'fontsize', 12, 'fontw','bold');
% set(h1, 'color', [184/255 134/255 011/255]); 
whitebg([210/255 180/255 140/255])
 
