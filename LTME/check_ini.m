close all; clear; clc; 
lon=ncread('grid_sumjin_v1970.nc','lon_rho');
lat=ncread('grid_sumjin_v1970.nc','lat_rho');
mask=ncread('grid_sumjin_v1970.nc','mask_rho');
h=ncread('grid_sumjin_v1970.nc','h');

figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on;

figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading interp;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on;


close all; clear; clc;   % -v2
lon=ncread('grid_sumjin_v1970_fix.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix.nc','h');
salt=ncread('gy_v11_s_ini_ts_const.nc','salt');
temp=ncread('gy_v11_s_ini_ts_const.nc','temp');
u=ncread('gy_v11_s_ini_ts_const.nc','u');
v=ncread('gy_v11_s_ini_ts_const.nc','v');

figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on;

figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading interp;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on;