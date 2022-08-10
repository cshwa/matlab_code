close all; clear; clc; 
lon=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m_v3.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix_3m_v3.nc','h');

figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on; caxis([0 40]);

%present
lonp=ncread('ocean_avg_0001.nc','lon_rho');
latp=ncread('ocean_avg_0001.nc','lat_rho');
maskp=ncread('ocean_avg_0001.nc','mask_rho');
hp=ncread('ocean_avg_0001.nc','h');

figure; pcolor(lonp,latp,hp.*(maskp./maskp)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on; caxis([0 40]);

hp_mask = hp.*(maskp./maskp);
hp_mask(find(isnan(hp_mask)==1)) = 0;

figure; colormap(jet); pcolor(lon,lat,hp_mask - h.*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on; caxis([-5 0]);xlim([127.58 127.86]); ylim([34.77 35]);
