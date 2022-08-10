clc;clear all;close all

temp=ncread('roms_ini_auto_fennel_woa13jan.nc','temp');
mask=ncread('../grid/roms_grd_auto_rdrg2_new4_smooth.nc','mask_rho');
temp=permute(temp,[3 2 1]);

mask=permute(mask,[2 1]);

pcolor(squeeze(temp(40,:,:)).*mask);shading flat;caxis([0 15])