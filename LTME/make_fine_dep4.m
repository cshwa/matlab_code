
close all; clear all; clc;

Xi=[127.5700:0.0005:128.1900]; % 0.0005 - 100m resolution
Yi=[34.5920:0.0005:35.1000];

load gjb_dep_deg_1970.dat
lon=gjb_dep_deg_1970(:,1);
lat=gjb_dep_deg_1970(:,2);
z=gjb_dep_deg_1970(:,3);
[Y X]=meshgrid(Yi,Xi);
% [lat lon]=meshgrid(lat1,lon1);
% h=griddata(lon,lat,z,X,Y);
h=kriging(lon,lat,z,X,Y);

save('interp_dep_1970_fine_fix.mat','h');

v = variogram([lon lat],z);
pcolor(X,Y,h); colorbar; shading flat;


close all; clear all; clc;
Xi=[127.5700:0.0005:128.1900]; % 0.0005 - 100m resolution
Yi=[34.5920:0.0005:35.1000];
load gjb_dep_deg.mat
% lon1=double(ncread('depth_sumjin_v11.nc','lon'));
% lat1=double(ncread('depth_sumjin_v11.nc','lat'));
% z=double(ncread('depth_sumjin_v11.nc','z'));


[Y X]=meshgrid(Yi,Xi);
% [lat lon]=meshgrid(lat1,lon1);
h=griddata(lon,lat,z,X,Y);

save('interp_dep_present_fine2.mat','h');
close all; clear;


close all; clear; clc;
load interp_dep_present_fine2.mat
Xi=[127.5700:0.0005:128.1900]; % 0.0005 - 100m resolution
Yi=[34.5920:0.0005:35.1000];
[Y X]=meshgrid(Yi,Xi);


h(h >= 0)=NaN;

pcolor(X,Y,h*-1); colorbar; shading flat;

