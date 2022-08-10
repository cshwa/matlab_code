clc; clear all; close all

%% 과거 자료에서 새로 매립된 곳 추가하여 수심 변경
% data=load('.\80년대지형\gjb_dep_deg_1970.dat');
load gjb_dep_deg.mat

% lon=data(:,1);
% lat=data(:,2);
% z=data(:,3);

% Xi=[min(lon)-0.005:(max(lon)-min(lon))/49:max(lon)+0.005]; % 0.0005 - 100m resolution
% Yi=[min(lat)-0.005:(35.2-34.853)/199:max(lat)+0.005];

Xi=[127.5700:0.0005:128.1900]; % 0.0005 - 100m resolution
Yi=[34.5920:0.0005:35.1];
% 
% [xi yi]=meshgrid(Yi,Xi);
% 
% z_inter=griddata(lon, lat, z, xi,yi);

% Xi = lon; Yi = lat;
[X,Y]=meshgrid(Xi,Yi);

% X=ncread('grid_sumjin_v11.nc','lon_rho');
% Y=ncread('grid_sumjin_v11.nc','lat_rho');

Z=griddata(lon,lat,z,X,Y);
[a b] = size(Z);

% Z(isnan(Z)==1)=0;
% pcolor(X,Y,Z);shading flat;

%%
nw=netcdf('depth_sumjin_present_csh2.nc','clobber');
nw('Lon') = length(Xi);
nw('Lat') = length(Yi);
nw('one') = 1;

nw{'lon'} = ncfloat('one','Lon');
nw{'lon'}.long_name = ncchar('Longitude');
nw{'lon'}.units = ncchar('degree_east');
nw{'lon'}(:)=Xi;

nw{'lat'} = ncfloat('one','Lat');
nw{'lat'}.long_name = ncchar('Latitude');
nw{'lat'}.units = ncchar('degree_north');
nw{'lat'}(:)=Yi;

nw{'z'} = ncfloat('Lat','Lon');
nw{'z'}.long_name = ncchar('z');
nw{'z'}.units = ncchar('meter');
nw{'z'}(:)=Z;
% nw{'z'}(:)=10;

close(nw);



