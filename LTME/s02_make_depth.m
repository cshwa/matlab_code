clc; clear all; close all

%% 과거 자료에서 새로 매립된 곳 추가하여 수심 변경
data=load('..\01_2dm(basemap)\depth\GY_depth_silver.dat');

lon=data(:,1);
lat=data(:,2);
z=data(:,3);

% Xi=[min(lon)-0.005:(max(lon)-min(lon))/49:max(lon)+0.005]; % 0.0005 - 100m resolution
% Yi=[min(lat)-0.005:(35.2-34.853)/199:max(lat)+0.005];
% 
% Xi=[127.7394-0.005:0.0005:127.81+0.005]; % 0.0005 - 100m resolution
% Yi=[34.853-0.005:0.0005:35.2+0.005];

Xi = lon; Yi = lat;
[X,Y]=meshgrid(Xi,Yi);
Z=griddata(lon,lat,z,X,Y);
[a b] = size(Z);

Z(isnan(Z)==1)=0;
pcolor(X,Y,Z);shading flat;

%%
nw=netcdf('depth_sumjin_v11.nc','clobber');
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



