clc; clear all; close all

%% 과거 자료에서 새로 매립된 곳 추가하여 수심 변경

load merge_1970_depth_gjb_fix.mat
lon=X; lat=Y; z=h_merge_1;
z(z<0)=0; % land(z<0) to be zero
X=ncread('grid_sumjin_v1970.nc','lon_rho');
Y=ncread('grid_sumjin_v1970.nc','lat_rho');

Z=griddata(lon,lat,z,X,Y);
[a b] = size(Z);

Z(isnan(Z)==1)=0;
% pcolor(X,Y,Z);shading flat;

%%
nw=netcdf('depth_sumjin_v1970_csh2.nc','clobber');
nw('Lon') = length(X(:,1));
nw('Lat') = length(Y(1,:));
nw('one') = 1;

nw{'lon'} = ncfloat('one','Lon');
nw{'lon'}.long_name = ncchar('Longitude');
nw{'lon'}.units = ncchar('degree_east');
nw{'lon'}(:)=X(:,1);

nw{'lat'} = ncfloat('one','Lat');
nw{'lat'}.long_name = ncchar('Latitude');
nw{'lat'}.units = ncchar('degree_north');
nw{'lat'}(:)=Y(1,:);

nw{'z'} = ncfloat('Lat','Lon');
nw{'z'}.long_name = ncchar('z');
nw{'z'}.units = ncchar('meter');
nw{'z'}(:)=Z;
% nw{'z'}(:)=10;

close(nw);



