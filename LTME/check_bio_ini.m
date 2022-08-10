close all; clear; clc;   % -v2
ininame=['gy_1980_tsbio_ini_v3.nc'];
grdname=['grid_sumjin_v1970_fix_3m_v3.nc'];
 
lon=ncread(grdname,'lon_rho');
lat=ncread(grdname,'lat_rho');
lonu=ncread(grdname,'lon_u');
latu=ncread(grdname,'lat_u');
mask=ncread(grdname,'mask_rho');
masku=ncread(grdname,'mask_u');
salt=ncread(ininame,'salt');
temp=ncread(ininame,'temp');
nh=ncread(ininame,'NH4');
no=ncread(ininame,'NO3');
do=ncread(ininame,'oxygen');
chl=ncread(ininame,'chlorophyll');
u=ncread(ininame,'u');
v=ncread(ininame,'v');
ub=ncread(ininame,'ubar');
vb=ncread(ininame,'vbar');

load('boundary_bio_input_for_ini.mat');
N=20;
temp_pre = ones(size(temp,3),size(temp,2),size(temp,1));
for i = 1:N
    chl_pre(:,:,i) = temp_pre(:,:,i) .* chl_in(1,i); 
end

for i = 1:N
   nh4_pre(:,:,i) = temp_pre(:,:,i) .* nh4_in(1,i); 
end

no3_in_m=squeeze(mean(no3_in,1));
for i = 1:N
   no3_pre(:,:,i) = temp_pre(:,:,i) .* no3_in_m(i,1); 
end

do_in_m=squeeze(mean(do_in(:,:,1,:),1));
for i = 1:N
    do_pre(:,:,i) = temp_pre(:,:,i) .* do_in_m(i,1); 
end

find((no(:,:,1) - no3_in_m(1,1)) ~= 0)
find((no(:,:,20) - no3_in_m(20,1)) ~= 0)

find((nh(:,:,1) - nh4_in(1,1)) ~= 0)
find((nh(:,:,20) - nh4_in(1,20)) ~= 0)

find((chl(:,:,1) - chl_in(1,1)) ~= 0)
find((chl(:,:,20) - chl_in(1,20)) ~= 0)

find((do(:,:,1) - do_in_m(1,1)) ~= 0)
find((do(:,:,20) - do_in_m(20,1)) ~= 0)
