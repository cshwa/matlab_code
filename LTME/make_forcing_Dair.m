clc;close all;clear all;

for year=2012
grdfile='../02_grid_depth/smoothing/grid_gy_v11_s.nc';

rawdata=[num2str(year),'-1.nc'];
nc1=netcdf(rawdata);
yy=yeardays(year);

longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
time=[0.5:1:yy-0.5];
Dair=nc1{'d2m'}(:);
Dsfactor=nc1{'d2m'}.scale_factor(:);
Daoffset=nc1{'d2m'}.add_offset(:);
size_D=size(Dair);
close(nc1);
Dair=Dair.*Dsfactor+Daoffset-273.15;


D_daily=ones(yy,size_D(2),size_D(3))*NaN;

for i=1:1:yy;
    aa=sum(Dair(1+4*(i-1):4*i,:,:));
    bb=sum(Dair((yy*4+1)+4*(i-1):(yy*4)+4*i,:,:));
    D_daily(i,:,:)=(aa+bb)/8;
end

clear Dair nc1;

grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_D=size(lat_rho);
inp_D_daily=ones(yy,size_inp_D(1),size_inp_D(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:yy;
    Zi=squeeze(D_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_D_daily(ii,:,:)=Z;
end

clear D_daily grd Z lat_u lon_u Xi Yi Zi

fname=['frc_ecmwf_',num2str(year),'Dair_',num2str(yy),'.nc'];
% fname=['easttest_',num2str(year),'Dair_',num2str(yy),'.nc'];
create_roms_forcing_D(fname,inp_D_daily,time,yy,yy)
end