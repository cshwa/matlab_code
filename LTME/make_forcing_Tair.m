clc;close all;clear all;

for year=2012
grdfile='../02_input_depth/result/grid_h_v01.nc';

rawdata=[num2str(year),'-2.nc'];
nc1=netcdf(rawdata);
yy=yeardays(year);

longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
time=[0.5:1:yy-0.5];
Tair=nc1{'t2m'}(:);
Tsfactor=nc1{'t2m'}.scale_factor(:);
Taoffset=nc1{'t2m'}.add_offset(:);
size_T=size(Tair);
close(nc1);
Tair=Tair.*Tsfactor+Taoffset-273.15;


T_daily=ones(yy,size_T(2),size_T(3))*NaN;

for i=1:1:yy;
    aa=sum(Tair(1+4*(i-1):4*i,:,:));
    bb=sum(Tair((yy*4+1)+4*(i-1):(yy*4)+4*i,:,:));
    T_daily(i,:,:)=(aa+bb)/8;
end

clear Tair nc1;

grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_T=size(lat_rho);
inp_T_daily=ones(yy,size_inp_T(1),size_inp_T(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:yy;
    Zi=squeeze(T_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_T_daily(ii,:,:)=Z;
end

clear T_daily grd Z lat_u lon_u Xi Yi Zi

fname=['frc_ecmwf_',num2str(year),'Tair_',num2str(yy),'.nc'];
% fname=['easttest_',num2str(year),'Tair_',num2str(yy),'.nc'];
create_roms_forcing_T(fname,inp_T_daily,time,yy,yy)
end