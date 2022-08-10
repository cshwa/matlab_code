clc;close all;clear all;

for year=2014
grdfile='E:\SAR_data\sar_grid.nc';


rawdata=['E:\SAR_data\ECMWF_Interim_v10_',num2str(year),'.nc'];
nc1=netcdf(rawdata);
yy=yeardays(year);
longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
time=[0.5:1:yy-0.5];
V_wind=nc1{'v10'}(:);
Vsfactor=nc1{'v10'}.scale_factor(:);
Vaoffset=nc1{'v10'}.add_offset(:);
size_V=size(V_wind);
close(nc1);
V_wind=V_wind.*Vsfactor+Vaoffset;



V_daily=ones(yy,size_V(2),size_V(3))*NaN;

for i=1:1:yy;
    aa=sum(V_wind(1+4*(i-1):4*i,:,:));
    bb=sum(V_wind((yy*4+1)+4*(i-1):(yy*4)+4*i,:,:));
    V_daily(i,:,:)=(aa+bb)/8;
end

clear V_wind nc1;



grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_V=size(lat_rho);
inp_V_daily=ones(yy,size_inp_V(1),size_inp_V(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:yy;
    Zi=squeeze(V_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_V_daily(ii,:,:)=Z;
end

clear V_daily grd Z lat_v lon_v Xi Yi Zi

fname=['frc_ecmwf_',num2str(year),'Vwind_',num2str(yy),'.nc'];
% fname=['easttest_',num2str(year),'Vwind_',num2str(yy),'.nc'];
create_roms_forcing_V(fname,inp_V_daily,time,yy,yy)
end