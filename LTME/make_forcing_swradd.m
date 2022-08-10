clc;close all;clear all;

for year=2012
grdfile='../02_input_depth/result/grid_h_v01.nc';
nc1=netcdf(rawdata);
yy=yeardays(year);

longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
time=[0.5:1:yy-0.5];
swrad=nc1{'ssrd'}(:);
SWsfactor=nc1{'ssrd'}.scale_factor(:);
SWaoffset=nc1{'ssrd'}.add_offset(:);
size_SW=size(swrad);
close(nc1);
swrad=(swrad.*SWsfactor+SWaoffset);

SW_daily=ones(yy,size_SW(2),size_SW(3))*NaN;

for i=1:1:yy;
    aa=sum(swrad(1+4*(i-1):4*i,:,:));
    bb=sum(swrad((yy*4+1)+4*(i-1):(yy*4)+4*i,:,:));
    SW_daily(i,:,:)=(aa+bb)/8;
end

clear swrad nc1;

grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_SW=size(lat_rho);
inp_SW_daily=ones(yy,size_inp_SW(1),size_inp_SW(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:yy;
    Zi=squeeze(SW_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_SW_daily(ii,:,:)=Z;
end

clear SW_daily grd Z lat_u lon_u Xi Yi Zi
fname=['frc_ecmwf_',num2str(year),'swradd_',num2str(yy),'.nc'];
% fname=['easttest_',num2str(year),'Tair_',num2str(yy),'.nc'];
create_roms_forcing_SWD(fname,inp_T_daily,time,yy,yy)
end