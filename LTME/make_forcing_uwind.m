clc;close all;clear all;

for year=2012
grdfile='../02_input_depth/result/grid_h_v01.nc';


rawdata=[num2str(year),'-1.nc'];
nc1=netcdf(rawdata);
yy=yeardays(year);
time=[0.5:1:yy-0.5];
longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
U_wind=nc1{'u10'}(:);
Usfactor=nc1{'u10'}.scale_factor(:);
Uaoffset=nc1{'u10'}.add_offset(:);
size_U=size(U_wind);
close(nc1);
U_wind=U_wind.*Usfactor+Uaoffset;


U_daily=ones(yy,size_U(2),size_U(3))*NaN;

for i=1:1:yy;
    aa=sum(U_wind(1+4*(i-1):4*i,:,:));
    bb=sum(U_wind((yy*4+1)+4*(i-1):(yy*4)+4*i,:,:));
    U_daily(i,:,:)=(aa+bb)/8;
end

clear U_wind nc1;



grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_U=size(lat_rho);
inp_U_daily=ones(yy,size_inp_U(1),size_inp_U(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:yy;
    Zi=squeeze(U_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_U_daily(ii,:,:)=Z;
end

clear U_daily grd Z lat_u lon_u Xi Yi Zi

fname=['frc_ecmwf_',num2str(year),'Uwind_',num2str(yy),'.nc'];
% fname=['easttest_',num2str(year),'Uwind_',num2str(yy),'.nc'];
create_roms_forcing_U(fname,inp_U_daily,time,yy,yy)
end