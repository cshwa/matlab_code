clc;close all;clear all;
cd /data1/temp/ECMWF_interim/u10
for year=2001:2010
clearvars -except year
grdfile='/home/cshwa/Long/02_grid_depth/smoothing/grid_gy_v11_s.nc';
% grdfile='/home/cshwa/SAR/sar_grid.nc';


rawdata=['/data1/temp/ECMWF_interim/u10/ECMWF_Interim_u10_',num2str(year),'.nc'];
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
cd /home/cshwa/Long
save([num2str(year),'_u10.mat']);
% save('2014_u10.mat','-v7.3');
end

close all; clear; clc;
% for year = 2001:2009
for year = 2010:2010
    clearvars -except year
    load([num2str(year),'_u10.mat']);
varname = 'u10'
grdfile='/home/cshwa/Long/02_grid_depth/smoothing/grid_gy_v11_s.nc';
tinterval=1;
% load([num2str(year),'_u10.mat']);
data_daily=permute(inp_U_daily,[3 2 1]);
surface_forcing_ECMWF(grdfile, year, varname, tinterval,data_daily)
% create_roms_forcing_U(fname,inp_U_daily,time,yy,yy)
end
