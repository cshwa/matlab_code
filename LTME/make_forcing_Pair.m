clc;close all;clear all;

for year=2012
grdfile='../02_input_depth/result/grid_h_v01.nc';

rawdata=[num2str(year),'-2.nc'];
nc1=netcdf(rawdata);
yy=yeardays(year);

longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
time=[0.5:1:yy-0.5];
Pair=nc1{'msl'}(:);
Psfactor=nc1{'msl'}.scale_factor(:);
Paoffset=nc1{'msl'}.add_offset(:);
size_P=size(Pair);
close(nc1);
Pair=(Pair.*Psfactor+Paoffset)./100;


P_daily=ones(yy,size_P(2),size_P(3))*NaN;

for i=1:1:yy;
    aa=sum(Pair(1+4*(i-1):4*i,:,:));
    bb=sum(Pair((yy*4+1)+4*(i-1):(yy*4)+4*i,:,:));
    P_daily(i,:,:)=(aa+bb)/8;
end

clear Pair nc1;



grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_P=size(lat_rho);
inp_P_daily=ones(yy,size_inp_P(1),size_inp_P(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:yy;
    Zi=squeeze(P_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_P_daily(ii,:,:)=Z;
end

clear P_daily grd Z lat_u lon_ Xi Yi Zi

fname=['frc_ecmwf_',num2str(year),'Pair_',num2str(yy),'.nc'];
% fname=['easttest_',num2str(year),'Pair_',num2str(yy),'.nc'];
create_roms_forcing_P(fname,inp_P_daily,time,yy,yy)
end