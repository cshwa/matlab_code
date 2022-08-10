clc;close all;clear all;

for year=2002:2009
grdfile='/home/cshwa/Long/02_grid_depth/smoothing/grid_gy_v11_s.nc';

rawdata=['/data1/temp/ECMWF_interim/ssrd/ECMWF_Interim_ssrd_',num2str(year),'.nc'];
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
% 
for i=1:1:yy;
%     aa=sum(swrad(1+8*(i-1):8*i,:,:));
%     aa=sum(swrad(1+4*(i-1):4*i,:,:));
%     bb=sum(swrad((yy*4+1)+4*(i-1):(yy*4)+4*i,:,:));
%     SW_daily(i,:,:)=(aa+bb)/8;
%     SW_daily(i,:,:)=(swrad(1+(i-1)*4,:,:)+swrad(yy*4+(i-1)*4,:,:))/(24*60*60);

    SW_daily(i,:,:)=(swrad(4+(i-1)*8,:,:)+swrad(8+(i-1)*8,:,:))./(24*60*60);
%       SW_daily(i,:,:)=(aa+bb)/8;
end

% i_c=1;
% for i = 1:2:yy*2
%     SW_temp(1+(i_c-1)*2,:,:) = swrad(i,:,:)/43200;
%     SW_temp(2+(i_c-1)*2,:,:) = swrad(i+1,:,:)/43200;
%     i_c = i_c+1;
% end
% for i=1:1:yy
%     SW_daily(i,:,:)=(SW_temp(1+(i-1)*2,:,:)+SW_temp(2+(i-1)*2,:,:))/2;
% end
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
% save('2001_swradd_v2.mat','-v7.3');
save([num2str(year),'_swradd.mat'],'-v7.3');
end
% create_roms_forcing_SWD(fname,inp_T_daily,time,yy,yy)


close all; clear; clc;
for year = 2001:2010
    clearvars -except year
    load([num2str(year),'_swradd.mat']);
varname = 'ssrd'
grdfile='/home/cshwa/Long/02_grid_depth/smoothing/grid_gy_v11_s.nc';
tinterval=1;
% load([num2str(year),'_Tair.mat']);
data_daily=permute(inp_SW_daily,[3 2 1]);
surface_forcing_ECMWF(grdfile, year, varname, tinterval,data_daily)
% create_roms_forcing_U(fname,inp_U_daily,time,yy,yy)
end
