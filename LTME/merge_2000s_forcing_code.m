clc;close all;clear all;

for year=2001:2010
clearvars -except year
grdfile='/home/cshwa/Long/02_grid_depth/smoothing/grid_gy_v11_s.nc';
% grdfile='/home/cshwa/SAR/sar_grid.nc';

rawdata=['/data1/temp/ECMWF_interim/airT/ECMWF_Interim_airT_',num2str(year),'.nc'];
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
cd /home/cshwa/Long
save([num2str(year),'_Tair.mat'],'-v7.3');
end
% create_roms_forcing_T(fname,inp_T_daily,time,yy,yy)


clc;close all;clear all;
var = 'ssrd';

for year=2001:2010
clearvars -except year
grdfile='/home/cshwa/Long/02_grid_depth/smoothing/grid_gy_v11_s.nc';
% grdfile='/home/cshwa/SAR/sar_grid.nc';

rawdata=['/data1/temp/ECMWF_interim/ssrd/ECMWF_Interim_ssrd_',num2str(year),'.nc'];
%     rawdata = '/data1/temp/ECMWF_interim/ssrd/ECMWF_Interim_ssrd_2014.nc';
    nc1=netcdf(rawdata);
    ly = (( mod(year,4) == 0 & mod(year,100) ~= 0 ) | mod(year,400) == 0);
    if ly,    yy = 366;
    else     yy = 365;
        
    end
    
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
    i_c=1;
    for i=1:4:yy*4;
        SW_three(1+(i_c-1)*8,:,:)=swrad(i,:,:)/10800/4;  % 12�� ���� 00�ñ��� ���� -00��
        SW_three(3+(i_c-1)*8,:,:)=swrad(i+1,:,:)/10800/2; % ���ʿ�
        SW_three(5+(i_c-1)*8,:,:)=swrad(i+2,:,:)/10800/4; % 00�� ���� 12�ñ��� ���� -12��
        SW_three(7+(i_c-1)*8,:,:)=swrad(i+3,:,:)/10800/2; % ���ʿ�
        i_c=i_c+1;
    end
    i_c=1;
    for i=1+yy*4:4:yy*8;
        SW_three(2+(i_c-1)*8,:,:)=swrad(i,:,:)/10800/1; % ���ʿ�
        SW_three(4+(i_c-1)*8,:,:)=swrad(i+1,:,:)/10800/3; % ���ʿ�
        SW_three(6+(i_c-1)*8,:,:)=swrad(i+2,:,:)/10800/1; % ���ʿ�
        SW_three(8+(i_c-1)*8,:,:)=swrad(i+3,:,:)/10800/3; % ���ʿ�
        i_c=i_c+1;
    end
    
    for i=1:1:yy;
        %     aa=sum(swrad(1+4*(i-1):4*i,:,:));
        %     bb=sum(swrad((yy*4+1)+4*(i-1):(yy*4)+4*i,:,:));
        %     SW_daily(i,:,:)=(aa+bb)/8;
        SW_daily(i,:,:)=(SW_three(1+(i-1)*8,:,:)+SW_three(5+(i-1)*8,:,:))/2;
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
    
%     save('ECMWF_post2_2014_daily.mat','-v7.3');
%     clear SW_daily grd Z lat_u lon_u Xi Yi Zi
%     fname=[outpath,'eastsea_',num2str(year),'_swrad.nc'];
%     create_roms_forcing_SWD(fname,inp_SW_daily,time,yy,yy)

cd /home/cshwa/Long
save([num2str(year),'_ssrd.mat'],'-v7.3');
end

clc;close all;clear all;

for year=2001:2010
clearvars -except year
grdfile='/home/cshwa/Long/02_grid_depth/smoothing/grid_gy_v11_s.nc';
% grdfile='/home/cshwa/SAR/sar_grid.nc';

rawdata=['/data1/temp/ECMWF_interim/msl/ECMWF_Interim_msl_',num2str(year),'.nc'];
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
cd /home/cshwa/Long
save([num2str(year),'_msl.mat'],'-v7.3');
end


clc;close all;clear all;

for year=2001:2010
clearvars -except year
grdfile='/home/cshwa/Long/02_grid_depth/smoothing/grid_gy_v11_s.nc';
% grdfile='/home/cshwa/SAR/sar_grid.nc';
grdfile='sar_grid_1m_reduced.nc';

rawdata=['/data1/temp/ECMWF_interim/dewt/ECMWF_Interim_dewt_',num2str(year),'.nc'];
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

fname=['frc_ecmwf_',num2str(year),'Dair_',num2str(yy),'_v2.nc'];
% fname=['easttest_',num2str(year),'Dair_',num2str(yy),'.nc'];
% save('2014_dewt.mat','-v7.3');
cd /home/cshwa/Long
save([num2str(year),'_dewt.mat'],'-v7.3');
end
