clc;clear all;close all

title='sumjin_initial';
%
% Common parameters
%
% coastaline - lon, lat load
load data\GY_coastline_lonlat.mat;  % 모델 영역의 해안선
temp_ann_data='temp_GY.nc';
salt_ann_data='salt_GY.nc';

%--- 기존 모델 결과의 TS 사용하기 -------------------------------------------
% ncload F:\ROMS\Sumjin\grid_v11\6\daily_mean\ocean_avg_0020.nc
%--------------------------------------------------------------------------
% %
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%

ininame=['gy_1970_fix_ini.nc'];     % 생성될 nc 파일 이름
grdname=['grid_sumjin_v1970_fix.nc']; % 사용될 grid 자료 - 앞 단계에서 만들어 진 것

scoord = [1 1 1 20]; 
theta_s=scoord(1);
theta_b=scoord(2);
hc=scoord(3);
N=scoord(4);
tini=0;
disp(' ')
disp([' Making initial file: ',ininame])
disp(' ')
disp([' Title: ',title])
%
% Initial file
%

create_inifile(ininame,grdname,title,...
               theta_s,theta_b,hc,N,...
               tini,'clobber');
%
% Horizontal and vertical interp/extrapolations 
%
disp(' ')
disp(' Interpolations / extrapolations')
disp(' ')
disp(' Temperature...')
ext_tracers_ini_pohang(ininame,grdname,temp_ann_data,...
            'temperature','temp',1,1,'r',tini);
disp(' ')
disp(' Salinity...')
ext_tracers_ini_pohang(ininame,grdname,salt_ann_data,...
             'salinity','salt',1,1,'r',tini);
         

%
% Geostrophy
%
%  disp(' ')
%  disp(' Compute geostrophic currents')
%  geost_currents(ininame,grdname,oaname,frcname,zref,obc,0)
%
% Initial file
%
% if (insitu2pot)
  disp(' ')
  disp(' Compute potential temperature from in-situ...')
  getpot(ininame,grdname)
% end
%%
figure(1)
hold on
plot(lon,lat,'k','linewidth',2);
nc=netcdf(ininame);
ng=netcdf(grdname);
lon_rho=ng{'lon_rho'}(:);
lat_rho=ng{'lat_rho'}(:);
temp=nc{'temp'}(:);
pcolor(lon_rho,lat_rho,squeeze(temp(end,:,:)));shading flat;colorbar;

% 
figure(2)
hold on
plot(lon,lat,'k','linewidth',2);
salt=nc{'salt'}(:);
pcolor(lon_rho,lat_rho,squeeze(salt(1,:,:)));shading flat;colorbar;



