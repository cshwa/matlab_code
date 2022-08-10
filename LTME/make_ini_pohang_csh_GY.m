clc;close all;clear all;

title='sumjin_initial_from_cshwa';

%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%

ininame=['gy_1980_tsbio_ini.nc'];
grdname=['grid_sumjin_v1970_fix_3m.nc'];

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

create_inifile_no3_nh4_do_chl_csh(ininame,grdname,title,...
               theta_s,theta_b,hc,N,...
               1,1,tini,'clobber');
           


%%% load file        
nc=netcdf('1980_ini_his_8785.nc','w');
temp=nc{'temp'}(:);
salt=nc{'salt'}(:);
u=nc{'u'}(:);
v=nc{'v'}(:);
ub=nc{'ubar'}(:);
vb=nc{'vbar'}(:);
close(nc);

load('KOEM_interped_climate.mat','do','nh4','no3','chla');

%make it input form and unit
DO=permute(squeeze(do(:,:,:,1)),[3 2 1])*0.7*44.661;  %% mg/L to millimole_oxygen meter-3
NH4=permute(squeeze(nh4(:,:,:,1)),[3 2 1])/1000*1000/14; %% ug/L to millimole_N meter-3
NO3=permute(squeeze(no3(:,:,:,1)),[3 2 1])/1000*1000/14; %% ug/L to millimole_N meter-3
Chla=permute(squeeze(chla(:,:,:,1)),[3 2 1]); %% ug/L == mg/m^3

nc1=netcdf(ininame,'w');
nc1{'temp'}(:) = temp;
nc1{'salt'}(:) = salt;
nc1{'oxygen'}(:) = DO;
nc1{'NH4'}(:) = NH4;
nc1{'NO3'}(:) = NO3;
nc1{'chlorophyll'}(:) = Chla;
nc1{'u'}(:) = u;
nc1{'v'}(:) = v;
nc1{'ubar'}(:) = ub;
nc1{'vbar'}(:) = vb;
close(nc1);
                     
return

%%% confirm

close all; clear; clc;   % -v2
ininame=['gy_1980_tsbio_ini.nc'];
grdname=['grid_sumjin_v1970_fix_3m.nc'];
 
lon=ncread(grdname,'lon_rho');
lat=ncread(grdname,'lat_rho');
lonu=ncread(grdname,'lon_u');
latu=ncread(grdname,'lat_u');
mask=ncread(grdname,'mask_rho');
masku=ncread(grdname,'mask_u');
salt=ncread(ininame,'salt');
temp=ncread(ininame,'temp');
nh=ncread(ininame,'NH4');
no=ncread(ininame,'NO3');
do=ncread(ininame,'oxygen');
chl=ncread(ininame,'chlorophyll');
u=ncread(ininame,'u');
v=ncread(ininame,'v');
ub=ncread(ininame,'ubar');
vb=ncread(ininame,'vbar');

figure; pcolor(lon,lat,squeeze(temp(:,:,20)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'^oC');
grid on;

figure; pcolor(lon,lat,squeeze(salt(:,:,20)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'PSU');
grid on;

figure; pcolor(lon,lat,squeeze(nh(:,:,20)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'mM/m^3');
grid on;

figure; pcolor(lon,lat,squeeze(no(:,:,20)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'mM/m^3');
grid on;

figure; pcolor(lon,lat,squeeze(chl(:,:,20)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'mg/m^3');
grid on;

figure; pcolor(lon,lat,squeeze(do(:,:,20)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'mg/m^3');
grid on;


figure; pcolor(lonu,latu,squeeze(u(:,:,20)).*(masku./masku)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'ms^-^1');
grid on;

figure; pcolor(lonu,latu,ub.*(masku./masku)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'ms^-^1');
grid on;

