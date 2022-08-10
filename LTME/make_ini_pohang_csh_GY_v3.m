clc;close all;clear all;

title='sumjin_initial_from_cshwa';

%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%

ininame=['gy_1980_tsbio_ini_v3.nc'];
grdname=['grid_sumjin_v1970_fix_3m_v3.nc'];

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

create_inifile_no3_nh4_do_chl_csh2(ininame,grdname,title,...
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

load('boundary_bio_input_for_ini.mat');

temp_pre = ones(size(temp,3),size(temp,2),size(temp,1));
for i = 1:N
    chl_pre(:,:,i) = temp_pre(:,:,i) .* chl_in(1,i); 
end

for i = 1:N
   nh4_pre(:,:,i) = temp_pre(:,:,i) .* nh4_in(1,i); 
end

no3_in_m=squeeze(mean(no3_in,1));
for i = 1:N
   no3_pre(:,:,i) = temp_pre(:,:,i) .* no3_in_m(i,1); 
end

do_in_m=squeeze(mean(do_in(:,:,1,:),1));
for i = 1:N
    do_pre(:,:,i) = temp_pre(:,:,i) .* do_in_m(i,1); 
end

do_ini = permute(do_pre,[3 2 1]);
no3_ini=permute(no3_pre,[3 2 1]);
nh4_ini =permute(nh4_pre,[3 2 1]);
chl_ini =permute(chl_pre,[3 2 1]);


%make it input form and unit

nc1=netcdf(ininame,'w');
nc1{'temp'}(:) = temp;
nc1{'salt'}(:) = salt;
nc1{'oxygen'}(:) = do_ini;
nc1{'NH4'}(:) = nh4_ini;
nc1{'NO3'}(:) = no3_ini;
nc1{'chlorophyll'}(:) = chl_ini;
nc1{'u'}(:) = u;
nc1{'v'}(:) = v;
nc1{'ubar'}(:) = ub;
nc1{'vbar'}(:) = vb;
% nc1{'phytoplankton'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho');
% nc1{'phytoplankton'}.long_name = ncchar('phytoplankton');
% nc1{'phytoplankton'}.units = ncchar('mMol N m-3');
% nc1{'phytoplankton'}.fields = ncchar('PHYTO, scalar, series');
% nc1{'phytoplankton'}(:)=0.05;
% 
% nc1{'zooplankton'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nc1{'zooplankton'}.long_name = ncchar('zooplankton');
% nc1{'zooplankton'}.units = ncchar('mMol N m-3');
% nc1{'zooplankton'}.fields = ncchar('ZOO, scalar, series');
% nc1{'zooplankton'}(:)=0.01;
% 
% nc1{'LdetritusN'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nc1{'LdetritusN'}.long_name = ncchar('large fraction nitrogen detritus concentration');
% nc1{'LdetritusN'}.units = ncchar('mMol N m-3');
% nc1{'LdetritusN'}.fields = ncchar('LdetritusN, scalar, series');
% nc1{'LdetritusN'}(:)=0;
% 
% nc1{'SdetritusN'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nc1{'SdetritusN'}.long_name = ncchar('small fraction nitrogen detritus concentration');
% nc1{'SdetritusN'}.units = ncchar('mMol N m-3');
% nc1{'SdetritusN'}.fields = ncchar('SdetritusN, scalar, series');
% nc1{'SdetritusN'}(:)=0;
% 
% nc1{'LdetritusC'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nc1{'LdetritusC'}.long_name = ncchar('large fraction carbon detritus concentration');
% nc1{'LdetritusC'}.units = ncchar('mMol C m-3');
% nc1{'LdetritusC'}.fields = ncchar('LdetritusC, scalar, series');
% nc1{'LdetritusC'}(:)=0;
% 
% nc1{'SdetritusC'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nc1{'SdetritusC'}.long_name = ncchar('small fraction small detritus concentration');
% nc1{'SdetritusC'}.units = ncchar('mMol C m-3');
% nc1{'SdetritusC'}.fields = ncchar('SdetritusC, scalar, series');
% nc1{'SdetritusC'}(:)=0;
% 
% nc1{'TIC'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nc1{'TIC'}.long_name = ncchar('total inorganic carbon');
% nc1{'TIC'}.units = ncchar('mMol Carobn m-3');
% nc1{'TIC'}.fields = ncchar('TIC, scalar, series');
% nc1{'TIC'}(:)=1;
% 
% nc1{'alkalinity'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nc1{'alkalinity'}.long_name = ncchar('total alkalinity');
% nc1{'alkalinity'}.units = ncchar('milliequivalents meter-3');
% nc1{'alkalinity'}.fields = ncchar('alkalinity, scalar, series');
% nc1{'alkalinity'}(:)=2;
close(nc1);
                     
return

%%% confirm

close all; clear; clc;   % -v2
ininame=['gy_1980_tsbio_ini_v3.nc'];
grdname=['grid_sumjin_v1970_fix_3m_v3.nc'];
 
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

figure; pcolor(lon,lat,squeeze(temp(:,:,1)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'^oC');
grid on;

figure; pcolor(lon,lat,squeeze(salt(:,:,1)).*(mask./mask)); hh=colorbar; shading flat;
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

