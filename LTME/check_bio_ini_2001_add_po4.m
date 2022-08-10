close all; clear; clc; 
do=ncread('monthly_2016_01_fix_po4.nc','oxygen');
nh4=ncread('monthly_2016_01_fix_po4.nc','NH4');
no3=ncread('monthly_2016_01_fix_po4.nc','NO3');
chl=ncread('monthly_2016_01_fix_po4.nc','chlorophyll');
phy=ncread('monthly_2016_01_fix_po4.nc','phytoplankton');
zoo=ncread('monthly_2016_01_fix_po4.nc','zooplankton');
detS=ncread('monthly_2016_01_fix_po4.nc','SdetritusN');
detL=ncread('monthly_2016_01_fix_po4.nc','LdetritusN');
lon=ncread('monthly_2016_01_fix_po4.nc','lon_rho');
lat=ncread('monthly_2016_01_fix_po4.nc','lat_rho');
mask=ncread('monthly_2016_01_fix_po4.nc','mask_rho');

zoo(zoo > 0.2) = 0.2;

figure; pcolor(lon,lat,squeeze(zoo(:,:,20)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'zoo (umol/m^3)');
grid on; caxis([0 0.2]);

nc=netcdf('mon_2016_01_fix_p_z_02.nc','w');
nc{'zooplankton'}(:) = permute(zoo,[3 2 1]);
close(nc);

close all; clear; clc; 
zoo=ncread('mon_2016_01_fix_p_z_02.nc','zooplankton');
lon=ncread('mon_2016_01_fix_p_z_02.nc','lon_rho');
lat=ncread('mon_2016_01_fix_p_z_02.nc','lat_rho');
mask=ncread('mon_2016_01_fix_p_z_02.nc','mask_rho');

figure; pcolor(lon,lat,squeeze(zoo(:,:,20)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'zoo (umol/m^3)');
grid on; caxis([0 0.5]);


find(do <= 0)
find(nh4 <= 0)
find(no3 <= 0)
find(chl <= 0)
find(zoo <= 0)
find(detS <= 0)
find(detL <= 0)


find(isnan(do) == 1)
find(isnan(nh4) == 1)
find(isnan(no3) == 1)
find(isnan(chl) == 1)
find(isnan(zoo) == 1)
find(isnan(detS) == 1)
find(isnan(detL) == 1)

close all; clear; clc; 
load('D:\장기생태\Dynamic\KOEM\boundary_bio_input_2001.mat','po4_input*');
po4_ini=nanmean(squeeze(po4_input(:,:,1)),1)

xnum=size(po4_input,1)
ynum=size(po4_input_we,1)

for i = 1:20
po4_ini_3d(20-i+1,:,:) = po4_ini(20-i+1) .* ones(ynum,xnum);
end

file_n = 'monthly_2016_01_fix_po4.nc'
nc1=netcdf(file_n,'w');
nc1{'tPO4'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
nc1{'tPO4'}.long_name = ncchar('Phosphate');
nc1{'tPO4'}.units = ncchar('milimole PO4 m-');
nc1{'tPO4'}.fields = ncchar('tPO4, scalar, series');
nc1{'tPO4'}.missing_value=ncfloat(1.0000e+37);
nc1{'tPO4'}(:)= po4_ini_3d;

nc1{'SdetritusP'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
nc1{'SdetritusP'}.long_name = ncchar('small fraction nitrogen detritus concentration');
nc1{'SdetritusP'}.units = ncchar('mMol N m-3');
nc1{'SdetritusP'}.fields = ncchar('SdetritusP, scalar, series');
nc1{'SdetritusP'}.missing_value=ncfloat(1.0000e+37);
nc1{'SdetritusP'}(:)=0;

nc1{'LdetritusP'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
nc1{'LdetritusP'}.long_name = ncchar('large fraction nitrogen detritus concentration');
nc1{'LdetritusP'}.units = ncchar('mMol N m-3');
nc1{'LdetritusP'}.fields = ncchar('LdetritusP, scalar, series');
nc1{'LdetritusP'}.missing_value=ncfloat(1.0000e+37);
nc1{'LdetritusP'}(:)=0;

close(nc1);

%check
tPO4=ncread('monthly_2016_01_fix_po4.nc','tPO4');
SdetritusP=ncread('monthly_2016_01_fix_po4.nc','SdetritusP');
LdetritusP=ncread('monthly_2016_01_fix_po4.nc','LdetritusP');

find(SdetritusP ~= 0)
find(LdetritusP ~= 0)

