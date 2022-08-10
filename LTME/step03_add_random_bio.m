clc;clear all;close all

filename='.\NWP_1_10_test6\data\roms_bndy_auto_NWP_2003_fennel_woa13.nc';

nc=netcdf(filename,'write');
time=nc{'temp_time'}(:);
cycle=nc{'temp_time'}.cycle_length(:);

model_gridfile='g:\auto_fennel\grid\roms_grd_auto_rdrg2_new8_smooth.nc';
ng=netcdf(model_gridfile,'r');
lon_rho=ng{'lon_rho'}(:);
N=40; % layer
[y,x]=size(lon_rho);
close(ng);

% NO3_north=nc{'NO3_north'}(:);
% NO3_south=nc{'NO3_south'}(:);
% NO3_east=nc{'NO3_east'}(:);
% NO3_west=nc{'NO3_west'}(:);
% NO3_north(NO3_north<0)=0;
% NO3_south(NO3_south<0)=0;
% NO3_eat(NO3_east<0)=0;
% NO3_west(NO3_west<0)=0;
% 
% nc{'NO3_south'}(:)=NO3_north;
% nc{'NO3_north'}(:)=NO3_south;
% nc{'NO3_east'}(:)=NO3_east;
% nc{'NO3_west'}(:)=NO3_west;

% oxygen_north=nc{'oxygen_north'}(:);
oxygen_south=nc{'oxygen_south'}(:);
oxygen_east=nc{'oxygen_east'}(:);
% oxygen_west=nc{'oxygen_west'}(:);
% oxygen_north=oxygen_north*44.661;
oxygen_south=oxygen_south*44.661;
oxygen_east=oxygen_east*44.661;
% oxygen_west=oxygen_west*44.661;
% chlo_north=nc{'chlo_north'}(:);
chlo_south=nc{'chlo_south'}(:);
chlo_east=nc{'chlo_east'}(:);
% chlo_west=nc{'chlo_west'}(:);

nc{'oxygen_south'}(:)=oxygen_south;
% nc{'oxygen_north'}(:)=oxygen_north;
nc{'oxygen_east'}(:)=oxygen_east;
% nc{'oxygen_west'}(:)=oxygen_west;


nc{'NH4_south'}(:)=0;
% nc{'NH4_north'}(:)=1;
nc{'NH4_east'}(:)=0;
% nc{'NH4_west'}(:)=1;

nc('zoop_time') = length(time);

nc{'zoop_time'} = ncdouble('zoop_time') ;
nc{'zoop_time'}.long_name = ncchar('time for zoop');
nc{'zoop_time'}.units = ncchar('day');
nc{'zoop_time'}.calendar = ncchar('360.0 days in every year');
nc{'zoop_time'}.cycle_length = cycle;
nc{'zoop_time'}(:)=time;

nc{'zoop_south'} = ncdouble('zoop_time','s_rho','xi_rho') ;
nc{'zoop_south'}.long_name = ncchar('southern boundary zooplankton');
nc{'zoop_south'}.units = ncchar('millimole nitrogen meter-3');
nc{'zoop_south'}.coordinates = ncchar('zoop_time s_rho lon_rho');
nc{'zoop_south'}(:)=chlo_south*0.2;

nc{'zoop_east'} = ncdouble('zoop_time','s_rho','eta_rho') ;
nc{'zoop_east'}.long_name = ncchar('eastern boundary zooplankton');
nc{'zoop_east'}.units = ncchar('millimole nitrogen meter-3');
nc{'zoop_east'}.coordinates = ncchar('zoop_time s_rho lat_rho');
nc{'zoop_east'}(:)=chlo_east*0.2;
% 
% nc{'zoop_north'} = ncdouble('zoop_time','s_rho','xi_rho') ;
% nc{'zoop_north'}.long_name = ncchar('northern boundary zooplankton');
% nc{'zoop_north'}.units = ncchar('millimole nitrogen meter-3');
% nc{'zoop_north'}.coordinates = ncchar('zoop_time s_rho lon_rho');
% nc{'zoop_north'}(:)=chlo_north*0.2;
% 
% nc{'zoop_west'} = ncdouble('zoop_time','s_rho','eta_rho') ;
% nc{'zoop_west'}.long_name = ncchar('western boundary zooplankton');
% nc{'zoop_west'}.units = ncchar('millimole nitrogen meter-3');
% nc{'zoop_west'}.coordinates = ncchar('zoop_time s_rho lat_rho');
% nc{'zoop_west'}(:)=chlo_west*0.2;


nc('phyt_time') = length(time);

nc{'phyt_time'} = ncdouble('phyt_time') ;
nc{'phyt_time'}.long_name = ncchar('time for phyt');
nc{'phyt_time'}.units = ncchar('day');
nc{'phyt_time'}.calendar = ncchar('360.0 days in every year');
nc{'phyt_time'}.cycle_length = cycle;
nc{'phyt_time'}(:)=time;

nc{'phyt_south'} = ncdouble('phyt_time','s_rho','xi_rho') ;
nc{'phyt_south'}.long_name = ncchar('southern boundary phytoplankton');
nc{'phyt_south'}.units = ncchar('millimole nitrogen meter-3');
nc{'phyt_south'}.coordinates = ncchar('phyt_time s_rho lon_rho');
nc{'phyt_south'}(:)=chlo_south*0.5;

nc{'phyt_east'} = ncdouble('phyt_time','s_rho','eta_rho') ;
nc{'phyt_east'}.long_name = ncchar('eastern boundary phytoplankton');
nc{'phyt_east'}.units = ncchar('millimole nitrogen meter-3');
nc{'phyt_east'}.coordinates = ncchar('phyt_time s_rho lat_rho');
nc{'phyt_east'}(:)=chlo_east*0.5;

% nc{'phyt_north'} = ncdouble('phyt_time','s_rho','xi_rho') ;
% nc{'phyt_north'}.long_name = ncchar('northern boundary phytoplankton');
% nc{'phyt_north'}.units = ncchar('millimole nitrogen meter-3');
% nc{'phyt_north'}.coordinates = ncchar('phyt_time s_rho lon_rho');
% nc{'phyt_north'}(:)=chlo_north*0.5;
% 
% nc{'phyt_west'} = ncdouble('phyt_time','s_rho','eta_rho') ;
% nc{'phyt_west'}.long_name = ncchar('western boundary phytoplankton');
% nc{'phyt_west'}.units = ncchar('millimole nitrogen meter-3');
% nc{'phyt_west'}.coordinates = ncchar('phyt_time s_rho lat_rho');
% nc{'phyt_west'}(:)=chlo_west*0.5;

nc('LDeN_time') = length(time);

nc{'LDeN_time'} = ncdouble('LDeN_time') ;
nc{'LDeN_time'}.long_name = ncchar('time for LDeN');
nc{'LDeN_time'}.units = ncchar('day');
nc{'LDeN_time'}.calendar = ncchar('360.0 days in every year');
nc{'LDeN_time'}.cycle_length = cycle;
nc{'LDeN_time'}(:)=time;

nc{'LDeN_south'} = ncdouble('LDeN_time','s_rho','xi_rho') ;
nc{'LDeN_south'}.long_name = ncchar('southern boundary LDeN');
nc{'LDeN_south'}.units = ncchar('millimole nitrogen meter-3');
nc{'LDeN_south'}.coordinates = ncchar('LDeN_time s_rho lon_rho');
nc{'LDeN_south'}(:)=chlo_south*0.5*0.35/2;

nc{'LDeN_east'} = ncdouble('LDeN_time','s_rho','eta_rho') ;
nc{'LDeN_east'}.long_name = ncchar('eastern boundary LDeN');
nc{'LDeN_east'}.units = ncchar('millimole nitrogen meter-3');
nc{'LDeN_east'}.coordinates = ncchar('LDeN_time s_rho lat_rho');
nc{'LDeN_east'}(:)=chlo_east*0.5*0.35/2;

nc('LDeP_time') = length(time);

nc{'LDeP_time'} = ncdouble('LDeP_time') ;
nc{'LDeP_time'}.long_name = ncchar('time for LDeP');
nc{'LDeP_time'}.units = ncchar('day');
nc{'LDeP_time'}.calendar = ncchar('360.0 days in every year');
nc{'LDeP_time'}.cycle_length = cycle;
nc{'LDeP_time'}(:)=time;

nc{'LDeP_south'} = ncdouble('LDeP_time','s_rho','xi_rho') ;
nc{'LDeP_south'}.long_name = ncchar('southern boundary LDeP');
nc{'LDeP_south'}.units = ncchar('millimole phosphorus meter-3');
nc{'LDeP_south'}.coordinates = ncchar('LDeP_time s_rho lon_rho');
nc{'LDeP_south'}(:)=chlo_south*0.5*0.35/2;

nc{'LDeP_east'} = ncdouble('LDeP_time','s_rho','eta_rho') ;
nc{'LDeP_east'}.long_name = ncchar('eastern boundary LDeP');
nc{'LDeP_east'}.units = ncchar('millimole phosphorus meter-3');
nc{'LDeP_east'}.coordinates = ncchar('LDeP_time s_rho lat_rho');
nc{'LDeP_east'}(:)=chlo_east*0.5*0.35/2;

% nc{'LDeN_north'} = ncdouble('LDeN_time','s_rho','xi_rho') ;
% nc{'LDeN_north'}.long_name = ncchar('northern boundary LDeN');
% nc{'LDeN_north'}.units = ncchar('millimole nitrogen meter-3');
% nc{'LDeN_north'}.coordinates = ncchar('LDeN_time s_rho lon_rho');
% nc{'LDeN_north'}(:)=chlo_north*0.5*0.35;
% 
% 
% nc{'LDeN_west'} = ncdouble('LDeN_time','s_rho','eta_rho') ;
% nc{'LDeN_west'}.long_name = ncchar('western boundary LDeN');
% nc{'LDeN_west'}.units = ncchar('millimole nitrogen meter-3');
% nc{'LDeN_west'}.coordinates = ncchar('LDeN_time s_rho lat_rho');
% nc{'LDeN_west'}(:)=chlo_west*0.5*0.35;

nc('SDeN_time') = length(time);

nc{'SDeN_time'} = ncdouble('SDeN_time') ;
nc{'SDeN_time'}.long_name = ncchar('time for SDeN');
nc{'SDeN_time'}.units = ncchar('day');
nc{'SDeN_time'}.calendar = ncchar('360.0 days in every year');
nc{'SDeN_time'}.cycle_length = cycle;
nc{'SDeN_time'}(:)=time;

nc{'SDeN_south'} = ncdouble('SDeN_time','s_rho','xi_rho') ;
nc{'SDeN_south'}.long_name = ncchar('southern boundary SDeN');
nc{'SDeN_south'}.units = ncchar('millimole nitrogen meter-3');
nc{'SDeN_south'}.coordinates = ncchar('SDeN_time s_rho lon_rho');
nc{'SDeN_south'}(:)=chlo_south*0.5*0.35/2;

nc{'SDeN_east'} = ncdouble('SDeN_time','s_rho','eta_rho') ;
nc{'SDeN_east'}.long_name = ncchar('eastern boundary SDeN');
nc{'SDeN_east'}.units = ncchar('millimole nitrogen meter-3');
nc{'SDeN_east'}.coordinates = ncchar('SDeN_time s_rho lat_rho');
nc{'SDeN_east'}(:)=chlo_east*0.5*0.35/2;

nc('SDeP_time') = length(time);

nc{'SDeP_time'} = ncdouble('SDeP_time') ;
nc{'SDeP_time'}.long_name = ncchar('time for SDeP');
nc{'SDeP_time'}.units = ncchar('day');
nc{'SDeP_time'}.calendar = ncchar('360.0 days in every year');
nc{'SDeP_time'}.cycle_length = cycle;
nc{'SDeP_time'}(:)=time;

nc{'SDeP_south'} = ncdouble('SDeP_time','s_rho','xi_rho') ;
nc{'SDeP_south'}.long_name = ncchar('southern boundary SDeP');
nc{'SDeP_south'}.units = ncchar('millimole phosphate meter-3');
nc{'SDeP_south'}.coordinates = ncchar('SDeP_time s_rho lon_rho');
nc{'SDeP_south'}(:)=chlo_south*0.5*0.35/2;

nc{'SDeP_east'} = ncdouble('SDeP_time','s_rho','eta_rho') ;
nc{'SDeP_east'}.long_name = ncchar('eastern boundary SDeP');
nc{'SDeP_east'}.units = ncchar('millimole phosphate meter-3');
nc{'SDeP_east'}.coordinates = ncchar('SDeP_time s_rho lat_rho');
nc{'SDeP_east'}(:)=chlo_east*0.5*0.35/2;

% nc{'SDeN_north'} = ncdouble('SDeN_time','s_rho','xi_rho') ;
% nc{'SDeN_north'}.long_name = ncchar('northern boundary SDeN');
% nc{'SDeN_north'}.units = ncchar('millimole nitrogen meter-3');
% nc{'SDeN_north'}.coordinates = ncchar('SDeN_time s_rho lon_rho');
% nc{'SDeN_north'}(:)=chlo_north*0.5*0.35;
% 
% nc{'SDeN_west'} = ncdouble('SDeN_time','s_rho','eta_rho') ;
% nc{'SDeN_west'}.long_name = ncchar('western boundary SDeN');
% nc{'SDeN_west'}.units = ncchar('millimole nitrogen meter-3');
% nc{'SDeN_west'}.coordinates = ncchar('SDeN_time s_rho lat_rho');
% nc{'SDeN_west'}(:)=chlo_west*0.5*0.35;
close(nc);