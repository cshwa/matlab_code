%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS boundary file
%
%  Extrapole and interpole temperature and salinity from a
%  climatology to get boundary conditions for
%  ROMS (boundary netcdf file) .
%  Get the velocities and sea surface elevation via a 
%  geostrophic computation.
%
%  Data input format (netcdf):
%     temperature(T, Z, Y, X)
%     T : time [Months]
%     Z : Depth [m]
%     Y : Latitude [degree north]
%     X : Longitude [degree east]
%
%  Data source : IRI/LDEO climate Data Library (World Ocean Atlas 1998)
%    http://ingrid.ldgo.columbia.edu/
%    http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NODC/.WOA98/
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    1-Sep-2006 by Pierrick Penven
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
for i = 3:3
    yeart=i;
% bryname=['bndy_auto_NWP_1_10_test6_clim_',num2str(i),'regime.nc'];
bryname=['bndy_auto_NWP_1_10_test6_clim_16to20.nc']; % 3regime
grdname=['D:\장기생태\Dynamic\result\2013\Input\grid_gy_v11_s.nc'];
makebry=1;
makeZbry=0;
makepisces=0;
makeplot=1;
insitu2pot=1;   %1: transform in-situ temperature to potential temperature
frcname=['sumjin_v11_tide_larged2_csh_v2.nc'];
ROMS_title  = 'sumjin_OBC_from_cshwa';
% zref = -10;
theta_s = 1;
theta_b = 1;
hc      =1;
N=20;

obc = [1 1 0 0]; % open boundaries (1=open , [S E N W])

%
%  Data climatologies file names:
%
%    temp_month_data : monthly temperature climatology
%    temp_ann_data   : annual temperature climatology
%    salt_month_data : monthly salinity climatology
%    salt_ann_data   : annual salinity climatology

bry_time=[15:30:365];
bry_cycle=365; cycle = 365;
%% from make_biovar_boundary_and_ini_06to15_v2.m
load(['D:\장기생태\Dynamic\KOEM\boundary_bio_input_to06to15_',num2str(i),'regime_v7.mat']);
do_tem=do_input;
do_sn = permute(do_tem,[3 2 1]);
clearvars do_tem
do_tem=do_input_we;
do_we = permute(do_tem,[3 2 1]);

no3_sn=permute(no3_input,[3 2 1]);
no3_we =permute(no3_input_we,[3 2 1]);

po4_sn=permute(po4_input,[3 2 1]);
po4_we =permute(po4_input_we,[3 2 1]);

nh4_sn=permute(nh4_input_sn,[3 2 1]);
nh4_we =permute(nh4_input_we,[3 2 1]);

chl_sn=permute(chl_input_sn,[3 2 1]);
chl_we =permute(chl_input_we,[3 2 1]);

nc1=netcdf(bryname,'w');

nc1('NO3_time') = length(bry_time);
nc1{'NO3_time'} = ncdouble('NO3_time') ;
nc1{'NO3_time'}.long_name = ncchar('time for Nitrate');
nc1{'NO3_time'}.long_name = 'time for Nitrate';
nc1{'NO3_time'}.units = ncchar('day');
nc1{'NO3_time'}.units = 'day';
nc1{'NO3_time'}.calendar = ncchar('360.0 days in every year');
nc1{'NO3_time'}.calendar = '360.0 days in every year';
nc1{'NO3_time'}.cycle_length = cycle;
nc1{'NO3_time'}(:)=bry_time;

nc1{'NO3_east'} = ncdouble('NO3_time','s_rho','eta_rho') ;
nc1{'NO3_east'}.long_name = ncchar('eastern boundary Nitrate');
nc1{'NO3_east'}.long_name = 'northern boundary Nitrate';
nc1{'NO3_east'}.units = ncchar('mMol N m-3');
nc1{'NO3_east'}.units = 'mMol N m-3';
nc1{'NO3_east'}.coordinates = ncchar('lon_rho s_rho NO3_time');
nc1{'NO3_east'}.coordinates = 'lon_rho s_rho NO3_time';
% nc1{'NO3_east'}(:)=permute(d3_no3_e_t,[3 2 1]);
nc1{'NO3_east'}(:)=no3_we;

nc1{'NO3_south'} = ncdouble('NO3_time','s_rho','xi_rho') ;
nc1{'NO3_south'}.long_name = ncchar('southern boundary Nitrate');
nc1{'NO3_south'}.long_name = 'southern boundary Nitrate';
nc1{'NO3_south'}.units = ncchar('mMol N m-3');
nc1{'NO3_south'}.units = 'mMol N m-3';
nc1{'NO3_south'}.coordinates = ncchar('lon_rho s_rho NO3_time');
nc1{'NO3_south'}.coordinates = 'lon_rho s_rho NO3_time';
% nc1{'NO3_south'}(:)=permute(d3_no3_s_t,[3 2 1]);
nc1{'NO3_south'}(:)=no3_sn;

nc1('NH4_time') = length(bry_time);
nc1{'NH4_time'} = ncdouble('NH4_time') ;
nc1{'NH4_time'}.long_name = ncchar('time for Ammonium');
nc1{'NH4_time'}.long_name = 'time for Ammonium';
nc1{'NH4_time'}.units = ncchar('day');
nc1{'NH4_time'}.units = 'day';
nc1{'NH4_time'}.calendar = ncchar('360.0 days in every year');
nc1{'NH4_time'}.calendar = '360.0 days in every year';
nc1{'NH4_time'}.cycle_length = cycle;
nc1{'NH4_time'}(:)=bry_time;

nc1{'NH4_east'} = ncdouble('NH4_time','s_rho','eta_rho') ;
nc1{'NH4_east'}.long_name = ncchar('eastern boundary Ammonium');
nc1{'NH4_east'}.long_name = 'eastern boundary Ammonium';
nc1{'NH4_east'}.units = ncchar('mMol N m-3');
nc1{'NH4_east'}.units = 'mMol N m-3';
nc1{'NH4_east'}.coordinates = ncchar('lon_rho s_rho NH4_time');
nc1{'NH4_east'}.coordinates = 'lon_rho s_rho NH4_time';
% nc1{'NH4_east'}(:)=permute(d3_nh4_e_t,[3 2 1]);
nc1{'NH4_east'}(:)=nh4_we;

nc1{'NH4_south'} = ncdouble('NH4_time','s_rho','xi_rho') ;
nc1{'NH4_south'}.long_name = ncchar('southern boundary Ammonium');
nc1{'NH4_south'}.long_name = 'southern boundary Ammonium';
nc1{'NH4_south'}.units = ncchar('mMol N m-3');
nc1{'NH4_south'}.units = 'mMol N m-3';
nc1{'NH4_south'}.coordinates = ncchar('lon_rho s_rho NH4_time');
nc1{'NH4_south'}.coordinates = 'lon_rho s_rho NH4_time';
% nc1{'NH4_south'}(:)=permute(d3_nh4_s_t,[3 2 1]);
nc1{'NH4_south'}(:)=nh4_sn;

nc1('tPO4_time') = length(bry_time);
nc1{'tPO4_time'} = ncdouble('tPO4_time') ;
nc1{'tPO4_time'}.long_name = ncchar('time for Phosphate');
nc1{'tPO4_time'}.long_name = 'time for Phosphate';
nc1{'tPO4_time'}.units = ncchar('day');
nc1{'tPO4_time'}.units = 'day';
nc1{'tPO4_time'}.calendar = ncchar('360.0 days in every year');
nc1{'tPO4_time'}.calendar = '360.0 days in every year';
nc1{'tPO4_time'}.cycle_length = cycle;
nc1{'tPO4_time'}(:)=bry_time;

nc1{'tPO4_south'} = ncdouble('tPO4_time','s_rho','xi_rho') ;
nc1{'tPO4_south'}.long_name = ncchar('southern boundary Phosphate');
nc1{'tPO4_south'}.long_name = 'southern boundary Phosphate';
nc1{'tPO4_south'}.units = ncchar('milimole PO4 m-3');
nc1{'tPO4_south'}.units = 'milimole PO4 m-3';
nc1{'tPO4_south'}.coordinates = ncchar('lon_rho s_rho tPO4_time');
nc1{'tPO4_south'}.coordinates = 'lon_rho s_rho tPO4_time';
nc1{'tPO4_south'}(:)=po4_sn;

nc1{'tPO4_east'} = ncdouble('tPO4_time','s_rho','eta_rho') ;
nc1{'tPO4_east'}.long_name = ncchar('eastern boundary Phosphate');
nc1{'tPO4_east'}.long_name = 'eastern boundary Phosphate';
nc1{'tPO4_east'}.units = ncchar('milimole PO4 m-3');
nc1{'tPO4_east'}.units = 'milimole PO4 m-3';
nc1{'tPO4_east'}.coordinates = ncchar('lon_rho s_rho tPO4_time');
nc1{'tPO4_east'}.coordinates = 'lon_rho s_rho tPO4_time';
nc1{'tPO4_east'}(:)=po4_we;

nc1('oxygen_time') = length(bry_time);
nc1{'oxygen_time'} = ncdouble('oxygen_time') ;
nc1{'oxygen_time'}.long_name = ncchar('time for oxygen');
nc1{'oxygen_time'}.long_name = 'time for oxygen';
nc1{'oxygen_time'}.units = ncchar('day');
nc1{'oxygen_time'}.units = 'day';
nc1{'oxygen_time'}.calendar = ncchar('360.0 days in every year');
nc1{'oxygen_time'}.calendar = '360.0 days in every year';
nc1{'oxygen_time'}.cycle_length = cycle;
nc1{'oxygen_time'}(:)=bry_time;
% 
% have to oxygen tobe close
nc1{'oxygen_south'} = ncdouble('oxygen_time','s_rho','xi_rho') ;
nc1{'oxygen_south'}.long_name = ncchar('southern boundary oxygen');
nc1{'oxygen_south'}.long_name = 'southern boundary oxygen';
nc1{'oxygen_south'}.units = ncchar('mg/L');
nc1{'oxygen_south'}.units = 'mg/L';
nc1{'oxygen_south'}.coordinates = ncchar('lon_rho s_rho oxygen_time');
nc1{'oxygen_south'}.coordinates = 'lon_rho s_rho oxygen_time';
% nc1{'oxygen_south'}(:)=permute(d3_nh4_s_t.*0,[3 2 1]);
nc1{'oxygen_south'}(:)=do_sn

nc1{'oxygen_east'} = ncdouble('oxygen_time','s_rho','eta_rho') ;
nc1{'oxygen_east'}.long_name = ncchar('eastern boundary oxygen');
nc1{'oxygen_east'}.long_name = 'eastern boundary oxygen';
nc1{'oxygen_east'}.units = ncchar('mg/L');
nc1{'oxygen_east'}.units = 'mg/L';
nc1{'oxygen_east'}.coordinates = ncchar('lon_rho s_rho oxygen_time');
nc1{'oxygen_east'}.coordinates = 'lon_rho s_rho oxygen_time';
% nc1{'oxygen_east'}(:)=permute(d3_nh4_s_t.*0,[3 2 1]);
nc1{'oxygen_east'}(:)=do_we;

nc1('chlo_time') = length(bry_time);
nc1{'chlo_time'} = ncdouble('chlo_time') ;
nc1{'chlo_time'}.long_name = ncchar('time for chlo');
nc1{'chlo_time'}.long_name = 'time for chlo';
nc1{'chlo_time'}.units = ncchar('day');
nc1{'chlo_time'}.units = 'day';
nc1{'chlo_time'}.calendar = ncchar('360.0 days in every year');
nc1{'chlo_time'}.calendar = '360.0 days in every year';
nc1{'chlo_time'}.cycle_length = cycle;
nc1{'chlo_time'}(:)=bry_time;

nc1{'chlo_east'} = ncdouble('chlo_time','s_rho','eta_rho') ;
nc1{'chlo_east'}.long_name = ncchar('eastern boundary chlo');
nc1{'chlo_east'}.long_name = 'eastern boundary chlo';
nc1{'chlo_east'}.units = ncchar('mg/L');
nc1{'chlo_east'}.units = 'mg/L';
nc1{'chlo_east'}.coordinates = ncchar('lon_rho s_rho chlo_time');
nc1{'chlo_east'}.coordinates = 'lon_rho s_rho chlo_time';
% nc1{'chlo_east'}(:)=permute(d3_chl_e_t,[3 2 1]);
nc1{'chlo_east'}(:)=chl_we;

nc1{'chlo_south'} = ncdouble('chlo_time','s_rho','xi_rho') ;
nc1{'chlo_south'}.long_name = ncchar('southern boundary chlo');
nc1{'chlo_south'}.long_name = 'southern boundary chlo';
nc1{'chlo_south'}.units = ncchar('mg/L');
nc1{'chlo_south'}.units = 'mg/L';
nc1{'chlo_south'}.coordinates = ncchar('lon_rho s_rho chlo_time');
nc1{'chlo_south'}.coordinates = 'lon_rho s_rho chlo_time';
% nc1{'chlo_south'}(:)=permute(d3_chl_s_t,[3 2 1]);
nc1{'chlo_south'}(:)=chl_sn;

% 
% nc1('zoop_time') = length(bry_time);
% nc1{'zoop_time'} = ncdouble('zoop_time') ;
% nc1{'zoop_time'}.long_name = ncchar('time for zoop');
% nc1{'zoop_time'}.long_name = 'time for zoop';
% nc1{'zoop_time'}.units = ncchar('day');
% nc1{'zoop_time'}.units = 'day';
% nc1{'zoop_time'}.calendar = ncchar('360.0 days in every year');
% nc1{'zoop_time'}.calendar = '360.0 days in every year';
% nc1{'zoop_time'}.cycle_length = cycle;
% nc1{'zoop_time'}(:)=time;
% 
% nc1{'zoop_north'} = ncdouble('zoop_time','s_rho','xi_rho') ;
% nc1{'zoop_north'}.long_name = ncchar('northern boundary zoop');
% nc1{'zoop_north'}.long_name = 'northern boundary zoop';
% nc1{'zoop_north'}.units = ncchar('mg/L');
% nc1{'zoop_north'}.units = 'mg/L';
% nc1{'zoop_north'}.coordinates = ncchar('lon_rho s_rho zoop_time');
% nc1{'zoop_north'}.coordinates = 'lon_rho s_rho zoop_time';
% nc1{'zoop_north'}(:)=0;
% 
% nc1('phyt_time') = length(bry_time);
% nc1{'phyt_time'} = ncdouble('phyt_time') ;
% nc1{'phyt_time'}.long_name = ncchar('time for phyt');
% nc1{'phyt_time'}.long_name = 'time for phyt';
% nc1{'phyt_time'}.units = ncchar('day');
% nc1{'phyt_time'}.units = 'day';
% nc1{'phyt_time'}.calendar = ncchar('360.0 days in every year');
% nc1{'phyt_time'}.calendar = '360.0 days in every year';
% nc1{'phyt_time'}.cycle_length = cycle;
% nc1{'phyt_time'}(:)=time;
% 
% nc1{'phyt_north'} = ncdouble('phyt_time','s_rho','xi_rho') ;
% nc1{'phyt_north'}.long_name = ncchar('northern boundary phyt');
% nc1{'phyt_north'}.long_name = 'northern boundary phyt';
% nc1{'phyt_north'}.units = ncchar('mg/L');
% nc1{'phyt_north'}.units = 'mg/L';
% nc1{'phyt_north'}.coordinates = ncchar('lon_rho s_rho phyt_time');
% nc1{'phyt_north'}.coordinates = 'lon_rho s_rho phyt_time';
% nc1{'phyt_north'}(:)=0;
% 
% nc1('LDeN_time') = length(bry_time);
% nc1{'LDeN_time'} = ncdouble('LDeN_time') ;
% nc1{'LDeN_time'}.long_name = ncchar('time for LDeN');
% nc1{'LDeN_time'}.long_name = 'time for LDeN';
% nc1{'LDeN_time'}.units = ncchar('day');
% nc1{'LDeN_time'}.units = 'day';
% nc1{'LDeN_time'}.calendar = ncchar('360.0 days in every year');
% nc1{'LDeN_time'}.calendar = '360.0 days in every year';
% nc1{'LDeN_time'}.cycle_length = cycle;
% nc1{'LDeN_time'}(:)=time;
% 
% nc1{'LDeN_north'} = ncdouble('LDeN_time','s_rho','xi_rho') ;
% nc1{'LDeN_north'}.long_name = ncchar('northern boundary LDeN');
% nc1{'LDeN_north'}.long_name = 'northern boundary LDeN';
% nc1{'LDeN_north'}.units = ncchar('mg/L');
% nc1{'LDeN_north'}.units = 'mg/L';
% nc1{'LDeN_north'}.coordinates = ncchar('lon_rho s_rho LDeN_time');
% nc1{'LDeN_north'}.coordinates = 'lon_rho s_rho LDeN_time';
% nc1{'LDeN_north'}(:)=0;
% 
% nc1('SDeN_time') = length(bry_time);
% nc1{'SDeN_time'} = ncdouble('SDeN_time') ;
% nc1{'SDeN_time'}.long_name = ncchar('time for SDeN');
% nc1{'SDeN_time'}.long_name = 'time for SDeN';
% nc1{'SDeN_time'}.units = ncchar('day');
% nc1{'SDeN_time'}.units = 'day';
% nc1{'SDeN_time'}.calendar = ncchar('360.0 days in every year');
% nc1{'SDeN_time'}.calendar = '360.0 days in every year';
% nc1{'SDeN_time'}.cycle_length = cycle;
% nc1{'SDeN_time'}(:)=time;
% 
% nc1{'SDeN_north'} = ncdouble('SDeN_time','s_rho','xi_rho') ;
% nc1{'SDeN_north'}.long_name = ncchar('northern boundary SDeN');
% nc1{'SDeN_north'}.long_name = 'northern boundary SDeN';
% nc1{'SDeN_north'}.units = ncchar('mg/L');
% nc1{'SDeN_north'}.units = 'mg/L';
% nc1{'SDeN_north'}.coordinates = ncchar('lon_rho s_rho SDeN_time');
% nc1{'SDeN_north'}.coordinates = 'lon_rho s_rho SDeN_time';
% nc1{'SDeN_north'}(:)=0;

close(nc1);
                     
clearvars -except i
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
