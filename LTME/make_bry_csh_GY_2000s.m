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
for i = 2001:2001
    yeart=i;
bryname=['bndy_',num2str(yeart),'_Gwangyang.nc'];
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

ss = cumsum(eomday(yeart,1:12));
cycle=ss(end);

bry_time=[15:30:cycle];
bry_cycle=cycle;

%
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%
disp(' ')
disp([' Making the file: ',bryname])
disp(' ')
disp([' Title: pohang'])
%
% Read in the grid
%
disp(' ')
disp(' Read in the grid...')
nc=netcdf(grdname);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
hmax=max(max(nc{'h'}(:)));
result=close(nc);
%
% Create the boundary file
%
if (makebry)
  disp(' ')
  disp(' Create the boundary file...')
  create_bryfile(bryname,grdname,ROMS_title,obc,...
                 theta_s,theta_b,hc,N,...
                 bry_time,bry_cycle,'clobber');
end


%%% load file        
load interp_bnd_2001.mat

nc=netcdf(['D:\장기생태\Dynamic\07_boundary_ts\Gwangyang_bry_ykang\Gwangyang_add_uv_bry\','Gwangyang_new_Y1980.nc'],'w');

zeta_s=nc{'zeta_south'}(:);
zeta_e=nc{'zeta_east'}(:);

ubar_s=nc{'ubar_south'}(:);
ubar_e=nc{'ubar_east'}(:);

vbar_s=nc{'vbar_south'}(:);
vbar_e=nc{'vbar_east'}(:);

u_s=nc{'u_south'}(:);
u_e=nc{'u_east'}(:);

v_s=nc{'v_south'}(:);
v_e=nc{'v_east'}(:);

close(nc);


nc1=netcdf(bryname,'w');
nc1{'temp_south'}(:)=permute(d3_temp_s_t,[3 2 1]);
nc1{'temp_east'}(:)=permute(d3_temp_e_t,[3 2 1]);

nc1{'salt_south'}(:)=permute(d3_salt_s_t,[3 2 1]);
nc1{'salt_east'}(:)=permute(d3_salt_e_t,[3 2 1]);

nc1{'zeta_south'}(:)=zeta_s;
nc1{'zeta_east'}(:)=zeta_e;


nc1{'ubar_south'}(:)=ubar_s;
nc1{'ubar_east'}(:)=ubar_e;

nc1{'vbar_south'}(:)=vbar_s;
nc1{'vbar_east'}(:)=vbar_e;

nc1{'u_south'}(:)=u_s;
nc1{'u_east'}(:)=u_e;

nc1{'v_south'}(:)=v_s;
nc1{'v_east'}(:)=v_e;

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
nc1{'NO3_east'}(:)=permute(d3_no3_e_t,[3 2 1]);

nc1{'NO3_south'} = ncdouble('NO3_time','s_rho','xi_rho') ;
nc1{'NO3_south'}.long_name = ncchar('southern boundary Nitrate');
nc1{'NO3_south'}.long_name = 'southern boundary Nitrate';
nc1{'NO3_south'}.units = ncchar('mMol N m-3');
nc1{'NO3_south'}.units = 'mMol N m-3';
nc1{'NO3_south'}.coordinates = ncchar('lon_rho s_rho NO3_time');
nc1{'NO3_south'}.coordinates = 'lon_rho s_rho NO3_time';
nc1{'NO3_south'}(:)=permute(d3_no3_s_t,[3 2 1]);

nc('NH4_time') = length(bry_time);
nc{'NH4_time'} = ncdouble('NH4_time') ;
nc{'NH4_time'}.long_name = ncchar('time for Ammonium');
nc{'NH4_time'}.long_name = 'time for Ammonium';
nc{'NH4_time'}.units = ncchar('day');
nc{'NH4_time'}.units = 'day';
nc{'NH4_time'}.calendar = ncchar('360.0 days in every year');
nc{'NH4_time'}.calendar = '360.0 days in every year';
nc{'NH4_time'}.cycle_length = cycle;
nc{'NH4_time'}(:)=bry_time;

nc{'NH4_east'} = ncdouble('NH4_time','s_rho','eta_rho') ;
nc{'NH4_east'}.long_name = ncchar('eastern boundary Ammonium');
nc{'NH4_east'}.long_name = 'eastern boundary Ammonium';
nc{'NH4_east'}.units = ncchar('mMol N m-3');
nc{'NH4_east'}.units = 'mMol N m-3';
nc{'NH4_east'}.coordinates = ncchar('lon_rho s_rho NH4_time');
nc{'NH4_east'}.coordinates = 'lon_rho s_rho NH4_time';
nc{'NH4_east'}(:)=permute(d3_nh4_e_t,[3 2 1]);

nc{'NH4_south'} = ncdouble('NH4_time','s_rho','xi_rho') ;
nc{'NH4_south'}.long_name = ncchar('southern boundary Ammonium');
nc{'NH4_south'}.long_name = 'southern boundary Ammonium';
nc{'NH4_south'}.units = ncchar('mMol N m-3');
nc{'NH4_south'}.units = 'mMol N m-3';
nc{'NH4_south'}.coordinates = ncchar('lon_rho s_rho NH4_time');
nc{'NH4_south'}.coordinates = 'lon_rho s_rho NH4_time';
nc{'NH4_south'}(:)=permute(d3_nh4_s_t,[3 2 1]);

% nc('tPO4_time') = length(bry_time);
% 
% nc{'tPO4_time'} = ncdouble('PO4_time') ;
% nc{'tPO4_time'}.long_name = ncchar('time for Phosphate');
% nc{'tPO4_time'}.long_name = 'time for Phosphate';
% nc{'tPO4_time'}.units = ncchar('day');
% nc{'tPO4_time'}.units = 'day';
% nc{'tPO4_time'}.calendar = ncchar('360.0 days in every year');
% nc{'tPO4_time'}.calendar = '360.0 days in every year';
% nc{'tPO4_time'}.cycle_length = cycle;
% nc{'tPO4_time'}(:)=time;
% 
% nc{'tPO4_north'} = ncdouble('tPO4_time','s_rho','xi_rho') ;
% nc{'tPO4_north'}.long_name = ncchar('northern boundary Phosphate');
% nc{'tPO4_north'}.long_name = 'northern boundary Phosphate';
% nc{'tPO4_north'}.units = ncchar('milimole PO4 m-3');
% nc{'tPO4_north'}.units = 'milimole PO4 m-3';
% nc{'tPO4_north'}.coordinates = ncchar('lon_rho s_rho tPO4_time');
% nc{'tPO4_north'}.coordinates = 'lon_rho s_rho tPO4_time';
% nc{'tPO4_north'}(:)=po4;

nc('oxygen_time') = length(bry_time);
nc{'oxygen_time'} = ncdouble('oxygen_time') ;
nc{'oxygen_time'}.long_name = ncchar('time for oxygen');
nc{'oxygen_time'}.long_name = 'time for oxygen';
nc{'oxygen_time'}.units = ncchar('day');
nc{'oxygen_time'}.units = 'day';
nc{'oxygen_time'}.calendar = ncchar('360.0 days in every year');
nc{'oxygen_time'}.calendar = '360.0 days in every year';
nc{'oxygen_time'}.cycle_length = cycle;
nc{'oxygen_time'}(:)=bry_time;
% 
% have to oxygen tobe close
nc{'oxygen_south'} = ncdouble('oxygen_time','s_rho','xi_rho') ;
nc{'oxygen_south'}.long_name = ncchar('southern boundary oxygen');
nc{'oxygen_south'}.long_name = 'southern boundary oxygen';
nc{'oxygen_south'}.units = ncchar('mg/L');
nc{'oxygen_south'}.units = 'mg/L';
nc{'oxygen_south'}.coordinates = ncchar('lon_rho s_rho oxygen_time');
nc{'oxygen_south'}.coordinates = 'lon_rho s_rho oxygen_time';
nc{'oxygen_south'}(:)=permute(d3_nh4_s_t.*0,[3 2 1]);

nc{'oxygen_east'} = ncdouble('oxygen_time','s_rho','eta_rho') ;
nc{'oxygen_east'}.long_name = ncchar('eastern boundary oxygen');
nc{'oxygen_east'}.long_name = 'eastern boundary oxygen';
nc{'oxygen_east'}.units = ncchar('mg/L');
nc{'oxygen_east'}.units = 'mg/L';
nc{'oxygen_east'}.coordinates = ncchar('lon_rho s_rho oxygen_time');
nc{'oxygen_east'}.coordinates = 'lon_rho s_rho oxygen_time';
nc{'oxygen_east'}(:)=permute(d3_nh4_s_t.*0,[3 2 1]);

nc('chlo_time') = length(bry_time);
nc{'chlo_time'} = ncdouble('chlo_time') ;
nc{'chlo_time'}.long_name = ncchar('time for chlo');
nc{'chlo_time'}.long_name = 'time for chlo';
nc{'chlo_time'}.units = ncchar('day');
nc{'chlo_time'}.units = 'day';
nc{'chlo_time'}.calendar = ncchar('360.0 days in every year');
nc{'chlo_time'}.calendar = '360.0 days in every year';
nc{'chlo_time'}.cycle_length = cycle;
nc{'chlo_time'}(:)=bry_time;

nc{'chlo_east'} = ncdouble('chlo_time','s_rho','eta_rho') ;
nc{'chlo_east'}.long_name = ncchar('eastern boundary chlo');
nc{'chlo_east'}.long_name = 'eastern boundary chlo';
nc{'chlo_east'}.units = ncchar('mg/L');
nc{'chlo_east'}.units = 'mg/L';
nc{'chlo_east'}.coordinates = ncchar('lon_rho s_rho chlo_time');
nc{'chlo_east'}.coordinates = 'lon_rho s_rho chlo_time';
nc{'chlo_east'}(:)=permute(d3_chl_e_t,[3 2 1]);

nc{'chlo_south'} = ncdouble('chlo_time','s_rho','xi_rho') ;
nc{'chlo_south'}.long_name = ncchar('southern boundary chlo');
nc{'chlo_south'}.long_name = 'southern boundary chlo';
nc{'chlo_south'}.units = ncchar('mg/L');
nc{'chlo_south'}.units = 'mg/L';
nc{'chlo_south'}.coordinates = ncchar('lon_rho s_rho chlo_time');
nc{'chlo_south'}.coordinates = 'lon_rho s_rho chlo_time';
nc{'chlo_south'}(:)=permute(d3_chl_s_t,[3 2 1]);

% 
% nc('zoop_time') = length(bry_time);
% nc{'zoop_time'} = ncdouble('zoop_time') ;
% nc{'zoop_time'}.long_name = ncchar('time for zoop');
% nc{'zoop_time'}.long_name = 'time for zoop';
% nc{'zoop_time'}.units = ncchar('day');
% nc{'zoop_time'}.units = 'day';
% nc{'zoop_time'}.calendar = ncchar('360.0 days in every year');
% nc{'zoop_time'}.calendar = '360.0 days in every year';
% nc{'zoop_time'}.cycle_length = cycle;
% nc{'zoop_time'}(:)=time;
% 
% nc{'zoop_north'} = ncdouble('zoop_time','s_rho','xi_rho') ;
% nc{'zoop_north'}.long_name = ncchar('northern boundary zoop');
% nc{'zoop_north'}.long_name = 'northern boundary zoop';
% nc{'zoop_north'}.units = ncchar('mg/L');
% nc{'zoop_north'}.units = 'mg/L';
% nc{'zoop_north'}.coordinates = ncchar('lon_rho s_rho zoop_time');
% nc{'zoop_north'}.coordinates = 'lon_rho s_rho zoop_time';
% nc{'zoop_north'}(:)=0;
% 
% nc('phyt_time') = length(bry_time);
% nc{'phyt_time'} = ncdouble('phyt_time') ;
% nc{'phyt_time'}.long_name = ncchar('time for phyt');
% nc{'phyt_time'}.long_name = 'time for phyt';
% nc{'phyt_time'}.units = ncchar('day');
% nc{'phyt_time'}.units = 'day';
% nc{'phyt_time'}.calendar = ncchar('360.0 days in every year');
% nc{'phyt_time'}.calendar = '360.0 days in every year';
% nc{'phyt_time'}.cycle_length = cycle;
% nc{'phyt_time'}(:)=time;
% 
% nc{'phyt_north'} = ncdouble('phyt_time','s_rho','xi_rho') ;
% nc{'phyt_north'}.long_name = ncchar('northern boundary phyt');
% nc{'phyt_north'}.long_name = 'northern boundary phyt';
% nc{'phyt_north'}.units = ncchar('mg/L');
% nc{'phyt_north'}.units = 'mg/L';
% nc{'phyt_north'}.coordinates = ncchar('lon_rho s_rho phyt_time');
% nc{'phyt_north'}.coordinates = 'lon_rho s_rho phyt_time';
% nc{'phyt_north'}(:)=0;
% 
% nc('LDeN_time') = length(bry_time);
% nc{'LDeN_time'} = ncdouble('LDeN_time') ;
% nc{'LDeN_time'}.long_name = ncchar('time for LDeN');
% nc{'LDeN_time'}.long_name = 'time for LDeN';
% nc{'LDeN_time'}.units = ncchar('day');
% nc{'LDeN_time'}.units = 'day';
% nc{'LDeN_time'}.calendar = ncchar('360.0 days in every year');
% nc{'LDeN_time'}.calendar = '360.0 days in every year';
% nc{'LDeN_time'}.cycle_length = cycle;
% nc{'LDeN_time'}(:)=time;
% 
% nc{'LDeN_north'} = ncdouble('LDeN_time','s_rho','xi_rho') ;
% nc{'LDeN_north'}.long_name = ncchar('northern boundary LDeN');
% nc{'LDeN_north'}.long_name = 'northern boundary LDeN';
% nc{'LDeN_north'}.units = ncchar('mg/L');
% nc{'LDeN_north'}.units = 'mg/L';
% nc{'LDeN_north'}.coordinates = ncchar('lon_rho s_rho LDeN_time');
% nc{'LDeN_north'}.coordinates = 'lon_rho s_rho LDeN_time';
% nc{'LDeN_north'}(:)=0;
% 
% nc('SDeN_time') = length(bry_time);
% nc{'SDeN_time'} = ncdouble('SDeN_time') ;
% nc{'SDeN_time'}.long_name = ncchar('time for SDeN');
% nc{'SDeN_time'}.long_name = 'time for SDeN';
% nc{'SDeN_time'}.units = ncchar('day');
% nc{'SDeN_time'}.units = 'day';
% nc{'SDeN_time'}.calendar = ncchar('360.0 days in every year');
% nc{'SDeN_time'}.calendar = '360.0 days in every year';
% nc{'SDeN_time'}.cycle_length = cycle;
% nc{'SDeN_time'}(:)=time;
% 
% nc{'SDeN_north'} = ncdouble('SDeN_time','s_rho','xi_rho') ;
% nc{'SDeN_north'}.long_name = ncchar('northern boundary SDeN');
% nc{'SDeN_north'}.long_name = 'northern boundary SDeN';
% nc{'SDeN_north'}.units = ncchar('mg/L');
% nc{'SDeN_north'}.units = 'mg/L';
% nc{'SDeN_north'}.coordinates = ncchar('lon_rho s_rho SDeN_time');
% nc{'SDeN_north'}.coordinates = 'lon_rho s_rho SDeN_time';
% nc{'SDeN_north'}(:)=0;

close(nc1);
                     
clearvars -except i
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
