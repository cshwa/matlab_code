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
for i = 1980:1989
    yeart=i;
bryname=['GY_OBC_',num2str(yeart),'.nc'];
grdname=['grid_sumjin_v1970_fix_3m.nc'];
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
nc=netcdf(['./'1980_ini_his_8785.nc','w');
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
nc1=netcdf('gy_1980_tsbio_ini.nc','w');
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



clearvars -except i
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
