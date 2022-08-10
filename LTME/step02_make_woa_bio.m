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

%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
addpath(genpath(['D:\romsplot']));
for year=2001:2012;
clearvars -except year
close all
clc
warning off
ROMS_title='ECSandYS fennel';
bryname=['.\NWP_1_10_test6\data\roms_bndy_auto_NWP_',num2str(year),'_fennel_woa13.nc'];
% bryname=['.\NWP_1_10_test6\data\roms_bndy_auto_NWP_1_10_test6_clim_fennel_woa13.nc'];
grdname=['g:\auto_fennel\grid\roms_grd_auto_rdrg2_new8_smooth.nc'];
Zbryname=['.\NWP_1_10_test6\data\roms_bndy_auto_NWP_',num2str(year),'_Zfennel_woa13.nc'];
% Zbryname=['.\NWP_1_10_test6\data\roms_bndy_auto_NWP_1_10_test6_clim_Zfennel_woa13.nc'];
makebry=1;
makeZbry=1;
makepisces=0;
makeplot=1;
insitu2pot=1;   %1: transform in-situ temperature to potential temperature

Vtransform=2;
Vstretching=4;
theta_s = 5;
theta_b = 0.4;
hc      =500;
N=40;

%
% Objective analysis decorrelation scale [m]
% (if Roa=0: simple extrapolation method; crude but much less costly)
%
% Roa=300e3;
Roa=0;
%

obc = [1 1 0 0]; % open boundaries (1=open , [S E N W])

%
%  Data climatologies file names:
%
%    temp_month_data : monthly temperature climatology
%    temp_ann_data   : annual temperature climatology
%    salt_month_data : monthly salinity climatology
%    salt_ann_data   : annual salinity climatology
%
woa_dir=['d:\OneDrive - SNU\woa2013v2\'];
NO3_ann_data=[woa_dir,'WOA13_nitrate_annual_1deg.nc'];
NO3_month_data=[woa_dir,'WOA13_nitrate_monthly_1deg.nc'];
NH4_ann_data=[woa_dir,'WOA13_nitrate_annual_1deg.nc'];
NH4_month_data=[woa_dir,'WOA13_nitrate_monthly_1deg.nc'];
oxy_ann_data=[woa_dir,'WOA13_oxygen_annual_1deg.nc'];
oxy_month_data=[woa_dir,'WOA13_oxygen_monthly_1deg.nc'];
chl_ann_data=[woa_dir,'chla_seawifs_annual_nwp_climate_ver_woa13.nc'];
chl_month_data=[woa_dir,'chla_seawifs_month_nwp_climate_ver_woa13.nc'];
tPO4_ann_data=[woa_dir,'WOA13_phosphate_annual_1deg.nc'];
tPO4_month_data=[woa_dir,'WOA13_phosphate_monthly_1deg.nc'];

bry_time=[15:30:365];
% chl_bry_time=[45:90:365];
bry_cycle=yeardays(year);

%
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%
disp(' ')
disp([' adding the ecosystem file: ',bryname])

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
  create_bryfile_fennel(bryname,grdname,ROMS_title,obc,...
                 theta_s,theta_b,hc,N,...
                 bry_time,bry_cycle,'clobber',Vtransform,Vstretching);
end
%
% Create the boundary file in Z-coordinates
%
if (makeZbry)
  disp(' ')
  disp(' Create the boundary Z-file...')
%
% get Z
%
  nc=netcdf(NO3_ann_data);
  Z=nc{'depth'}(:);
  kmax=max(find(Z<hmax));
  Z=Z(1:kmax);
  close(nc);
  create_bry_Z_fennel(Zbryname,grdname,ROMS_title,obc,...
                Z,bry_time,bry_cycle,'clobber');
  disp(' ')
  disp(' Horizontal extrapolations')
%
% Loop on the lateral boundaries 
%
  for obcndx=1:4
    if obc(obcndx)==1
      if obcndx==1
        disp(' Processing southern boundary...')
	suffix='_south';
      elseif obcndx==2
        disp(' Processing eastern boundary...')
	suffix='_east';
      elseif obcndx==3
        disp(' Processing northern boundary...')
	suffix='_north';
      elseif obcndx==4
        disp(' Processing western boundary...')
	suffix='_west';
      end
      disp('  NO3...')
      bry_interp_2009(Zbryname,lon,lat,NO3_month_data,NO3_ann_data,...
               'n_an',['NO3',suffix],obcndx,Roa);
      disp('  NH4...')
      bry_interp_2009(Zbryname,lon,lat,NH4_month_data,NH4_ann_data,...
               'n_an',['NH4',suffix],obcndx,Roa);     
      disp('  Oxygen...')
      bry_interp_2009(Zbryname,lon,lat,oxy_month_data,oxy_ann_data,...
               'o_an',['oxygen',suffix],obcndx,Roa);
      disp('  Chlorophyll a...')
      bry_interp_2009(Zbryname,lon,lat,chl_month_data,chl_ann_data,...
               'chlorophyll',['chlo',suffix],obcndx,Roa);            
      disp('  PO4...')
      bry_interp_2009(Zbryname,lon,lat,tPO4_month_data,tPO4_ann_data,...
               'p_an',['tPO4',suffix],obcndx,Roa);             
           
    end
  end
end
%
% Vertical interpolations 
%
if (makebry)
  disp(' ')
  disp(' Vertical interpolations')
%
% Loop on the lateral boundaries 
%
  for obcndx=1:4
    if obc(obcndx)==1
      if obcndx==1
        disp(' Processing southern boundary...')
	suffix='_south';
      elseif obcndx==2
        disp(' Processing eastern boundary...')
	suffix='_east';
      elseif obcndx==3
        disp(' Processing northern boundary...')
	suffix='_north';
      elseif obcndx==4
        disp(' Processing western boundary...')
	suffix='_west';
      end
      disp('  NO3...')
      vinterp_bry(bryname,grdname,Zbryname,['NO3',suffix],obcndx,Vtransform,Vstretching);
      disp('  NH4...')
      vinterp_bry(bryname,grdname,Zbryname,['NH4',suffix],obcndx,Vtransform,Vstretching);
      disp('  chlorophyll...')
      vinterp_bry(bryname,grdname,Zbryname,['chlo',suffix],obcndx,Vtransform,Vstretching);
      disp('  oxygen...')
      vinterp_bry(bryname,grdname,Zbryname,['oxygen',suffix],obcndx,Vtransform,Vstretching);
      disp('  tPO4...')
      vinterp_bry(bryname,grdname,Zbryname,['tPO4',suffix],obcndx,Vtransform,Vstretching);
      if (insitu2pot)
        disp(' ')
        disp('  Compute potential temperature from in-situ...')
        getpot_bry(bryname,grdname,obcndx,Vtransform,Vstretching)
      end
%
% Geostrophy
%
%       disp(' ')
%       disp('  Compute geostrophic currents')
%       geost_currents_bry(bryname,grdname,Zbryname,frcname,zref,obcndx)
    end
  end
%
% Remove avg SSH
%
  rmavgssh(bryname,grdname,obc)

%
%% Compute bry for pisces variables
%
 if makepisces
   disp('====================================== ')
   disp('Compute boundary for Pisces tracer')
   make_bry_pisces
 end
end
%
% Make a few plots
%
if makeplot==1
  disp(' ')
  disp(' Make a few plots...')
  figure
  test_bry(bryname,grdname,Vtransform,Vstretching,'NO3',12,obc)
  figure
  test_bry(bryname,grdname,Vtransform,Vstretching,'NH4',12,obc)  
  figure
  test_bry(bryname,grdname,Vtransform,Vstretching,'chlo',12,obc)
  figure
  test_bry(bryname,grdname,Vtransform,Vstretching,'oxygen',12,obc)    
  figure
  test_bry(bryname,grdname,Vtransform,Vstretching,'tPO4',12,obc)   
end
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end