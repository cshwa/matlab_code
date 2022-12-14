%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS climatology file
%
%  Extrapole and interpole temperature and salinity from a
%  Climatology to get boundary and initial conditions for
%  ROMS (climatology and initial netcdf files) .
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
%  Data source : IRI/LDEO Climate Data Library (World Ocean Atlas 1998)
%    http://ingrid.ldgo.columbia.edu/
%    http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NODC/.WOA98/
%
%  Pierrick Penven, IRD, 2002.
%
%  Version of 10-Oct-2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
romstools_param
%
%  Title 
%
title='Climatology';
%
%  Switches for selecting what to process (1=ON)
%
makeclim=1; %1: process boundary data
makeoa=1;   %1: process oa data
makeini=1;  %1: process initial data
%
%  Grid file name - Climatology file name
%  Initial file name - OA file name
%
grdname='roms_grd.nc';
clmname='roms_clm.nc';
ininame='roms_ini.nc';
oaname ='roms_oa.nc';
%
%  Day of initialisation
%
tini=15;  
%
% Set times and cycles: monthly climatology for all data
%
time=[15:30:345];    % time 
cycle=360;           % cycle
%
%  Data climatology file name:
%
datafile='../../OCCAM_AGULHAS/Matlab/occam_agulhas_ext_clm.nc';
%
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%
disp(' ')
disp([' Making the clim: ',clmname])
disp(' ')
disp([' Title: ',title])
%
% Read in the grid
%
disp(' ')
disp(' Read in the grid...')
nc=netcdf(grdname);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
hmax=max(max(nc{'h'}(:)));
result=close(nc);
%
% Create the climatology file
%
if (makeclim)
  disp(' ')
  disp(' Create the climatology file...')
  create_climfile(clmname,grdname,title,...
                  theta_s,theta_b,hc,N,...
                  time,cycle,'clobber');
end
%
% Create the OA file
%
if (makeoa)
  disp(' ')
  disp(' Create the OA file...')
  nc=netcdf(datafile);
  Z=nc{'Z'}(:);
  close(nc)
  create_oafile(oaname,grdname,title,Z,...
                time,cycle,'clobber');
%
% Horizontal extrapolations 
%
  disp(' ')
  disp(' Horizontal interpolations')
  disp(' ')
  disp(' Temperature...')
  int_data_occam(oaname,datafile,'temp','lon_rho','lat_rho','temp','lonT','latT',1);
  disp(' ')
  disp(' Salinity...')
  int_data_occam(oaname,datafile,'salt','lon_rho','lat_rho','salt','lonT','latT',1);
  disp(' ')
  disp(' U...')
  int_data_occam(oaname,datafile,'u','lon_u','lat_u','u','lonU','latU',1);
  disp(' ')
  disp(' V...')
  int_data_occam(oaname,datafile,'v','lon_v','lat_v','v','lonU','latU',1);
  disp(' ')
  disp(' UBAR...')
  int_data_occam(oaname,datafile,'ubar','lon_u','lat_u','ubar','lonU','latU',0);
  disp(' ')
  disp(' VBAR...')
  int_data_occam(oaname,datafile,'vbar','lon_v','lat_v','vbar','lonU','latU',0);
  disp(' ')
  disp(' ZETA...')
  int_data_occam(oaname,datafile,'zeta','lon_rho','lat_rho','zeta','lonT','latT',0);
  int_data_occam(oaname,datafile,'ssh','lon_rho','lat_rho','zeta','lonT','latT',0);
end
%
% Vertical interpolations 
%
if (makeclim)
  disp(' ')
  disp(' Vertical interpolations')
  disp(' ')
  disp(' Temperature...')
  vinterp_clm(clmname,grdname,oaname,'temp','tclm_time','Z',0,'r');
  disp(' ')
  disp(' Salinity...')
  vinterp_clm(clmname,grdname,oaname,'salt','sclm_time','Z',0,'r');
  disp(' ')
  disp(' U...')
  vinterp_clm(clmname,grdname,oaname,'u','uclm_time','Z',0,'u');
  disp(' ')
  disp(' V...')
  vinterp_clm(clmname,grdname,oaname,'v','vclm_time','Z',0,'v');
  disp(' ')
  disp(' ZETA...')
  rem_avg_zeta(clmname,grdname,oaname)
  disp(' ')
  disp(' UBAR & VBAR...')
  barotropic_currents(clmname,grdname,obc) 
end
%
% Initial file
%
if (makeini)
  disp(' ')
  disp(' Initial')
  create_inifile(ininame,grdname,title,...
                 theta_s,theta_b,hc,N,...
                 tini,'clobber');
  disp(' ')
  disp(' Temperature...')
  vinterp_clm(ininame,grdname,oaname,'temp','tclm_time','Z',tini,'r',1);
  disp(' ')
  disp(' Salinity...')
  vinterp_clm(ininame,grdname,oaname,'salt','sclm_time','Z',tini,'r',1);
end		 
%
% Make a few plots
%
disp(' ')
disp(' Make a few plots...')
test_clim(clmname,grdname,'temp',1)
figure
test_clim(clmname,grdname,'salt',1)
figure
test_clim(clmname,grdname,'temp',6)
figure
test_clim(clmname,grdname,'salt',6)
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
