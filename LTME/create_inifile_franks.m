function create_inifile(inifile,gridfile,title,...
                         theta_s,theta_b,hc,N,time,clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function nc=create_inifile(inifile,gridfile,theta_s,...
%                  theta_b,hc,N,ttime,stime,utime,... 
%                  cycle,clobber)
%
%   This function create the header of a Netcdf climatology 
%   file.
%
%   Input: 
% 
%   inifile      Netcdf initial file name (character string).
%   gridfile     Netcdf grid file name (character string).
%   theta_s      S-coordinate surface control parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   hc           Width (m) of surface or bottom boundary layer
%                where higher vertical resolution is required 
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)  
%   time         Initial time.(Real) 
%   clobber      Switch to allow or not writing over an existing
%                file.(character string) 
%
%   Output
%
%   nc       Output netcdf object.
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
%  Copyright (c) 2001-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp([' Creating the file : ',inifile])
%
%  Read the grid file
%
nc=netcdf(gridfile);
h=nc{'h'}(:);  
mask=nc{'mask_rho'}(:);
close(nc);
hmin=min(min(h(mask==1)));
if hc > hmin
  error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
end
[Mp,Lp]=size(h);
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the initial file
%
type = 'INITIAL file' ; 
history = 'ROMS' ;
nc = netcdf(inifile,clobber);
result = redef(nc);
%
%  Create dimensions
%
nc('xi_u') = L;
nc('xi_v') = Lp;
nc('xi_rho') = Lp;
nc('eta_u') = Mp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
nc('s_rho') = N;
nc('s_w') = Np;
nc('tracer') = 2;
nc('time') = 0;
nc('one') = 1;
%
%  Create variables
%
nc{'tstart'} = ncdouble('one') ;
nc{'tend'} = ncdouble('one') ;
nc{'theta_s'} = ncdouble('one') ;
nc{'theta_b'} = ncdouble('one') ;
nc{'Tcline'} = ncdouble('one') ;
nc{'hc'} = ncdouble('one') ;
nc{'sc_r'} = ncdouble('s_rho') ;
nc{'Cs_r'} = ncdouble('s_rho') ;
nc{'ocean_time'} = ncdouble('time') ;
nc{'scrum_time'} = ncdouble('time') ;
nc{'u'} = ncdouble('time','s_rho','eta_u','xi_u') ;
nc{'v'} = ncdouble('time','s_rho','eta_v','xi_v') ;
nc{'ubar'} = ncdouble('time','eta_u','xi_u') ;
nc{'vbar'} = ncdouble('time','eta_v','xi_v') ;
nc{'zeta'} = ncdouble('time','eta_rho','xi_rho') ;
nc{'temp'} = ncdouble('time','s_rho','eta_rho','xi_rho') ;
nc{'salt'} = ncdouble('time','s_rho','eta_rho','xi_rho') ;
nc{'NO3'} = ncdouble('time','s_rho','eta_rho','xi_rho') ;
nc{'phytoplankton'} = ncdouble('time','s_rho','eta_rho','xi_rho') ;
nc{'zooplankton'} = ncdouble('time','s_rho','eta_rho','xi_rho') ;
nc{'detritus'} = ncdouble('time','s_rho','eta_rho','xi_rho') ;
%
%  Create attributes
%
nc{'tstart'}.long_name = ncchar('start processing day');
nc{'tstart'}.long_name = 'start processing day';
nc{'tstart'}.units = ncchar('day');
nc{'tstart'}.units = 'day';
%
nc{'tend'}.long_name = ncchar('end processing day');
nc{'tend'}.long_name = 'end processing day';
nc{'tend'}.units = ncchar('day');
nc{'tend'}.units = 'day';
%
nc{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
nc{'theta_s'}.long_name = 'S-coordinate surface control parameter';
nc{'theta_s'}.units = ncchar('nondimensional');
nc{'theta_s'}.units = 'nondimensional';
%
nc{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
nc{'theta_b'}.long_name = 'S-coordinate bottom control parameter';
nc{'theta_b'}.units = ncchar('nondimensional');
nc{'theta_b'}.units = 'nondimensional';
%
nc{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
nc{'Tcline'}.long_name = 'S-coordinate surface/bottom layer width';
nc{'Tcline'}.units = ncchar('meter');
nc{'Tcline'}.units = 'meter';
%
nc{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
nc{'hc'}.long_name = 'S-coordinate parameter, critical depth';
nc{'hc'}.units = ncchar('meter');
nc{'hc'}.units = 'meter';
%
nc{'sc_r'}.long_name = ncchar('S-coordinate at RHO-points');
nc{'sc_r'}.long_name = 'S-coordinate at RHO-points';
nc{'sc_r'}.units = ncchar('nondimensional');
nc{'sc_r'}.units = 'nondimensional';
nc{'sc_r'}.valid_min = -1;
nc{'sc_r'}.valid_max = 0;
%
nc{'Cs_r'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nc{'Cs_r'}.long_name = 'S-coordinate stretching curves at RHO-points';
nc{'Cs_r'}.units = ncchar('nondimensional');
nc{'Cs_r'}.units = 'nondimensional';
nc{'Cs_r'}.valid_min = -1;
nc{'Cs_r'}.valid_max = 0;
%
nc{'ocean_time'}.long_name = ncchar('time since initialization');
nc{'ocean_time'}.long_name = 'time since initialization';
nc{'ocean_time'}.units = ncchar('second');
nc{'ocean_time'}.units = 'second';
%
nc{'scrum_time'}.long_name = ncchar('time since initialization');
nc{'scrum_time'}.long_name = 'time since initialization';
nc{'scrum_time'}.units = ncchar('second');
nc{'scrum_time'}.units = 'second';
%
nc{'u'}.long_name = ncchar('u-momentum component');
nc{'u'}.long_name = 'u-momentum component';
nc{'u'}.units = ncchar('meter second-1');
nc{'u'}.units = 'meter second-1';
%
nc{'v'}.long_name = ncchar('v-momentum component');
nc{'v'}.long_name = 'v-momentum component';
nc{'v'}.units = ncchar('meter second-1');
nc{'v'}.units = 'meter second-1';
%
nc{'ubar'}.long_name = ncchar('vertically integrated u-momentum component');
nc{'ubar'}.long_name = 'vertically integrated u-momentum component';
nc{'ubar'}.units = ncchar('meter second-1');
nc{'ubar'}.units = 'meter second-1';
%
nc{'vbar'}.long_name = ncchar('vertically integrated v-momentum component');
nc{'vbar'}.long_name = 'vertically integrated v-momentum component';
nc{'vbar'}.units = ncchar('meter second-1');
nc{'vbar'}.units = 'meter second-1';
%
nc{'zeta'}.long_name = ncchar('free-surface');
nc{'zeta'}.long_name = 'free-surface';
nc{'zeta'}.units = ncchar('meter');
nc{'zeta'}.units = 'meter';
%
nc{'temp'}.long_name = ncchar('potential temperature');
nc{'temp'}.long_name = 'potential temperature';
nc{'temp'}.units = ncchar('Celsius');
nc{'temp'}.units = 'Celsius';
%
nc{'salt'}.long_name = ncchar('salinity');
nc{'salt'}.long_name = 'salinity';
nc{'salt'}.units = ncchar('PSU');
nc{'salt'}.units = 'PSU';

nc{'salt'}.long_name = ncchar('salinity');
nc{'salt'}.long_name = 'salinity';
nc{'salt'}.units = ncchar('PSU');
nc{'salt'}.units = 'PSU';

nc{'NO3'}.long_name = ncchar('nitrate');
nc{'NO3'}.long_name = 'nitrate';
nc{'NO3'}.units = ncchar('millie mole N meter-3');
nc{'NO3'}.units = 'millie mole N meter-3';

nc{'phytoplankton'}.long_name = ncchar('phytoplankton');
nc{'phytoplankton'}.long_name = 'phytoplankton';
nc{'phytoplankton'}.units = ncchar('millie mole N meter-3');
nc{'phytoplankton'}.units = 'millie mole N meter-3';

nc{'zooplankton'}.long_name = ncchar('zooplankton');
nc{'zooplankton'}.long_name = 'zooplankton';
nc{'zooplankton'}.units = ncchar('millie mole N meter-3');
nc{'zooplankton'}.units = 'millie mole N meter-3';

nc{'detritus'}.long_name = ncchar('detritus');
nc{'detritus'}.long_name = 'detritus';
nc{'detritus'}.units = ncchar('millie mole N meter-3');
nc{'detritus'}.units = 'millie mole N meter-3';
%
% Create global attributes
%
nc.title = ncchar(title);
nc.title = title;
nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(inifile);
nc.clim_file = inifile;
nc.grd_file = ncchar(gridfile);
nc.grd_file = gridfile;
nc.type = ncchar(type);
nc.type = type;
nc.history = ncchar(history);
nc.history = history;
%
% Leave define mode
%
result = endef(nc);
%
% Compute S coordinates
%
cff1=1./sinh(theta_s);
cff2=0.5/tanh(0.5*theta_s);
sc=((1:N)-N-0.5)/N;
Cs=(1.-theta_b)*cff1*sinh(theta_s*sc)...
    +theta_b*(cff2*tanh(theta_s*(sc+0.5))-0.5);
%
% Write variables
%
nc{'tstart'}(:) =  time; 
nc{'tend'}(:) =  time; 
nc{'theta_s'}(:) =  theta_s; 
nc{'theta_b'}(:) =  theta_b; 
nc{'Tcline'}(:) =  hc; 
nc{'hc'}(:) =  hc; 
nc{'sc_r'}(:) =  sc; 
nc{'Cs_r'}(:) =  Cs; 
nc{'scrum_time'}(1) =  time*24*3600; 
nc{'ocean_time'}(1) =  time*24*3600; 
nc{'u'}(:) =  0; 
nc{'v'}(:) =  0; 
nc{'zeta'}(:) =  0; 
nc{'ubar'}(:) =  0; 
nc{'vbar'}(:) =  0; 
nc{'temp'}(:) =  0; 
nc{'salt'}(:) =  0; 
nc{'NO3'}(:) =  0; 
nc{'phytoplankton'}(:) =  0; 
nc{'zooplankton'}(:) =  0; 
nc{'detritus'}(:) =  0; 
%
% Synchronize on disk
%
close(nc);
return


