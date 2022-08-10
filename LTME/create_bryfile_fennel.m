function create_bryfile(bryname,grdname,title,obc,...
                        theta_s,theta_b,hc,N,...
                        time,cycle,clobber);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function create_bryfile(bryname,grdname,title,obc...
%                          theta_s,theta_b,hc,N,...
%                          time,cycle,clobber);
%
%   This function create the header of a Netcdf climatology 
%   file.
%
%   Input:
%
%   bryname      Netcdf climatology file name (character string).
%   grdname      Netcdf grid file name (character string).
%   obc          open boundaries flag (1=open , [S E N W]).
%   theta_s      S-coordinate surface control parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   hc           Width (m) of surface or bottom boundary layer 
%                where higher vertical resolution is required 
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)
%   time         time.(vector)
%   cycle        Length (days) for cycling the climatology.(Real)
%   clobber      Switch to allow or not writing over an existing
%                file.(character string)
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
disp([' Creating the file : ',bryname])
disp(' ')
%
%  Read the grid file and check the topography
%
nc = netcdf(grdname, 'nowrite');
h=nc{'h'}(:);
maskr=nc{'mask_rho'}(:);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
status=close(nc);
hmin=min(min(h(maskr==1)));
if hc > hmin
  error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
end
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the boundary file
%
type = 'BOUNDARY file' ; 
history = 'ROMS' ;
nc = netcdf(bryname,clobber);
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
nc('bry_time') = length(time);
nc('tclm_time') = length(time);
nc('temp_time') = length(time);
nc('sclm_time') = length(time);
nc('salt_time') = length(time);
nc('uclm_time') = length(time);
nc('vclm_time') = length(time);
nc('v2d_time')  = length(time);
nc('v3d_time')  = length(time);
nc('ssh_time')  = length(time);
nc('zeta_time') = length(time);
nc('NO3_time') = length(time);
nc('NH4_time') = length(time);
nc('chlo_time') = length(time);
nc('oxygen_time') = length(time);
nc('tPO4_time') = length(time);
nc('one') = 1;
%
%  Create variables and attributes
%
nc{'spherical'} = ncchar('one') ;
nc{'spherical'}.long_name = ncchar('grid type logical switch');
nc{'spherical'}.long_name = 'grid type logical switch';
nc{'spherical'}.flag_values = ncchar('T, F');
nc{'spherical'}.flag_values = 'T, F';
nc{'spherical'}.flag_meanings = ncchar('spherical Cartesian');
nc{'spherical'}.flag_meanings = 'spherical Cartesian';
%
nc{'Vtransform'} = ncint('one') ;
nc{'Vtransform'}.long_name = ncchar('vertical terrain-following transformation equation');
nc{'Vtransform'}.long_name = 'vertical terrain-following transformation equation';
%
nc{'Vstretching'} = ncint('one') ;
nc{'Vstretching'}.long_name = ncchar('vertical terrain-following stretching function');
nc{'Vstretching'}.long_name = 'vertical terrain-following stretching function';
%
nc{'tstart'} = ncdouble('one') ;
nc{'tstart'}.long_name = ncchar('start processing day');
nc{'tstart'}.long_name = 'start processing day';
nc{'tstart'}.units = ncchar('day');
nc{'tstart'}.units = 'day';
%
nc{'tend'} = ncdouble('one') ;
nc{'tend'}.long_name = ncchar('end processing day');
nc{'tend'}.long_name = 'end processing day';
nc{'tend'}.units = ncchar('day');
nc{'tend'}.units = 'day';
%
nc{'theta_s'} = ncdouble('one') ;
nc{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
nc{'theta_s'}.long_name = 'S-coordinate surface control parameter';
nc{'theta_s'}.units = ncchar('nondimensional');
nc{'theta_s'}.units = 'nondimensional';
%
nc{'theta_b'} = ncdouble('one') ;
nc{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
nc{'theta_b'}.long_name = 'S-coordinate bottom control parameter';
nc{'theta_b'}.units = ncchar('nondimensional');
nc{'theta_b'}.units = 'nondimensional';
%
nc{'Tcline'} = ncdouble('one') ;
nc{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
nc{'Tcline'}.long_name = 'S-coordinate surface/bottom layer width';
nc{'Tcline'}.units = ncchar('meter');
nc{'Tcline'}.units = 'meter';
%
nc{'hc'} = ncdouble('one') ;
nc{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
nc{'hc'}.long_name = 'S-coordinate parameter, critical depth';
nc{'hc'}.units = ncchar('meter');
nc{'hc'}.units = 'meter';
%
nc{'s_rho'} = ncdouble('s_rho') ;
nc{'s_rho'}.long_name = ncchar('S-coordinate at RHO-points');
nc{'s_rho'}.long_name = 'S-coordinate at RHO-points';
nc{'s_rho'}.valid_min = -1.;
nc{'s_rho'}.valid_max = 0.;
nc{'s_rho'}.positive = ncchar('up');
nc{'s_rho'}.positive = 'up';
nc{'s_rho'}.standard_name = ncchar('ocena_s_coordinate_g2');
nc{'s_rho'}.standard_name = 'ocena_s_coordinate_g2';
nc{'s_rho'}.formula_terms = ncchar('s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc');
nc{'s_rho'}.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc';
%
nc{'s_w'} = ncdouble('s_w') ;
nc{'s_w'}.long_name = ncchar('S-coordinate at W-points');
nc{'s_w'}.long_name = 'S-coordinate at W-points';
nc{'s_w'}.valid_min = -1. ;
nc{'s_w'}.valid_max = 0. ;
nc{'s_w'}.positive = ncchar('up');
nc{'s_w'}.positive = 'up';
nc{'s_w'}.standard_name = ncchar('ocena_s_coordinate_g2');
nc{'s_w'}.standard_name = 'ocena_s_coordinate_g2';
nc{'s_w'}.formula_terms = ncchar('s: s_w C: Cs_w eta: zeta depth: h depth_c: hc');
nc{'s_w'}.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc';
%
nc{'Cs_r'} = ncdouble('s_rho') ;
nc{'Cs_r'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nc{'Cs_r'}.long_name = 'S-coordinate stretching curves at RHO-points';
nc{'Cs_r'}.units = ncchar('nondimensional');
nc{'Cs_r'}.units = 'nondimensional';
nc{'Cs_r'}.valid_min = -1;
nc{'Cs_r'}.valid_max = 0;
%
nc{'Cs_w'} = ncdouble('s_w') ;
nc{'Cs_w'}.long_name = ncchar('S-coordinate stretching curves at W-points');
nc{'Cs_w'}.long_name = 'S-coordinate stretching curves at W-points';
nc{'Cs_w'}.units = ncchar('nondimensional');
nc{'Cs_w'}.units = 'nondimensional';
nc{'Cs_w'}.valid_min = -1;
nc{'Cs_w'}.valid_max = 0;
%
nc{'bry_time'} = ncdouble('bry_time') ;
nc{'bry_time'}.long_name = ncchar('time for boundary climatology');
nc{'bry_time'}.long_name = 'time for boundary climatology';
nc{'bry_time'}.units = ncchar('day');
nc{'bry_time'}.units = 'day';
nc{'bry_time'}.calendar = ncchar('360.0 days in every year');
nc{'bry_time'}.calendar = '360.0 days in every year';
nc{'bry_time'}.cycle_length = cycle;
%
nc{'tclm_time'} = ncdouble('tclm_time') ;
nc{'tclm_time'}.long_name = ncchar('time for temperature climatology');
nc{'tclm_time'}.long_name = 'time for temperature climatology';
nc{'tclm_time'}.units = ncchar('day');
nc{'tclm_time'}.units = 'day';
nc{'tclm_time'}.calendar = ncchar('360.0 days in every year');
nc{'tclm_time'}.calendar = '360.0 days in every year';
nc{'tclm_time'}.cycle_length = cycle;
%
nc{'temp_time'} = ncdouble('temp_time') ;
nc{'temp_time'}.long_name = ncchar('time for temperature climatology');
nc{'temp_time'}.long_name = 'time for temperature climatology';
nc{'temp_time'}.units = ncchar('day');
nc{'temp_time'}.units = 'day';
nc{'temp_time'}.calendar = ncchar('360.0 days in every year');
nc{'temp_time'}.calendar = '360.0 days in every year';
nc{'temp_time'}.cycle_length = cycle;
%
nc{'NO3_time'} = ncdouble('NO3_time') ;
nc{'NO3_time'}.long_name = ncchar('time for nitrate climatology');
nc{'NO3_time'}.long_name = 'time for nitrate climatology';
nc{'NO3_time'}.units = ncchar('day');
nc{'NO3_time'}.units = 'day';
nc{'NO3_time'}.calendar = ncchar('360.0 days in every year');
nc{'NO3_time'}.calendar = '360.0 days in every year';
nc{'NO3_time'}.cycle_length = cycle;
%
nc{'NH4_time'} = ncdouble('NH4_time') ;
nc{'NH4_time'}.long_name = ncchar('time for ammonium climatology');
nc{'NH4_time'}.long_name = 'time for ammonium climatology';
nc{'NH4_time'}.units = ncchar('day');
nc{'NH4_time'}.units = 'day';
nc{'NH4_time'}.calendar = ncchar('360.0 days in every year');
nc{'NH4_time'}.calendar = '360.0 days in every year';
nc{'NH4_time'}.cycle_length = cycle;
%
nc{'oxygen_time'} = ncdouble('oxygen_time') ;
nc{'oxygen_time'}.long_name = ncchar('time for oxygen climatology');
nc{'oxygen_time'}.long_name = 'time for oxygen climatology';
nc{'oxygen_time'}.units = ncchar('day');
nc{'oxygen_time'}.units = 'day';
nc{'oxygen_time'}.calendar = ncchar('360.0 days in every year');
nc{'oxygen_time'}.calendar = '360.0 days in every year';
nc{'oxygen_time'}.cycle_length = cycle;
%
nc{'chlo_time'} = ncdouble('chlo_time') ;
nc{'chlo_time'}.long_name = ncchar('time for chlorophyll climatology');
nc{'chlo_time'}.long_name = 'time for chlorophyll climatology';
nc{'chlo_time'}.units = ncchar('day');
nc{'chlo_time'}.units = 'day';
nc{'chlo_time'}.calendar = ncchar('360.0 days in every year');
nc{'chlo_time'}.calendar = '360.0 days in every year';
nc{'chlo_time'}.cycle_length = cycle;
%
nc{'tPO4_time'} = ncdouble('tPO4_time') ;
nc{'tPO4_time'}.long_name = ncchar('time for phosphate climatology');
nc{'tPO4_time'}.long_name = 'time for phosphate climatology';
nc{'tPO4_time'}.units = ncchar('day');
nc{'tPO4_time'}.units = 'day';
nc{'tPO4_time'}.calendar = ncchar('360.0 days in every year');
nc{'tPO4_time'}.calendar = '360.0 days in every year';
nc{'tPO4_time'}.cycle_length = cycle;
%
nc{'sclm_time'} = ncdouble('sclm_time') ;
nc{'sclm_time'}.long_name = ncchar('time for salinity climatology');
nc{'sclm_time'}.long_name = 'time for salinity climatology';
nc{'sclm_time'}.units = ncchar('day');
nc{'sclm_time'}.units = 'day';
nc{'sclm_time'}.calendar = ncchar('360.0 days in every year');
nc{'sclm_time'}.calendar = '360.0 days in every year';
nc{'sclm_time'}.cycle_length = cycle;
%
nc{'salt_time'} = ncdouble('salt_time') ;
nc{'salt_time'}.long_name = ncchar('time for salinity climatology');
nc{'salt_time'}.long_name = 'time for salinity climatology';
nc{'salt_time'}.units = ncchar('day');
nc{'salt_time'}.units = 'day';
nc{'salt_time'}.calendar = ncchar('360.0 days in every year');
nc{'salt_time'}.calendar = '360.0 days in every year';
nc{'salt_time'}.cycle_length = cycle;
%
nc{'uclm_time'} = ncdouble('uclm_time') ;
nc{'uclm_time'}.long_name = ncchar('time climatological u');
nc{'uclm_time'}.long_name = 'time climatological u';
nc{'uclm_time'}.units = ncchar('day');
nc{'uclm_time'}.units = 'day';
nc{'uclm_time'}.calendar = ncchar('360.0 days in every year');
nc{'uclm_time'}.calendar = '360.0 days in every year';
nc{'uclm_time'}.cycle_length = cycle;
%
nc{'vclm_time'} = ncdouble('vclm_time') ;
nc{'vclm_time'}.long_name = ncchar('time climatological v');
nc{'vclm_time'}.long_name = 'time climatological v';
nc{'vclm_time'}.units = ncchar('day');
nc{'vclm_time'}.units = 'day';
nc{'vclm_time'}.calendar = ncchar('360.0 days in every year');
nc{'vclm_time'}.calendar = '360.0 days in every year';
nc{'vclm_time'}.cycle_length = cycle;
%
nc{'v2d_time'} = ncdouble('v2d_time') ;
nc{'v2d_time'}.long_name = ncchar('time for 2D velocity climatology');
nc{'v2d_time'}.long_name = 'time for 2D velocity climatology';
nc{'v2d_time'}.units = ncchar('day');
nc{'v2d_time'}.units = 'day';
nc{'v2d_time'}.calendar = ncchar('360.0 days in every year');
nc{'v2d_time'}.calendar = '360.0 days in every year';
nc{'v2d_time'}.cycle_length = cycle;
%
nc{'v3d_time'} = ncdouble('v3d_time') ;
nc{'v3d_time'}.long_name = ncchar('time for 3D velocity climatology');
nc{'v3d_time'}.long_name = 'time for 3D velocity climatology';
nc{'v3d_time'}.units = ncchar('day');
nc{'v3d_time'}.units = 'day';
nc{'v3d_time'}.calendar = ncchar('360.0 days in every year');
nc{'v3d_time'}.calendar = '360.0 days in every year';
nc{'v3d_time'}.cycle_length = cycle;
%
nc{'ssh_time'} = ncdouble('ssh_time') ;
nc{'ssh_time'}.long_name = ncchar('time for sea surface height');
nc{'ssh_time'}.long_name = 'time for sea surface height';
nc{'ssh_time'}.units = ncchar('day');
nc{'ssh_time'}.units = 'day';
nc{'ssh_time'}.calendar = ncchar('360.0 days in every year');
nc{'ssh_time'}.calendar = '360.0 days in every year';
nc{'ssh_time'}.cycle_length = cycle;
%
nc{'zeta_time'} = ncdouble('zeta_time') ;
nc{'zeta_time'}.long_name = ncchar('time for sea surface height');
nc{'zeta_time'}.long_name = 'time for sea surface height';
nc{'zeta_time'}.units = ncchar('day');
nc{'zeta_time'}.units = 'day';
nc{'zeta_time'}.calendar = ncchar('360.0 days in every year');
nc{'zeta_time'}.calendar = '360.0 days in every year';
nc{'zeta_time'}.cycle_length = cycle;
%
if obc(1)==1
%
%   Southern boundary
%
  nc{'temp_south'} = ncdouble('temp_time','s_rho','xi_rho') ;
  nc{'temp_south'}.long_name = ncchar('southern boundary potential temperature');
  nc{'temp_south'}.long_name = 'southern boundary potential temperature';
  nc{'temp_south'}.units = ncchar('Celsius');
  nc{'temp_south'}.units = 'Celsius';
  nc{'temp_south'}.coordinates = ncchar('lon_rho s_rho temp_time');
  nc{'temp_south'}.coordinates = 'lon_rho s_rho temp_time';
%
  nc{'salt_south'} = ncdouble('salt_time','s_rho','xi_rho') ;
  nc{'salt_south'}.long_name = ncchar('southern boundary salinity');
  nc{'salt_south'}.long_name = 'southern boundary salinity';
  nc{'salt_south'}.units = ncchar('PSU');
  nc{'salt_south'}.units = 'PSU';
  nc{'salt_south'}.coordinates = ncchar('lon_rho s_rho salt_time');
  nc{'salt_south'}.coordinates = 'lon_rho s_rho salt_time';
%
  nc{'NO3_south'} = ncdouble('NO3_time','s_rho','xi_rho') ;
  nc{'NO3_south'}.long_name = ncchar('southern boundary Nitrate');
  nc{'NO3_south'}.long_name = 'southern boundary Nitrate';
  nc{'NO3_south'}.units = ncchar('milliemole N meter-3');
  nc{'NO3_south'}.units = 'milliemole N meter-3';  
  nc{'NO3_south'}.coordinates = ncchar('lon_rho s_rho NO3_time');
  nc{'NO3_south'}.coordinates = 'lon_rho s_rho NO3_time';

  nc{'NH4_south'} = ncdouble('NH4_time','s_rho','xi_rho') ;
  nc{'NH4_south'}.long_name = ncchar('southern boundary ammonium');
  nc{'NH4_south'}.long_name = 'southern boundary ammonium';
  nc{'NH4_south'}.units = ncchar('milliemole N meter-3');
  nc{'NH4_south'}.units = 'milliemole N meter-3';    
  nc{'NH4_south'}.coordinates = ncchar('lon_rho s_rho NH4_time');
  nc{'NH4_south'}.coordinates = 'lon_rho s_rho NH4_time';
  
  nc{'oxygen_south'} = ncdouble('oxygen_time','s_rho','xi_rho') ;
  nc{'oxygen_south'}.long_name = ncchar('southern boundary oxygen');
  nc{'oxygen_south'}.long_name = 'southern boundary oxygen';
  nc{'oxygen_south'}.units = ncchar('milliemole O meter-3');
  nc{'oxygen_south'}.units = 'milliemole O meter-3';
  nc{'oxygen_south'}.coordinates = ncchar('lon_rho s_rho oxygen_time');
  nc{'oxygen_south'}.coordinates = 'lon_rho s_rho oxygen_time';

  nc{'chlo_south'} = ncdouble('chlo_time','s_rho','xi_rho') ;
  nc{'chlo_south'}.long_name = ncchar('southern boundary chlorophyll');
  nc{'chlo_south'}.long_name = 'southern boundary chlorophyll';
  nc{'chlo_south'}.units = ncchar('milligram meter-3');
  nc{'chlo_south'}.units = 'milligram meter-3';
  nc{'chlo_south'}.coordinates = ncchar('lon_rho s_rho chlo_time');
  nc{'chlo_south'}.coordinates = 'lon_rho s_rho chlo_time';
  
  nc{'tPO4_south'} = ncdouble('tPO4_time','s_rho','xi_rho') ;
  nc{'tPO4_south'}.long_name = ncchar('southern boundary phosphate');
  nc{'tPO4_south'}.long_name = 'southern boundary phosphate';
  nc{'tPO4_south'}.units = ncchar('millimole P meter-3');
  nc{'tPO4_south'}.units = 'millimole P meter-3';
  nc{'tPO4_south'}.coordinates = ncchar('lon_rho s_rho tPO4_time');
  nc{'tPO4_south'}.coordinates = 'lon_rho s_rho tPO4_time';

  nc{'u_south'} = ncdouble('v3d_time','s_rho','xi_u') ;
  nc{'u_south'}.long_name = ncchar('southern boundary u-momentum component');
  nc{'u_south'}.long_name = 'southern boundary u-momentum component';
  nc{'u_south'}.units = ncchar('meter second-1');
  nc{'u_south'}.units = 'meter second-1';
  nc{'u_south'}.coordinates = ncchar('lon_u s_rho uclm_time');
  nc{'u_south'}.coordinates = 'lon_u s_rho u_time';
%
  nc{'v_south'} = ncdouble('v3d_time','s_rho','xi_rho') ;
  nc{'v_south'}.long_name = ncchar('southern boundary v-momentum component');
  nc{'v_south'}.long_name = 'southern boundary v-momentum component';
  nc{'v_south'}.units = ncchar('meter second-1');
  nc{'v_south'}.units = 'meter second-1';
  nc{'v_south'}.coordinates = ncchar('lon_v s_rho vclm_time');
  nc{'v_south'}.coordinates = 'lon_v s_rho vclm_time';
%
  nc{'ubar_south'} = ncdouble('v2d_time','xi_u') ;
  nc{'ubar_south'}.long_name = ncchar('southern boundary vertically integrated u-momentum component');
  nc{'ubar_south'}.long_name = 'southern boundary vertically integrated u-momentum component';
  nc{'ubar_south'}.units = ncchar('meter second-1');
  nc{'ubar_south'}.units = 'meter second-1';
  nc{'ubar_south'}.coordinates = ncchar('lon_u uclm_time');
  nc{'ubar_south'}.coordinates = 'lon_u uclm_time';
%
  nc{'vbar_south'} = ncdouble('v2d_time','xi_rho') ;
  nc{'vbar_south'}.long_name = ncchar('southern boundary vertically integrated v-momentum component');
  nc{'vbar_south'}.long_name = 'southern boundary vertically integrated v-momentum component';
  nc{'vbar_south'}.units = ncchar('meter second-1');
  nc{'vbar_south'}.units = 'meter second-1';
  nc{'vbar_south'}.coordinates = ncchar('lon_v vclm_time');
  nc{'vbar_south'}.coordinates = 'lon_v vclm_time';
%
  nc{'zeta_south'} = ncdouble('zeta_time','xi_rho') ;
  nc{'zeta_south'}.long_name = ncchar('southern boundary sea surface height');
  nc{'zeta_south'}.long_name = 'southern boundary sea surface height';
  nc{'zeta_south'}.units = ncchar('meter');
  nc{'zeta_south'}.units = 'meter';
  nc{'zeta_south'}.coordinates = ncchar('lon_rho zeta_time');
  nc{'zeta_south'}.coordinates = 'lon_rho zeta_time';
%
end
%
if obc(2)==1
%
%   Eastern boundary
%
  nc{'temp_east'} = ncdouble('temp_time','s_rho','eta_rho') ;
  nc{'temp_east'}.long_name = ncchar('eastern boundary potential temperature');
  nc{'temp_east'}.long_name = 'eastern boundary potential temperature';
  nc{'temp_east'}.units = ncchar('Celsius');
  nc{'temp_east'}.units = 'Celsius';
  nc{'temp_east'}.coordinates = ncchar('lat_rho s_rho temp_time');
  nc{'temp_east'}.coordinates = 'lat_rho s_rho temp_time';
%
  nc{'salt_east'} = ncdouble('salt_time','s_rho','eta_rho') ;
  nc{'salt_east'}.long_name = ncchar('eastern boundary salinity');
  nc{'salt_east'}.long_name = 'eastern boundary salinity';
  nc{'salt_east'}.units = ncchar('PSU');
  nc{'salt_east'}.units = 'PSU';
  nc{'salt_east'}.coordinates = ncchar('lat_rho s_rho salt_time');
  nc{'salt_east'}.coordinates = 'lat_rho s_rho salt_time';
%
  nc{'NO3_east'} = ncdouble('NO3_time','s_rho','eta_rho') ;
  nc{'NO3_east'}.long_name = ncchar('eastern boundary Nitrate');
  nc{'NO3_east'}.long_name = 'eastern boundary Nitrate';
  nc{'NO3_east'}.units = ncchar('milliemole N meter-3');
  nc{'NO3_east'}.units = 'milliemole N meter-3';  
  nc{'NO3_east'}.coordinates = ncchar('lon_rho s_rho NO3_time');
  nc{'NO3_east'}.coordinates = 'lon_rho s_rho NO3_time';

  nc{'NH4_east'} = ncdouble('NH4_time','s_rho','eta_rho') ;
  nc{'NH4_east'}.long_name = ncchar('eastern boundary ammonium');
  nc{'NH4_east'}.long_name = 'eastern boundary ammonium';
  nc{'NH4_east'}.units = ncchar('milliemole N meter-3');
  nc{'NH4_east'}.units = 'milliemole N meter-3';    
  nc{'NH4_east'}.coordinates = ncchar('lon_rho s_rho NH4_time');
  nc{'NH4_east'}.coordinates = 'lon_rho s_rho NH4_time';
  
  nc{'oxygen_east'} = ncdouble('oxygen_time','s_rho','eta_rho') ;
  nc{'oxygen_east'}.long_name = ncchar('eastern boundary oxygen');
  nc{'oxygen_east'}.long_name = 'eastern boundary oxygen';
  nc{'oxygen_east'}.units = ncchar('milliemole O meter-3');
  nc{'oxygen_east'}.units = 'milliemole O meter-3';
  nc{'oxygen_east'}.coordinates = ncchar('lon_rho s_rho oxygen_time');
  nc{'oxygen_east'}.coordinates = 'lon_rho s_rho oxygen_time';

  nc{'chlo_east'} = ncdouble('chlo_time','s_rho','eta_rho') ;
  nc{'chlo_east'}.long_name = ncchar('eastern boundary chlorophyll');
  nc{'chlo_east'}.long_name = 'eastern boundary chlorophyll';
  nc{'chlo_east'}.units = ncchar('milligram meter-3');
  nc{'chlo_east'}.units = 'milligram meter-3';
  nc{'chlo_east'}.coordinates = ncchar('lon_rho s_rho chlo_time');
  nc{'chlo_east'}.coordinates = 'lon_rho s_rho chlo_time';
  
  nc{'tPO4_east'} = ncdouble('tPO4_time','s_rho','eta_rho') ;
  nc{'tPO4_east'}.long_name = ncchar('eastern boundary phosphate');
  nc{'tPO4_east'}.long_name = 'eastern boundary phosphate';
  nc{'tPO4_east'}.units = ncchar('millimole P meter-3');
  nc{'tPO4_east'}.units = 'millimole P meter-3';
  nc{'tPO4_east'}.coordinates = ncchar('lon_rho s_rho tPO4_time');
  nc{'tPO4_east'}.coordinates = 'lon_rho s_rho tPO4_time';
  %
  nc{'u_east'} = ncdouble('v3d_time','s_rho','eta_rho') ;
  nc{'u_east'}.long_name = ncchar('eastern boundary u-momentum component');
  nc{'u_east'}.long_name = 'eastern boundary u-momentum component';
  nc{'u_east'}.units = ncchar('meter second-1');
  nc{'u_east'}.units = 'meter second-1';
  nc{'u_east'}.coordinates = ncchar('lat_u s_rho uclm_time');
  nc{'u_east'}.coordinates = 'lat_u s_rho u_time';
%
  nc{'v_east'} = ncdouble('v3d_time','s_rho','eta_v') ;
  nc{'v_east'}.long_name = ncchar('eastern boundary v-momentum component');
  nc{'v_east'}.long_name = 'eastern boundary v-momentum component';
  nc{'v_east'}.units = ncchar('meter second-1');
  nc{'v_east'}.units = 'meter second-1';
  nc{'v_east'}.coordinates = ncchar('lat_v s_rho vclm_time');
  nc{'v_east'}.coordinates = 'lat_v s_rho vclm_time';
%
  nc{'ubar_east'} = ncdouble('v2d_time','eta_rho') ;
  nc{'ubar_east'}.long_name = ncchar('eastern boundary vertically integrated u-momentum component');
  nc{'ubar_east'}.long_name = 'eastern boundary vertically integrated u-momentum component';
  nc{'ubar_east'}.units = ncchar('meter second-1');
  nc{'ubar_east'}.units = 'meter second-1';
  nc{'ubar_east'}.coordinates = ncchar('lat_u uclm_time');
  nc{'ubar_east'}.coordinates = 'lat_u uclm_time';
%
  nc{'vbar_east'} = ncdouble('v2d_time','eta_v') ;
  nc{'vbar_east'}.long_name = ncchar('eastern boundary vertically integrated v-momentum component');
  nc{'vbar_east'}.long_name = 'eastern boundary vertically integrated v-momentum component';
  nc{'vbar_east'}.units = ncchar('meter second-1');
  nc{'vbar_east'}.units = 'meter second-1';
  nc{'vbar_east'}.coordinates = ncchar('lat_v vclm_time');
  nc{'vbar_east'}.coordinates = 'lat_v vclm_time';
%
  nc{'zeta_east'} = ncdouble('zeta_time','eta_rho') ;
  nc{'zeta_east'}.long_name = ncchar('eastern boundary sea surface height');
  nc{'zeta_east'}.long_name = 'eastern boundary sea surface height';
  nc{'zeta_east'}.units = ncchar('meter');
  nc{'zeta_east'}.units = 'meter';
  nc{'zeta_east'}.coordinates = ncchar('lat_rho zeta_time');
  nc{'zeta_east'}.coordinates = 'lat_rho zeta_time';
%
end
%
if obc(3)==1
%
%   Northern boundary
%
  nc{'temp_north'} = ncdouble('temp_time','s_rho','xi_rho') ;
  nc{'temp_north'}.long_name = ncchar('northern boundary potential temperature');
  nc{'temp_north'}.long_name = 'northern boundary potential temperature';
  nc{'temp_north'}.units = ncchar('Celsius');
  nc{'temp_north'}.units = 'Celsius';
  nc{'temp_north'}.coordinates = ncchar('lon_rho s_rho temp_time');
  nc{'temp_north'}.coordinates = 'lon_rho s_rho temp_time';
%
  nc{'salt_north'} = ncdouble('salt_time','s_rho','xi_rho') ;
  nc{'salt_north'}.long_name = ncchar('northern boundary salinity');
  nc{'salt_north'}.long_name = 'northern boundary salinity';
  nc{'salt_north'}.units = ncchar('PSU');
  nc{'salt_north'}.units = 'PSU';
  nc{'salt_north'}.coordinates = ncchar('lon_rho s_rho salt_time');
  nc{'salt_north'}.coordinates = 'lon_rho s_rho salt_time';
%
  nc{'NO3_north'} = ncdouble('NO3_time','s_rho','xi_rho') ;
  nc{'NO3_north'}.long_name = ncchar('northern boundary Nitrate');
  nc{'NO3_north'}.long_name = 'northern boundary Nitrate';
  nc{'NO3_north'}.units = ncchar('milliemole N meter-3');
  nc{'NO3_north'}.units = 'milliemole N meter-3';  
  nc{'NO3_north'}.coordinates = ncchar('lon_rho s_rho NO3_time');
  nc{'NO3_north'}.coordinates = 'lon_rho s_rho NO3_time';

  nc{'NH4_north'} = ncdouble('NH4_time','s_rho','xi_rho') ;
  nc{'NH4_north'}.long_name = ncchar('northern boundary ammonium');
  nc{'NH4_north'}.long_name = 'northern boundary ammonium';
  nc{'NH4_north'}.units = ncchar('milliemole N meter-3');
  nc{'NH4_north'}.units = 'milliemole N meter-3';    
  nc{'NH4_north'}.coordinates = ncchar('lon_rho s_rho NH4_time');
  nc{'NH4_north'}.coordinates = 'lon_rho s_rho NH4_time';
  
  nc{'oxygen_north'} = ncdouble('oxygen_time','s_rho','xi_rho') ;
  nc{'oxygen_north'}.long_name = ncchar('northern boundary oxygen');
  nc{'oxygen_north'}.long_name = 'northern boundary oxygen';
  nc{'oxygen_north'}.units = ncchar('milliemole O meter-3');
  nc{'oxygen_north'}.units = 'milliemole O meter-3';
  nc{'oxygen_north'}.coordinates = ncchar('lon_rho s_rho oxygen_time');
  nc{'oxygen_north'}.coordinates = 'lon_rho s_rho oxygen_time';

  nc{'chlo_north'} = ncdouble('chlo_time','s_rho','xi_rho') ;
  nc{'chlo_north'}.long_name = ncchar('northern boundary chlorophyll');
  nc{'chlo_north'}.long_name = 'northern boundary chlorophyll';
  nc{'chlo_north'}.units = ncchar('milligram meter-3');
  nc{'chlo_north'}.units = 'milligram meter-3';
  nc{'chlo_north'}.coordinates = ncchar('lon_rho s_rho chlo_time');
  nc{'chlo_north'}.coordinates = 'lon_rho s_rho chlo_time';
  
  nc{'tPO4_north'} = ncdouble('tPO4_time','s_rho','xi_rho') ;
  nc{'tPO4_north'}.long_name = ncchar('northern boundary phosphate');
  nc{'tPO4_north'}.long_name = 'northern boundary phosphate';
  nc{'tPO4_north'}.units = ncchar('millimole P meter-3');
  nc{'tPO4_north'}.units = 'millimole P meter-3';
  nc{'tPO4_north'}.coordinates = ncchar('lon_rho s_rho tPO4_time');
  nc{'tPO4_north'}.coordinates = 'lon_rho s_rho tPO4_time';
  %
  nc{'u_north'} = ncdouble('v3d_time','s_rho','xi_u') ;
  nc{'u_north'}.long_name = ncchar('northern boundary u-momentum component');
  nc{'u_north'}.long_name = 'northern boundary u-momentum component';
  nc{'u_north'}.units = ncchar('meter second-1');
  nc{'u_north'}.units = 'meter second-1';
  nc{'u_north'}.coordinates = ncchar('lon_u s_rho uclm_time');
  nc{'u_north'}.coordinates = 'lon_u s_rho u_time';
%
  nc{'v_north'} = ncdouble('v3d_time','s_rho','xi_rho') ;
  nc{'v_north'}.long_name = ncchar('northern boundary v-momentum component');
  nc{'v_north'}.long_name = 'northern boundary v-momentum component';
  nc{'v_north'}.units = ncchar('meter second-1');
  nc{'v_north'}.units = 'meter second-1';
  nc{'v_north'}.coordinates = ncchar('lon_v s_rho vclm_time');
  nc{'v_north'}.coordinates = 'lon_v s_rho vclm_time';
%
  nc{'ubar_north'} = ncdouble('v2d_time','xi_u') ;
  nc{'ubar_north'}.long_name = ncchar('northern boundary vertically integrated u-momentum component');
  nc{'ubar_north'}.long_name = 'northern boundary vertically integrated u-momentum component';
  nc{'ubar_north'}.units = ncchar('meter second-1');
  nc{'ubar_north'}.units = 'meter second-1';
  nc{'ubar_north'}.coordinates = ncchar('lon_u uclm_time');
  nc{'ubar_north'}.coordinates = 'lon_u uclm_time';
%
  nc{'vbar_north'} = ncdouble('v2d_time','xi_rho') ;
  nc{'vbar_north'}.long_name = ncchar('northern boundary vertically integrated v-momentum component');
  nc{'vbar_north'}.long_name = 'northern boundary vertically integrated v-momentum component';
  nc{'vbar_north'}.units = ncchar('meter second-1');
  nc{'vbar_north'}.units = 'meter second-1';
  nc{'vbar_north'}.coordinates = ncchar('lon_v vclm_time');
  nc{'vbar_north'}.coordinates = 'lon_v vclm_time';

  nc{'zeta_north'} = ncdouble('zeta_time','xi_rho') ;
  nc{'zeta_north'}.long_name = ncchar('northern boundary sea surface height');
  nc{'zeta_north'}.long_name = 'northern boundary sea surface height';
  nc{'zeta_north'}.units = ncchar('meter');
  nc{'zeta_north'}.units = 'meter';
  nc{'zeta_north'}.coordinates = ncchar('lon_rho zeta_time');
  nc{'zeta_north'}.coordinates = 'lon_rho zeta_time';
%
end
%
if obc(4)==1
%
%   Western boundary
%
  nc{'temp_west'} = ncdouble('temp_time','s_rho','eta_rho') ;
  nc{'temp_west'}.long_name = ncchar('western boundary potential temperature');
  nc{'temp_west'}.long_name = 'western boundary potential temperature';
  nc{'temp_west'}.units = ncchar('Celsius');
  nc{'temp_west'}.units = 'Celsius';
  nc{'temp_west'}.coordinates = ncchar('lat_rho s_rho temp_time');
  nc{'temp_west'}.coordinates = 'lat_rho s_rho temp_time';
%
  nc{'salt_west'} = ncdouble('salt_time','s_rho','eta_rho') ;
  nc{'salt_west'}.long_name = ncchar('western boundary salinity');
  nc{'salt_west'}.long_name = 'western boundary salinity';
  nc{'salt_west'}.units = ncchar('PSU');
  nc{'salt_west'}.units = 'PSU';
  nc{'salt_west'}.coordinates = ncchar('lat_rho s_rho salt_time');
  nc{'salt_west'}.coordinates = 'lat_rho s_rho salt_time';
%
  nc{'NO3_west'} = ncdouble('NO3_time','s_rho','eta_rho') ;
  nc{'NO3_west'}.long_name = ncchar('western boundary Nitrate');
  nc{'NO3_west'}.long_name = 'western boundary Nitrate';
  nc{'NO3_west'}.units = ncchar('milliemole N meter-3');
  nc{'NO3_west'}.units = 'milliemole N meter-3';  
  nc{'NO3_west'}.coordinates = ncchar('lon_rho s_rho NO3_time');
  nc{'NO3_west'}.coordinates = 'lon_rho s_rho NO3_time';

  nc{'NH4_west'} = ncdouble('NH4_time','s_rho','eta_rho') ;
  nc{'NH4_west'}.long_name = ncchar('western boundary ammonium');
  nc{'NH4_west'}.long_name = 'western boundary ammonium';
  nc{'NH4_west'}.units = ncchar('milliemole N meter-3');
  nc{'NH4_west'}.units = 'milliemole N meter-3';    
  nc{'NH4_west'}.coordinates = ncchar('lon_rho s_rho NH4_time');
  nc{'NH4_west'}.coordinates = 'lon_rho s_rho NH4_time';
  
  nc{'oxygen_west'} = ncdouble('oxygen_time','s_rho','eta_rho') ;
  nc{'oxygen_west'}.long_name = ncchar('western boundary oxygen');
  nc{'oxygen_west'}.long_name = 'western boundary oxygen';
  nc{'oxygen_west'}.units = ncchar('milliemole O meter-3');
  nc{'oxygen_west'}.units = 'milliemole O meter-3';
  nc{'oxygen_west'}.coordinates = ncchar('lon_rho s_rho oxygen_time');
  nc{'oxygen_west'}.coordinates = 'lon_rho s_rho oxygen_time';

  nc{'chlo_west'} = ncdouble('chlo_time','s_rho','eta_rho') ;
  nc{'chlo_west'}.long_name = ncchar('western boundary chlorophyll');
  nc{'chlo_west'}.long_name = 'western boundary chlorophyll';
  nc{'chlo_west'}.units = ncchar('milligram meter-3');
  nc{'chlo_west'}.units = 'milligram meter-3';
  nc{'chlo_west'}.coordinates = ncchar('lon_rho s_rho chlo_time');
  nc{'chlo_west'}.coordinates = 'lon_rho s_rho chlo_time';
  
  nc{'tPO4_west'} = ncdouble('tPO4_time','s_rho','eta_rho') ;
  nc{'tPO4_west'}.long_name = ncchar('western boundary phosphate');
  nc{'tPO4_west'}.long_name = 'western boundary phosphate';
  nc{'tPO4_west'}.units = ncchar('millimole P meter-3');
  nc{'tPO4_west'}.units = 'millimole P meter-3';
  nc{'tPO4_west'}.coordinates = ncchar('lon_rho s_rho tPO4_time');
  nc{'tPO4_west'}.coordinates = 'lon_rho s_rho tPO4_time';
  %
  nc{'u_west'} = ncdouble('v3d_time','s_rho','eta_rho') ;
  nc{'u_west'}.long_name = ncchar('western boundary u-momentum component');
  nc{'u_west'}.long_name = 'western boundary u-momentum component';
  nc{'u_west'}.units = ncchar('meter second-1');
  nc{'u_west'}.units = 'meter second-1';
  nc{'u_west'}.coordinates = ncchar('lat_u s_rho uclm_time');
  nc{'u_west'}.coordinates = 'lat_u s_rho u_time';
%
  nc{'v_west'} = ncdouble('v3d_time','s_rho','eta_v') ;
  nc{'v_west'}.long_name = ncchar('western boundary v-momentum component');
  nc{'v_west'}.long_name = 'western boundary v-momentum component';
  nc{'v_west'}.units = ncchar('meter second-1');
  nc{'v_west'}.units = 'meter second-1';
  nc{'v_west'}.coordinates = ncchar('lat_v s_rho vclm_time');
  nc{'v_west'}.coordinates = 'lat_v s_rho vclm_time';
%
  nc{'ubar_west'} = ncdouble('v2d_time','eta_rho') ;
  nc{'ubar_west'}.long_name = ncchar('western boundary vertically integrated u-momentum component');
  nc{'ubar_west'}.long_name = 'western boundary vertically integrated u-momentum component';
  nc{'ubar_west'}.units = ncchar('meter second-1');
  nc{'ubar_west'}.units = 'meter second-1';
  nc{'ubar_west'}.coordinates = ncchar('lat_u uclm_time');
  nc{'ubar_west'}.coordinates = 'lat_u uclm_time';
%
  nc{'vbar_west'} = ncdouble('v2d_time','eta_v') ;
  nc{'vbar_west'}.long_name = ncchar('western boundary vertically integrated v-momentum component');
  nc{'vbar_west'}.long_name = 'western boundary vertically integrated v-momentum component';
  nc{'vbar_west'}.units = ncchar('meter second-1');
  nc{'vbar_west'}.units = 'meter second-1';
  nc{'vbar_west'}.coordinates = ncchar('lat_v vclm_time');
  nc{'vbar_west'}.coordinates = 'lat_v vclm_time';
%
  nc{'zeta_west'} = ncdouble('zeta_time','eta_rho') ;
  nc{'zeta_west'}.long_name = ncchar('western boundary sea surface height');
  nc{'zeta_west'}.long_name = 'western boundary sea surface height';
  nc{'zeta_west'}.units = ncchar('meter');
  nc{'zeta_west'}.units = 'meter';
  nc{'zeta_west'}.coordinates = ncchar('lat_rho zeta_time');
  nc{'zeta_west'}.coordinates = 'lat_rho zeta_time';
%
end
%
%
% Create global attributes
%
nc.title = ncchar(title);
nc.title = title;
nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(bryname);
nc.clim_file = bryname;
nc.grd_file = ncchar(grdname);
nc.grd_file = grdname;
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
sc_r=((1:N)-N-0.5)/N;
Cs_r=(1.-theta_b)*cff1*sinh(theta_s*sc_r)...
    +theta_b*(cff2*tanh(theta_s*(sc_r+0.5))-0.5);
sc_w=((0:N)-N)/N;
Cs_w=(1.-theta_b)*cff1*sinh(theta_s*sc_w)...
    +theta_b*(cff2*tanh(theta_s*(sc_w+0.5))-0.5);
%
% Write variables
%
nc{'spherical'}(:)='T';
nc{'Vtransform'}(:)=1;
nc{'Vstretching'}(:)=1;
nc{'tstart'}(:) =  min([min(time) min(time) min(time)]); 
nc{'tend'}(:) =  max([max(time) max(time) max(time)]); 
nc{'theta_s'}(:) =  theta_s; 
nc{'theta_b'}(:) =  theta_b; 
nc{'Tcline'}(:) =  hc; 
nc{'hc'}(:) =  hc; 
nc{'s_rho'}(:) = sc_r;
nc{'s_w'}(:) = sc_w;
nc{'Cs_r'}(:) =  Cs_r; 
nc{'Cs_w'}(:) = Cs_w;
nc{'tclm_time'}(:) =  time; 
nc{'temp_time'}(:) =  time;
nc{'NO3_time'}(:) =  time;
nc{'NH4_time'}(:) =  time;
nc{'tPO4_time'}(:) =  time;
nc{'chlo_time'}(:) =  time;
nc{'oxygen_time'}(:) =  time;
nc{'sclm_time'}(:) =  time; 
nc{'salt_time'}(:) =  time; 
nc{'uclm_time'}(:) =  time; 
nc{'vclm_time'}(:) =  time; 
nc{'v2d_time'}(:) =   time; 
nc{'v3d_time'}(:) =   time; 
nc{'ssh_time'}(:) =   time;
nc{'zeta_time'}(:) =  time;
nc{'bry_time'}(:) =  time; 
if obc(1)==1
  nc{'u_south'}(:) =  0; 
  nc{'v_south'}(:) =  0; 
  nc{'ubar_south'}(:) =  0; 
  nc{'vbar_south'}(:) =  0; 
  nc{'zeta_south'}(:) =  0; 
  nc{'temp_south'}(:) =  0; 
  nc{'salt_south'}(:) =  0;
  nc{'NO3_south'}(:) =  0;
  nc{'NH4_south'}(:) =  0;
  nc{'tPO4_south'}(:) =  0;
  nc{'oxygen_south'}(:) =  0;
  nc{'chlo_south'}(:) =  0;
end 
if obc(2)==1
  nc{'u_east'}(:) =  0; 
  nc{'v_east'}(:) =  0; 
  nc{'ubar_east'}(:) =  0; 
  nc{'vbar_east'}(:) =  0; 
  nc{'zeta_east'}(:) =  0; 
  nc{'temp_east'}(:) =  0; 
  nc{'salt_east'}(:) =  0;
  nc{'NO3_east'}(:) =  0;
  nc{'NH4_east'}(:) =  0;
  nc{'tPO4_east'}(:) =  0;
  nc{'oxygen_east'}(:) =  0;
  nc{'chlo_east'}(:) =  0;
end 
if obc(3)==1
  nc{'u_north'}(:) =  0; 
  nc{'v_north'}(:) =  0; 
  nc{'ubar_north'}(:) =  0; 
  nc{'vbar_north'}(:) =  0; 
  nc{'zeta_north'}(:) =  0; 
  nc{'temp_north'}(:) =  0; 
  nc{'salt_north'}(:) =  0;
  nc{'NO3_north'}(:) =  0;
  nc{'NH4_north'}(:) =  0;
  nc{'tPO4_north'}(:) =  0;
  nc{'oxygen_north'}(:) =  0;
  nc{'chlo_north'}(:) =  0;
end 
if obc(4)==1
  nc{'u_west'}(:) =  0; 
  nc{'v_west'}(:) =  0; 
  nc{'ubar_west'}(:) =  0; 
  nc{'vbar_west'}(:) =  0; 
  nc{'zeta_west'}(:) =  0; 
  nc{'temp_west'}(:) =  0; 
  nc{'salt_west'}(:) =  0;
  nc{'NO3_west'}(:) =  0;
  nc{'NH4_west'}(:) =  0;
  nc{'tPO4_west'}(:) =  0;
  nc{'oxygen_west'}(:) =  0;
  nc{'chlo_west'}(:) =  0;
end 
close(nc)
return

