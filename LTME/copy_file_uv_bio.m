function copy_file_uv_bio(newfile,oldfile,grdname,title,obc,...
                        theta_s,theta_b,hc,N,...
                        time,cycle,year,clobber)
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
disp([' Creating the file : ',newfile])
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
% if hc > hmin
%   error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
% end
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the boundary file
%
type = 'BOUNDARY file' ; 
history = 'ROMS' ;
ncid = netcdf.create(newfile, 'clobber');
%
%  Create dimensions
%
xi_u_dimID = netcdf.defDim(ncid, 'xi_u', L);
xi_v_dimID = netcdf.defDim(ncid, 'xi_v', Lp);
xi_rho_dimID = netcdf.defDim(ncid, 'xi_rho', Lp);
eta_u_dimID = netcdf.defDim(ncid, 'eta_u', Mp);
eta_v_dimID = netcdf.defDim(ncid, 'eta_v', M);
eta_rho_dimID = netcdf.defDim(ncid, 'eta_rho', Mp);
s_rho_dimID = netcdf.defDim(ncid, 's_rho', N);
s_w_dimID = netcdf.defDim(ncid, 's_w', Np);
tracer_dimID = netcdf.defDim(ncid, 'tracer', 2);
bry_time_dimID = netcdf.defDim(ncid, 'bry_time', length(time));
oxygen_time_dimID = netcdf.defDim(ncid, 'oxygen_time', length(time));
NO3_time_dimID = netcdf.defDim(ncid, 'NO3_time', length(time));
NH4_time_dimID = netcdf.defDim(ncid, 'NH4_time', length(time));
chlo_time_dimID = netcdf.defDim(ncid, 'chlo_time', length(time));
phyt_time_dimID = netcdf.defDim(ncid, 'phyt_time', length(time));
zoop_time_dimID = netcdf.defDim(ncid, 'zoop_time', length(time));
LDeN_time_dimID = netcdf.defDim(ncid, 'LDeN_time', length(time));
SDeN_time_dimID = netcdf.defDim(ncid, 'SDeN_time', length(time));
% tclm_time_dimID = netcdf.defDim(ncid, 'tclm_time', length(time));
temp_time_dimID = netcdf.defDim(ncid, 'temp_time', length(time));
% sclm_time_dimID = netcdf.defDim(ncid, 'sclm_time', length(time));
salt_time_dimID = netcdf.defDim(ncid, 'salt_time', length(time));
% uclm_time_dimID = netcdf.defDim(ncid, 'uclm_time', length(time));
% vclm_time_dimID = netcdf.defDim(ncid, 'vclm_time', length(time));
v2d_time_dimID = netcdf.defDim(ncid, 'v2d_time', length(time));
v3d_time_dimID = netcdf.defDim(ncid, 'v3d_time', length(time));
ssh_time_dimID = netcdf.defDim(ncid, 'ssh_time', length(time));
zeta_time_dimID = netcdf.defDim(ncid, 'zeta_time', length(time));
one_dimID = netcdf.defDim(ncid, 'one', 1);

%
%  Create variables and attributes
%
var_ID = netcdf.defVar(ncid, 'spherical', 'char', one_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'grid type logical switch');
netcdf.putAtt(ncid, var_ID, 'flag_values', 'T, F');
netcdf.putAtt(ncid, var_ID, 'flag_meanings', 'spherical Cartesian');
%
var_ID = netcdf.defVar(ncid, 'Vtransform', 'double', one_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'vertical terrain-following transformation equation');
%
var_ID = netcdf.defVar(ncid, 'Vstretching', 'double', one_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'vertical terrain-following stretching function');
%
var_ID = netcdf.defVar(ncid, 'tstart', 'double', one_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'start processing day');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
%
var_ID = netcdf.defVar(ncid, 'tend', 'double', one_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'end processing day');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
%
var_ID = netcdf.defVar(ncid, 'theta_s', 'double', one_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'S-coordinate surface control parameter');
netcdf.putAtt(ncid, var_ID, 'units', 'nondimensional');
%
var_ID = netcdf.defVar(ncid, 'theta_b', 'double', one_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'S-coordinate bottom control parameter');
netcdf.putAtt(ncid, var_ID, 'units', 'nondimensional');
%
var_ID = netcdf.defVar(ncid, 'Tcline', 'double', one_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'S-coordinate surface/bottom layer width');
netcdf.putAtt(ncid, var_ID, 'units', 'meter');
%
var_ID = netcdf.defVar(ncid, 'hc', 'double', one_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'S-coordinate parameter, critical depth');
netcdf.putAtt(ncid, var_ID, 'units', 'meter');
%
var_ID = netcdf.defVar(ncid, 's_rho', 'double', s_rho_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'S-coordinate at RHO-points');
netcdf.putAtt(ncid, var_ID, 'valid_min', '-1.');
netcdf.putAtt(ncid, var_ID, 'valid_max', '0');
netcdf.putAtt(ncid, var_ID, 'positive', 'up');
netcdf.putAtt(ncid, var_ID, 'standard_name', 'ocena_s_coordinate_g2');
netcdf.putAtt(ncid, var_ID, 'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc');
%
var_ID = netcdf.defVar(ncid, 's_w', 'double', s_w_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'S-coordinate at W-points');
netcdf.putAtt(ncid, var_ID, 'valid_min', '-1.');
netcdf.putAtt(ncid, var_ID, 'valid_max', '0');
netcdf.putAtt(ncid, var_ID, 'positive', 'up');
netcdf.putAtt(ncid, var_ID, 'standard_name', 'ocena_s_coordinate_g2');
netcdf.putAtt(ncid, var_ID, 'formula_terms', 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc');
%
var_ID = netcdf.defVar(ncid, 'Cs_r', 'double', s_rho_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'S-coordinate stretching curves at RHO-points');
netcdf.putAtt(ncid, var_ID, 'units', 'nondimensional');
netcdf.putAtt(ncid, var_ID, 'valid_min', '-1.');
netcdf.putAtt(ncid, var_ID, 'valid_max', '0');
%
var_ID = netcdf.defVar(ncid, 'Cs_w', 'double', s_w_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'S-coordinate stretching curves at W-points');
netcdf.putAtt(ncid, var_ID, 'units', 'nondimensional');
netcdf.putAtt(ncid, var_ID, 'valid_min', '-1.');
netcdf.putAtt(ncid, var_ID, 'valid_max', '0');
%
var_ID = netcdf.defVar(ncid, 'bry_time', 'double', bry_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for boundary climatology');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
% var_ID = netcdf.defVar(ncid, 'tclm_time', 'double', tclm_time_dimID);
% netcdf.putAtt(ncid, var_ID, 'long_name', 'time for temperature climatology');
% netcdf.putAtt(ncid, var_ID, 'units', 'day');
% netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'temp_time', 'double', temp_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for temperature climatology');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
% var_ID = netcdf.defVar(ncid, 'sclm_time', 'double', sclm_time_dimID);
% netcdf.putAtt(ncid, var_ID, 'long_name', 'time for salinity climatology');
% netcdf.putAtt(ncid, var_ID, 'units', 'day');
% netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'salt_time', 'double', salt_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for salinity climatology');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'oxygen_time', 'double', oxygen_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for oxygen');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'NH4_time', 'double', NH4_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for Ammonium');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'NO3_time', 'double', NO3_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for Nitrate');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'chlo_time', 'double', chlo_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for chlo');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'phyt_time', 'double', phyt_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for phyt');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'zoop_time', 'double', zoop_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for zoop');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'LDeN_time', 'double', LDeN_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for LDeN');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'SDeN_time', 'double', SDeN_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for SDeN');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
% var_ID = netcdf.defVar(ncid, 'uclm_time', 'double', uclm_time_dimID);
% netcdf.putAtt(ncid, var_ID, 'long_name', 'time climatological u');
% netcdf.putAtt(ncid, var_ID, 'units', 'day');
% netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
% var_ID = netcdf.defVar(ncid, 'vclm_time', 'double', vclm_time_dimID);
% netcdf.putAtt(ncid, var_ID, 'long_name', 'time climatological v');
% netcdf.putAtt(ncid, var_ID, 'units', 'day');
% netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'v2d_time', 'double', v2d_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for 2D velocity climatology');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
%
var_ID = netcdf.defVar(ncid, 'v3d_time', 'double', v3d_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for 3D velocity climatology');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
% var_ID = netcdf.defVar(ncid, 'ssh_time', 'double', ssh_time_dimID);
% netcdf.putAtt(ncid, var_ID, 'long_name', 'time for sea surface height');
% netcdf.putAtt(ncid, var_ID, 'units', 'day');
% netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
var_ID = netcdf.defVar(ncid, 'zeta_time', 'double', zeta_time_dimID);
netcdf.putAtt(ncid, var_ID, 'long_name', 'time for sea surface height');
netcdf.putAtt(ncid, var_ID, 'units', 'day');
netcdf.putAtt(ncid, var_ID, 'cycle_length', cycle);
% %
if obc(1)==1
%
%   Southern boundary
%
  var_ID = netcdf.defVar(ncid, 'temp_south', 'double', [xi_rho_dimID s_rho_dimID temp_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary potential temperature');
  netcdf.putAtt(ncid, var_ID, 'units', 'Celsius');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'temp_time s_rho lon_rho');
%
  var_ID = netcdf.defVar(ncid, 'salt_south', 'double', [xi_rho_dimID s_rho_dimID salt_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary salinity');
  netcdf.putAtt(ncid, var_ID, 'units', 'PSU');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'salt_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'oxygen_south', 'double', [xi_rho_dimID s_rho_dimID oxygen_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary oxygen');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'oxygen_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'NO3_south', 'double', [xi_rho_dimID s_rho_dimID NO3_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary Nitrate');
  netcdf.putAtt(ncid, var_ID, 'units', 'mMol N m-3');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'NO3_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'NH4_south', 'double', [xi_rho_dimID s_rho_dimID NH4_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary Ammonium');
  netcdf.putAtt(ncid, var_ID, 'units', 'mMol N m-3');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'NH4_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'chlo_south', 'double', [xi_rho_dimID s_rho_dimID chlo_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary chlo');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'chlo_time s_rho lon_rho');  
% 
  var_ID = netcdf.defVar(ncid, 'phyt_south', 'double', [xi_rho_dimID s_rho_dimID phyt_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary phyt');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'phyt_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'zoop_south', 'double', [xi_rho_dimID s_rho_dimID zoop_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary zoop');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'zoop_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'LDeN_south', 'double', [xi_rho_dimID s_rho_dimID LDeN_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary LDeN');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'LDeN_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'SDeN_south', 'double', [xi_rho_dimID s_rho_dimID SDeN_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary SDeN');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'SDeN_time s_rho lon_rho');     
%
  var_ID = netcdf.defVar(ncid, 'u_south', 'double', [xi_u_dimID s_rho_dimID v3d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary u-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'uclm_time s_rho lon_u');
%
  var_ID = netcdf.defVar(ncid, 'v_south', 'double', [xi_v_dimID s_rho_dimID v3d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary v-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'vclm_time s_rho lon_v');
%
  var_ID = netcdf.defVar(ncid, 'ubar_south', 'double', [xi_u_dimID v2d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary vertically integrated u-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'uclm_time lon_u');
%
  var_ID = netcdf.defVar(ncid, 'vbar_south', 'double', [xi_v_dimID v2d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary vertically integrated v-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'vclm_time lon_v');

  var_ID = netcdf.defVar(ncid, 'zeta_south', 'double', [xi_rho_dimID zeta_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'southern boundary sea surface height');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'zeta_time lon_rho');
% %
end
%
if obc(2)==1
%
%   Eastern boundary
%
  var_ID = netcdf.defVar(ncid, 'temp_east', 'double', [eta_rho_dimID s_rho_dimID temp_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary potential temperature');
  netcdf.putAtt(ncid, var_ID, 'units', 'Celsius');
    netcdf.putAtt(ncid, var_ID, 'coordinates', 'temp_time s_rho lat_rho');
%
  var_ID = netcdf.defVar(ncid, 'salt_east', 'double', [eta_rho_dimID s_rho_dimID salt_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary salinity');
  netcdf.putAtt(ncid, var_ID, 'units', 'PSU');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'salt_time s_rho lat_rho');
%
  var_ID = netcdf.defVar(ncid, 'oxygen_east', 'double', [eta_rho_dimID s_rho_dimID oxygen_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary oxygen');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'oxygen_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'NO3_east', 'double', [eta_rho_dimID s_rho_dimID NO3_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary Nitrate');
  netcdf.putAtt(ncid, var_ID, 'units', 'mMol N m-3');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'NO3_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'NH4_east', 'double', [eta_rho_dimID s_rho_dimID NH4_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary Ammonium');
  netcdf.putAtt(ncid, var_ID, 'units', 'mMol N m-3');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'NH4_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'chlo_east', 'double', [eta_rho_dimID s_rho_dimID chlo_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary chlo');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'chlo_time s_rho lat_rho');  
% 
  var_ID = netcdf.defVar(ncid, 'phyt_east', 'double', [eta_rho_dimID s_rho_dimID phyt_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary phyt');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'phyt_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'zoop_east', 'double', [eta_rho_dimID s_rho_dimID zoop_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary zoop');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'zoop_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'LDeN_east', 'double', [eta_rho_dimID s_rho_dimID LDeN_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary LDeN');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'LDeN_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'SDeN_east', 'double', [eta_rho_dimID s_rho_dimID SDeN_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary SDeN');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'SDeN_time s_rho lat_rho');     
%
  var_ID = netcdf.defVar(ncid, 'u_east', 'double', [eta_u_dimID s_rho_dimID v3d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary u-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'uclm_time s_rho lat_u');
%
  var_ID = netcdf.defVar(ncid, 'v_east', 'double', [eta_v_dimID s_rho_dimID v3d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary v-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'vclm_time s_rho lat_v');
%
  var_ID = netcdf.defVar(ncid, 'ubar_east', 'double', [eta_u_dimID v2d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary vertically integrated u-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'uclm_time lat_u');
%
  var_ID = netcdf.defVar(ncid, 'vbar_east', 'double', [eta_v_dimID v2d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary vertically integrated v-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'vclm_time lat_v');
%
  var_ID = netcdf.defVar(ncid, 'zeta_east', 'double', [eta_rho_dimID zeta_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'eastern boundary sea surface height');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'zeta_time lat_rho');
% %
end
%
if obc(3)==1
%
%   Northern boundary
%
  var_ID = netcdf.defVar(ncid, 'temp_north', 'double', [xi_rho_dimID s_rho_dimID temp_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary potential temperature');
  netcdf.putAtt(ncid, var_ID, 'units', 'Celsius');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'temp_time s_rho lon_rho');
%
  var_ID = netcdf.defVar(ncid, 'salt_north', 'double', [xi_rho_dimID s_rho_dimID salt_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary salinity');
  netcdf.putAtt(ncid, var_ID, 'units', 'PSU');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'salt_time s_rho lon_rho');
%  
  var_ID = netcdf.defVar(ncid, 'oxygen_north', 'double', [xi_rho_dimID s_rho_dimID oxygen_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary oxygen');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'oxygen_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'NO3_north', 'double', [xi_rho_dimID s_rho_dimID NO3_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary Nitrate');
  netcdf.putAtt(ncid, var_ID, 'units', 'mMol N m-3');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'NO3_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'NH4_north', 'double', [xi_rho_dimID s_rho_dimID NH4_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary Ammonium');
  netcdf.putAtt(ncid, var_ID, 'units', 'mMol N m-3');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'NH4_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'chlo_north', 'double', [xi_rho_dimID s_rho_dimID chlo_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary chlo');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'chlo_time s_rho lon_rho');  
% 
  var_ID = netcdf.defVar(ncid, 'phyt_north', 'double', [xi_rho_dimID s_rho_dimID phyt_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary phyt');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'phyt_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'zoop_north', 'double', [xi_rho_dimID s_rho_dimID zoop_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary zoop');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'zoop_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'LDeN_north', 'double', [xi_rho_dimID s_rho_dimID LDeN_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary LDeN');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'LDeN_time s_rho lon_rho');
% 
  var_ID = netcdf.defVar(ncid, 'SDeN_north', 'double', [xi_rho_dimID s_rho_dimID SDeN_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary SDeN');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'SDeN_time s_rho lon_rho');     
%
  var_ID = netcdf.defVar(ncid, 'u_north', 'double', [xi_u_dimID s_rho_dimID v3d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary u-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'uclm_time s_rho lon_u');
%
  var_ID = netcdf.defVar(ncid, 'v_north', 'double', [xi_v_dimID s_rho_dimID v3d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary v-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'vclm_time s_rho lon_v');
%
  var_ID = netcdf.defVar(ncid, 'ubar_north', 'double', [xi_u_dimID v2d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary vertically integrated u-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'uclm_time lon_u');
%
  var_ID = netcdf.defVar(ncid, 'vbar_north', 'double', [xi_v_dimID v2d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary vertically integrated v-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'vclm_time lon_v');
%
  var_ID = netcdf.defVar(ncid, 'zeta_north', 'double', [xi_rho_dimID zeta_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'northern boundary sea surface height');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'zeta_time lon_rho');
% %
end
%
if obc(4)==1
%
%   Western boundary
%
  var_ID = netcdf.defVar(ncid, 'temp_west', 'double', [eta_rho_dimID s_rho_dimID temp_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary potential temperature');
  netcdf.putAtt(ncid, var_ID, 'units', 'Celsius');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'temp_time s_rho lat_rho');
%
  var_ID = netcdf.defVar(ncid, 'salt_west', 'double', [eta_rho_dimID s_rho_dimID salt_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary salinity');
  netcdf.putAtt(ncid, var_ID, 'units', 'PSU');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'salt_time s_rho lat_rho');
%  
  var_ID = netcdf.defVar(ncid, 'oxygen_west', 'double', [eta_rho_dimID s_rho_dimID oxygen_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary oxygen');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'oxygen_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'NO3_west', 'double', [eta_rho_dimID s_rho_dimID NO3_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary Nitrate');
  netcdf.putAtt(ncid, var_ID, 'units', 'mMol N m-3');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'NO3_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'NH4_west', 'double', [eta_rho_dimID s_rho_dimID NH4_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary Ammonium');
  netcdf.putAtt(ncid, var_ID, 'units', 'mMol N m-3');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'NH4_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'chlo_west', 'double', [eta_rho_dimID s_rho_dimID chlo_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary chlo');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'chlo_time s_rho lat_rho');  
% 
  var_ID = netcdf.defVar(ncid, 'phyt_west', 'double', [eta_rho_dimID s_rho_dimID phyt_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary phyt');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'phyt_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'zoop_west', 'double', [eta_rho_dimID s_rho_dimID zoop_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary zoop');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'zoop_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'LDeN_west', 'double', [eta_rho_dimID s_rho_dimID LDeN_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary LDeN');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'LDeN_time s_rho lat_rho');
% 
  var_ID = netcdf.defVar(ncid, 'SDeN_west', 'double', [eta_rho_dimID s_rho_dimID SDeN_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary SDeN');
  netcdf.putAtt(ncid, var_ID, 'units', 'mg/L');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'SDeN_time s_rho lat_rho');       
%
  var_ID = netcdf.defVar(ncid, 'u_west', 'double', [eta_u_dimID s_rho_dimID v3d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary u-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'uclm_time s_rho lat_u');
%
  var_ID = netcdf.defVar(ncid, 'v_west', 'double', [eta_v_dimID s_rho_dimID v3d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary v-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'vclm_time s_rho lat_v');
%
  var_ID = netcdf.defVar(ncid, 'ubar_west', 'double', [eta_u_dimID v2d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary vertically integrated u-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'uclm_time lat_u');
%
  var_ID = netcdf.defVar(ncid, 'vbar_west', 'double', [eta_v_dimID v2d_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary vertically integrated v-momentum component');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter second-1');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'vclm_time lat_v');
%
  var_ID = netcdf.defVar(ncid, 'zeta_west', 'double', [eta_rho_dimID zeta_time_dimID]);
  netcdf.putAtt(ncid, var_ID, 'long_name', 'western boundary sea surface height');
  netcdf.putAtt(ncid, var_ID, 'units', 'meter');
  netcdf.putAtt(ncid, var_ID, 'coordinates', 'zeta_time lat_rho');
%
end
%
%
% Create global attributes
%
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'title', title);
netcdf.putAtt(ncid,varid,'date',datestr(date, 'yyyymmdd'));
netcdf.putAtt(ncid,varid,'clm_file', newfile);
netcdf.putAtt(ncid,varid,'grd_file', grdname);
netcdf.putAtt(ncid,varid,'type', type);
netcdf.putAtt(ncid,varid,'history', history);
%
% Leave define mode
%
netcdf.endDef(ncid);
netcdf.close(ncid);
% %
% % Compute S coordinates
% %
cff1=1./sinh(theta_s);
cff2=0.5/tanh(0.5*theta_s);
sc_r=((1:N)-N-0.5)/N;
Cs_r=(1.-theta_b)*cff1*sinh(theta_s*sc_r)...
    +theta_b*(cff2*tanh(theta_s*(sc_r+0.5))-0.5);
sc_w=((0:N)-N)/N;
Cs_w=(1.-theta_b)*cff1*sinh(theta_s*sc_w)...
    +theta_b*(cff2*tanh(theta_s*(sc_w+0.5))-0.5);

% % initializing S coordinates
% sc_r(1:N)=0;
% Cs_r(1:N)=0;
% sc_w(1:N+1)=0;
% Cs_w(1:N+1)=0;

%
% Write variables
%
nc = netcdf(newfile, 'write');
nk = netcdf(oldfile,'r');

nc{'spherical'}(:) = 'T';
nc{'Vtransform'}(:) = nk{'Vtransform'}(:);
nc{'Vstretching'}(:) = nk{'Vstretching'}(:);
nc{'tstart'}(:) =  nk{'tstart'}(:); 
nc{'tend'}(:) =  nk{'tend'}(:); 
nc{'theta_s'}(:) =  theta_s; 
nc{'theta_b'}(:) =  theta_b; 
nc{'Tcline'}(:) =  hc; 
nc{'hc'}(:) =  hc; 
nc{'s_rho'}(:) = sc_r;
nc{'s_w'}(:) = sc_w;
nc{'Cs_r'}(:) =  Cs_r; 
nc{'Cs_w'}(:) = Cs_w;
% nc{'tclm_time'}(:) =  time; 
nc{'temp_time'}(:) =  time; 
% nc{'sclm_time'}(:) =  time; 
nc{'salt_time'}(:) =  time; 
nc{'oxygen_time'}(:) =  time; 
nc{'NO3_time'}(:) =  time; 
nc{'NH4_time'}(:) =  time;
nc{'chlo_time'}(:) =  time;
nc{'phyt_time'}(:) =  time; 
nc{'zoop_time'}(:) =  time; 
nc{'LDeN_time'}(:) =  time; 
nc{'SDeN_time'}(:) =  time; 
% nc{'uclm_time'}(:) =  time; 
% nc{'vclm_time'}(:) =  time; 
nc{'v2d_time'}(:) =   time; 
nc{'v3d_time'}(:) =   time; 
% nc{'ssh_time'}(:) =   time;
nc{'zeta_time'}(:) =  time;
nc{'bry_time'}(:) =  time; 

load boundary_bio_input
do_tem=squeeze(do_input(:,:,year-1979,:));
do_sn = permute(do_tem,[3 2 1]);
clearvars do_tem
do_tem=squeeze(do_input_we(:,:,year-1979,:));
do_we = permute(do_tem,[3 2 1]);

no3_sn=permute(no3_input,[3 2 1]);
no3_we =permute(no3_input_we,[3 2 1]);

nh4_sn=permute(nh4_input_sn,[3 2 1]);
nh4_we =permute(nh4_input_we,[3 2 1]);

chl_sn=permute(chl_input_sn,[3 2 1]);
chl_we =permute(chl_input_we,[3 2 1]);

if obc(1)==1
  nc{'u_south'}(:) =  0; 
  nc{'v_south'}(:) =  0; 
  nc{'ubar_south'}(:) =  0; 
  nc{'vbar_south'}(:) =  0; 
  nc{'zeta_south'}(:) =  nk{'zeta_south'}(:); 
  nc{'temp_south'}(:) =  nk{'temp_south'}(:); 
  nc{'salt_south'}(:) =  nk{'salt_south'}(:);
    nc{'oxygen_south'}(:) =  do_sn(:); 
    nc{'NO3_south'}(:) =  no3_sn(:); 
    nc{'NH4_south'}(:) =  nh4_sn(:);
    nc{'chlo_south'}(:) =  chl_sn;
    nc{'phyt_south'}(:) =  0; 
    nc{'zoop_south'}(:) =  0; 
    nc{'LDeN_south'}(:) =  0; 
    nc{'SDeN_south'}(:) =  0;  
end 

if obc(2)==1
  nc{'u_east'}(:) =  0; 
  nc{'v_east'}(:) =  0; 
  nc{'ubar_east'}(:) =  0; 
  nc{'vbar_east'}(:) =  0; 
  nc{'zeta_east'}(:) =  nk{'zeta_east'}(:); 
  nc{'temp_east'}(:) =  nk{'temp_east'}(:); 
  nc{'salt_east'}(:) =  nk{'salt_east'}(:);
    nc{'oxygen_east'}(:) =  do_we(:); 
    nc{'NO3_east'}(:) =  no3_we(:); 
    nc{'NH4_east'}(:) =  nh4_we(:);
    nc{'chlo_east'}(:) =  chl_we;
    nc{'phyt_east'}(:) =  0; 
    nc{'zoop_east'}(:) =  0; 
    nc{'LDeN_east'}(:) =  0; 
    nc{'SDeN_east'}(:) =  0;  
end 

if obc(3)==1
  nc{'u_north'}(:) =  0; 
  nc{'v_north'}(:) =  0; 
  nc{'ubar_north'}(:) =  0; 
  nc{'vbar_north'}(:) =  0; 
  nc{'zeta_north'}(:) =  nk{'zeta_north'}(:); 
  nc{'temp_north'}(:) =  nk{'temp_north'}(:); 
  nc{'salt_north'}(:) =  nk{'salt_north'}(:);
    nc{'oxygen_north'}(:) =  do_sn(:); 
    nc{'NO3_north'}(:) =  no3_sn(:); 
    nc{'NH4_north'}(:) =  nh4_sn(:);
    nc{'chlo_north'}(:) =  chl_sn;
    nc{'phyt_north'}(:) =  0; 
    nc{'zoop_north'}(:) =  0; 
    nc{'LDeN_north'}(:) =  0; 
    nc{'SDeN_north'}(:) =  0;  
  
end 

if obc(4)==1
  nc{'u_west'}(:) =  0; 
  nc{'v_west'}(:) =  0; 
  nc{'ubar_west'}(:) =  0; 
  nc{'vbar_west'}(:) =  0; 
  nc{'zeta_west'}(:) =  nk{'zeta_west'}(:); 
  nc{'temp_west'}(:) =  nk{'temp_west'}(:); 
  nc{'salt_west'}(:) =  nk{'salt_west'}(:);
    nc{'oxygen_east'}(:) =  do_we(:); 
    nc{'NO3_east'}(:) =  no3_we(:); 
    nc{'NH4_east'}(:) =  nh4_we(:);
    nc{'chlo_east'}(:) =  chl_we;
    nc{'phyt_east'}(:) =  0; 
    nc{'zoop_east'}(:) =  0; 
    nc{'LDeN_east'}(:) =  0; 
    nc{'SDeN_east'}(:) =  0;  
end 

close(nc)
return