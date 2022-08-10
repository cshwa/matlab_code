function create_bry_Z(zbryname,grdname,title,obc,...
                      Z,time,cycle,clobber);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This function create the header of a Netcdf boundary
%   file (on a z grid).
%
%   Input:
%
%   zbryname     Netcdf climatology file name (character string).
%   grdname      Netcdf grid file name (character string).
%   obc          open boundaries flag (1=open , [S E N W]).
%   Z            Depth of vertical levels.(Vector)
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
disp([' Creating the file : ',zbryname])
disp(' ')
%
%  Read the grid file and check the topography
%
nc = netcdf(grdname, 'nowrite');
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
status=close(nc);
L=Lp-1;
M=Mp-1;
%
%  Create the boundary file
%
type = 'BOUNDARY Z (or OA) file' ; 
history = 'ROMS' ;
nc = netcdf(zbryname,clobber);
result = redef(nc);
%
%  Create dimensions
%
nc('xi_u') = L;
nc('xi_rho') = Lp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
nc('Z') = length(Z);
nc('bry_time') = length(time);
nc('one') = 1;
%
%  Create variables and attributes
%
nc{'Z'} = ncdouble('Z') ;
nc{'Z'}.long_name = ncchar('Depth');
nc{'Z'}.long_name = 'Depth';
nc{'Z'}.units = ncchar('m');
nc{'Z'}.units = 'm';
%
nc{'bry_time'} = ncdouble('bry_time') ;
nc{'bry_time'}.long_name = ncchar('time for temperature climatology');
nc{'bry_time'}.long_name = 'time for temperature climatology';
nc{'bry_time'}.units = ncchar('day');
nc{'bry_time'}.units = 'day';
nc{'bry_time'}.cycle_length = cycle;%
%
if obc(1)==1
%
%   Southern boundary
%
  nc{'temp_south'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'temp_south'}.long_name = ncchar('southern boundary potential temperature');
  nc{'temp_south'}.long_name = 'southern boundary potential temperature';
  nc{'temp_south'}.units = ncchar('Celsius');
  nc{'temp_south'}.units = 'Celsius';
%
  nc{'salt_south'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'salt_south'}.long_name = ncchar('southern boundary salinity');
  nc{'salt_south'}.long_name = 'southern boundary salinity';
  nc{'salt_south'}.units = ncchar('PSU');
  nc{'salt_south'}.units = 'PSU';
  
  nc{'NO3_south'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'NO3_south'}.long_name = ncchar('southern boundary Nitrate');
  nc{'NO3_south'}.long_name = 'southern boundary Nitrate';
  nc{'NO3_south'}.units = ncchar('milliemole N meter-3');
  nc{'NO3_south'}.units = 'milliemole N meter-3';  

  nc{'NH4_south'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'NH4_south'}.long_name = ncchar('southern boundary ammonium');
  nc{'NH4_south'}.long_name = 'southern boundary ammonium';
  nc{'NH4_south'}.units = ncchar('milliemole N meter-3');
  nc{'NH4_south'}.units = 'milliemole N meter-3';    

  nc{'oxygen_south'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'oxygen_south'}.long_name = ncchar('southern boundary oxygen');
  nc{'oxygen_south'}.long_name = 'southern boundary oxygen';
  nc{'oxygen_south'}.units = ncchar('milliemole O meter-3');
  nc{'oxygen_south'}.units = 'milliemole O meter-3';   

  nc{'chlo_south'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'chlo_south'}.long_name = ncchar('southern boundary chlorophyll');
  nc{'chlo_south'}.long_name = 'southern boundary chlorophyll';
  nc{'chlo_south'}.units = ncchar('milligram meter-3');
  nc{'chlo_south'}.units = 'milligram meter-3';
  
  nc{'tPO4_south'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'tPO4_south'}.long_name = ncchar('southern boundary phosphate');
  nc{'tPO4_south'}.long_name = 'southern boundary phosphate';
  nc{'tPO4_south'}.units = ncchar('millimole P meter-3');
  nc{'tPO4_south'}.units = 'millimole P meter-3';  
%
end
%
if obc(2)==1
%
%   Eastern boundary
%
  nc{'temp_east'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'temp_east'}.long_name = ncchar('eastern boundary potential temperature');
  nc{'temp_east'}.long_name = 'eastern boundary potential temperature';
  nc{'temp_east'}.units = ncchar('Celsius');
  nc{'temp_east'}.units = 'Celsius';
%
  nc{'salt_east'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'salt_east'}.long_name = ncchar('eastern boundary salinity');
  nc{'salt_east'}.long_name = 'eastern boundary salinity';
  nc{'salt_east'}.units = ncchar('PSU');
  nc{'salt_east'}.units = 'PSU';
  
  nc{'NO3_east'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'NO3_east'}.long_name = ncchar('eastern boundary Nitrate');
  nc{'NO3_east'}.long_name = 'eastern boundary Nitrate';
  nc{'NO3_east'}.units = ncchar('milliemole N meter-3');
  nc{'NO3_east'}.units = 'milliemole N meter-3';  

  nc{'NH4_east'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'NH4_east'}.long_name = ncchar('eastern boundary ammonium');
  nc{'NH4_east'}.long_name = 'eastern boundary ammonium';
  nc{'NH4_east'}.units = ncchar('milliemole N meter-3');
  nc{'NH4_east'}.units = 'milliemole N meter-3';    

  nc{'oxygen_east'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'oxygen_east'}.long_name = ncchar('eastern boundary oxygen');
  nc{'oxygen_east'}.long_name = 'eastern boundary oxygen';
  nc{'oxygen_east'}.units = ncchar('milliemole O meter-3');
  nc{'oxygen_east'}.units = 'milliemole O meter-3';   

  nc{'chlo_east'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'chlo_east'}.long_name = ncchar('eastern boundary chlorophyll');
  nc{'chlo_east'}.long_name = 'eastern boundary chlorophyll';
  nc{'chlo_east'}.units = ncchar('milligram meter-3');
  nc{'chlo_east'}.units = 'milligram meter-3';
  
  nc{'tPO4_east'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'tPO4_east'}.long_name = ncchar('eastern boundary phosphate');
  nc{'tPO4_east'}.long_name = 'eastern boundary phosphate';
  nc{'tPO4_east'}.units = ncchar('millimole P meter-3');
  nc{'tPO4_east'}.units = 'millimole P meter-3';    
%
end
%
if obc(3)==1
%
%   Northern boundary
%
  nc{'temp_north'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'temp_north'}.long_name = ncchar('northern boundary potential temperature');
  nc{'temp_north'}.long_name = 'northern boundary potential temperature';
  nc{'temp_north'}.units = ncchar('Celsius');
  nc{'temp_north'}.units = 'Celsius';
%
  nc{'salt_north'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'salt_north'}.long_name = ncchar('northern boundary salinity');
  nc{'salt_north'}.long_name = 'northern boundary salinity';
  nc{'salt_north'}.units = ncchar('PSU');
  nc{'salt_north'}.units = 'PSU';
  
  nc{'NO3_north'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'NO3_north'}.long_name = ncchar('northern boundary Nitrate');
  nc{'NO3_north'}.long_name = 'northern boundary Nitrate';
  nc{'NO3_north'}.units = ncchar('milliemole N meter-3');
  nc{'NO3_north'}.units = 'milliemole N meter-3';  

  nc{'NH4_north'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'NH4_north'}.long_name = ncchar('northern boundary ammonium');
  nc{'NH4_north'}.long_name = 'northern boundary ammonium';
  nc{'NH4_north'}.units = ncchar('milliemole N meter-3');
  nc{'NH4_north'}.units = 'milliemole N meter-3';    

  nc{'oxygen_north'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'oxygen_north'}.long_name = ncchar('northern boundary oxygen');
  nc{'oxygen_north'}.long_name = 'northern boundary oxygen';
  nc{'oxygen_north'}.units = ncchar('milliemole O meter-3');
  nc{'oxygen_north'}.units = 'milliemole O meter-3';   

  nc{'chlo_north'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'chlo_north'}.long_name = ncchar('northern boundary chlorophyll');
  nc{'chlo_north'}.long_name = 'northern boundary chlorophyll';
  nc{'chlo_north'}.units = ncchar('milligram meter-3');
  nc{'chlo_north'}.units = 'milligram meter-3';
  
  nc{'tPO4_north'} = ncdouble('bry_time','Z','xi_rho') ;
  nc{'tPO4_north'}.long_name = ncchar('northern boundary phosphate');
  nc{'tPO4_north'}.long_name = 'northern boundary phosphate';
  nc{'tPO4_north'}.units = ncchar('millimole P meter-3');
  nc{'tPO4_north'}.units = 'millimole P meter-3';    
%
end
%
if obc(4)==1
%
%   Western boundary
%
  nc{'temp_west'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'temp_west'}.long_name = ncchar('western boundary potential temperature');
  nc{'temp_west'}.long_name = 'western boundary potential temperature';
  nc{'temp_west'}.units = ncchar('Celsius');
  nc{'temp_west'}.units = 'Celsius';
%
  nc{'salt_west'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'salt_west'}.long_name = ncchar('western boundary salinity');
  nc{'salt_west'}.long_name = 'western boundary salinity';
  nc{'salt_west'}.units = ncchar('PSU');
  nc{'salt_west'}.units = 'PSU';

  nc{'NO3_west'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'NO3_west'}.long_name = ncchar('western boundary Nitrate');
  nc{'NO3_west'}.long_name = 'western boundary Nitrate';
  nc{'NO3_west'}.units = ncchar('milliemole N meter-3');
  nc{'NO3_west'}.units = 'milliemole N meter-3';  

  nc{'NH4_west'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'NH4_west'}.long_name = ncchar('western boundary ammonium');
  nc{'NH4_west'}.long_name = 'western boundary ammonium';
  nc{'NH4_west'}.units = ncchar('milliemole N meter-3');
  nc{'NH4_west'}.units = 'milliemole N meter-3';    

  nc{'oxygen_west'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'oxygen_west'}.long_name = ncchar('western boundary oxygen');
  nc{'oxygen_west'}.long_name = 'western boundary oxygen';
  nc{'oxygen_west'}.units = ncchar('milliemole O meter-3');
  nc{'oxygen_west'}.units = 'milliemole O meter-3';   

  nc{'chlo_west'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'chlo_west'}.long_name = ncchar('western boundary chlorophyll');
  nc{'chlo_west'}.long_name = 'western boundary chlorophyll';
  nc{'chlo_west'}.units = ncchar('milligram meter-3');
  nc{'chlo_west'}.units = 'milligram meter-3';
  
  nc{'tPO4_west'} = ncdouble('bry_time','Z','eta_rho') ;
  nc{'tPO4_west'}.long_name = ncchar('western boundary phosphate');
  nc{'tPO4_west'}.long_name = 'western boundary phosphate';
  nc{'tPO4_west'}.units = ncchar('millimole P meter-3');
  nc{'tPO4_west'}.units = 'millimole P meter-3';    
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
nc.clim_file = ncchar(zbryname);
nc.clim_file = zbryname;
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
% Write variables
%
nc{'Z'}(:) =  Z;
nc{'bry_time'}(:) =  time; 
if obc(1)==1
  nc{'temp_south'}(:) =  0; 
  nc{'salt_south'}(:) =  0;
  nc{'NO3_south'}(:) =  0; 
  nc{'NH4_south'}(:) =  0;  
  nc{'chlo_south'}(:) =  0; 
  nc{'oxygen_south'}(:) =  0;  
  nc{'tPO4_south'}(:) =  0;   
end 
if obc(2)==1 
  nc{'temp_east'}(:) =  0; 
  nc{'salt_east'}(:) =  0;
  nc{'NO3_east'}(:) =  0; 
  nc{'NH4_east'}(:) =  0;  
  nc{'chlo_east'}(:) =  0; 
  nc{'oxygen_east'}(:) =  0;  
  nc{'tPO4_east'}(:) =  0;   
end 
if obc(3)==1 
  nc{'temp_north'}(:) =  0; 
  nc{'salt_north'}(:) =  0;
  nc{'NO3_north'}(:) =  0; 
  nc{'NH4_north'}(:) =  0;  
  nc{'chlo_north'}(:) =  0; 
  nc{'oxygen_north'}(:) =  0;  
  nc{'tPO4_north'}(:) =  0;   
end 
if obc(4)==1 
  nc{'temp_west'}(:) =  0; 
  nc{'salt_west'}(:) =  0;
  nc{'NO3_west'}(:) =  0; 
  nc{'NH4_west'}(:) =  0;  
  nc{'chlo_west'}(:) =  0; 
  nc{'oxygen_west'}(:) =  0;  
  nc{'tPO4_west'}(:) =  0;   
end 
close(nc)
return


