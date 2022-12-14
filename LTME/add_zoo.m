function add_zoo(climfile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function add_zoo(climfile);
%
%  Add zooplancton (mMol N m-3) in a ROMS climatology file.
%  take the chlorophyll (mg C) from the climatology file and
%  multiply by the ratio chlorophyll / phytoplancton derived
%  from previous simulations (Gruber et al., 2005)
%
%  zoo = 0.2 * chla
%
%  input:
%    
%    climfile      : roms climatology file to process (netcdf)
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
%  Updated    2005 by Patrick Marchesiello (IRD)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
theta=0.2;
%
disp('Add_zoo: creating variable and attribute')
%
% open the clim file  
% 
nc=netcdf(climfile,'write');
time= nc{'chla_time'}(:);
cycle= nc{'chla_time'}.cycle_length(:);
tlen=length(time);
redef(nc);
nc('zoo_time') = tlen;
nc{'zoo_time'} = ncdouble('zoo_time') ;
nc{'zooplankton'} = ncdouble('zoo_time','s_rho','eta_rho','xi_rho') ;
%
nc{'zoo_time'}.long_name = ncchar('time for zooplankton');
nc{'zoo_time'}.long_name = 'time for zooplankton';
nc{'zoo_time'}.units = ncchar('day');
nc{'zoo_time'}.units = 'day';
if cycle~=0
  nc{'zoo_time'}.cycle_length = cycle;
end
%
nc{'zooplankton'}.long_name = ncchar('Zooplankton');
nc{'zooplankton'}.long_name = 'Zooplankton';
nc{'zooplankton'}.units = ncchar('mMol N m-3');
nc{'zooplankton'}.units = 'mMol N m-3';
nc{'zooplankton'}.fields = ncchar('ZOO, scalar, series');
nc{'zooplankton'}.fields = 'ZOO, scalar, series';
%
endef(nc);
%
% record the time
%
nc{'zoo_time'}(:)=time;
%
% loop on time
%
for l=1:tlen
  disp(['time index: ',num2str(l),' of total: ',num2str(tlen)])
  nc{'zooplankton'}(l,:,:,:)=theta*squeeze(nc{'chlorophyll'}(l,:,:,:));
end
close(nc);
return
