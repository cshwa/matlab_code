function getpot(Vtransform,Vstretching,clmname,grdname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get potential temperature of seawater from insitu
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
%  Copyright (c) 2004-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%  
%  Updated 04-Apr-2018 by Y.Y. Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% open the vertical coordinates parameter
% 
% vert_param;
% romstools_param;
%
% open the grid file  
% 
nc=netcdf(grdname);
h=nc{'h'}(:);
close(nc);

%
% open the clim file 
% 
nc=netcdf(clmname,'write');
theta_s = nc{'theta_s'}(:);
theta_b =  nc{'theta_b'}(:);
hc  =  nc{'hc'}(:);
N =  length(nc('s_rho'));
tlen =  length(nc('tclm_time'));
if tlen==0
  tlen=1;
end
%
% Get the sigma depths
%
P=-1e-4*1025*9.81*zlevs(Vtransform,Vstretching,h,0.*h,theta_s,theta_b,hc,N,'r');
%
% loop on time
%
for l=1:tlen
  disp(['   getpot: Time index: ',num2str(l),' of total: ',num2str(tlen)])
  T=squeeze(nc{'temp'}(l,:,:,:));
  S=squeeze(nc{'salt'}(l,:,:,:));
  nc{'temp'}(l,:,:,:)=theta(S,T,P);
end
close(nc);
return
