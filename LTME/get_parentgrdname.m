function handles=get_parentgrdname(h,handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get the name of the parent grid file from the GUI
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[filename, pathname] = uigetfile({'*.nc*', 'All netcdf-Files (*.nc*)'; ...
		'*.*','All Files (*.*)'},'PARENT GRID');
handles=reset_handle(h,handles);
if isequal([filename,pathname],[0,0])
  return
end
handles.parentgrid=fullfile(pathname,filename);
nc=netcdf(handles.parentgrid,'nowrite');
if isempty(nc)
  disp('Warning : this is not a netcdf file !!!')
  handles.parentgrid=[];
  return
end
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
h=nc{'h'}(:);
[handles.Mparent,handles.Lparent]=size(lon);
handles.lonmin=min(min(lon));
handles.lonmax=max(max(lon));
handles.latmin=min(min(lat));
handles.latmax=max(max(lat));
handles.hmin=min(min(h));
set(handles.edithmin,'String',num2str(handles.hmin));
close(nc);
%
handles.imin=2;
handles.imax=handles.Lparent-2;
handles.jmin=2;
handles.jmax=handles.Mparent-2;
%
handles=update_plot(h,handles);
%
return
