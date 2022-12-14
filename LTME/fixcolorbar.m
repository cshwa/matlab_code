function fixcolorbar(colpos,colaxis,vname,fsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function fixcolorbar(colpos,colaxis,vname,fsize)
%
% put a colorbar at a fixed position
%
% input:
%
%  colpos   position of the colorbar [left, bottom, width, height]
%  colaxis  = [cmin cmax]: values assigned to the first and 
%           last colors in the current colormap 
%           (default: fit the min and max values of the variable)  
%  vname    name of the variable (string)
%           (default: 'temp')
%  fsize    font size  
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
dc=(colaxis(2)-colaxis(1))/64;
subplot('position',colpos)
x=[0:1];
y=[colaxis(1):dc:colaxis(2)];
[X,Y]=meshgrid(x,y);
pcolor(Y,X,Y)
shading flat
set(gca,'YTick',y,'YTickLabel',[' '])
return
