function handles=get_rcoef(h,handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get the embedding refinment coefficient from the GUI
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
%  Updated    6-Apr-2007 by Pierrick Penven (display MMm instead of M)
%  Updated    6-Apr-2007 by Pierrick Penven (display LLm instead of L)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.rcoeff=round(str2num(get(handles.edit_rcoef,'String')));
set(handles.edit_rcoef,'String',num2str(handles.rcoeff));
handles.Lchild=1+handles.rcoeff*(handles.imax-handles.imin);
set(handles.editLchild,'String',num2str(handles.Lchild-1));
handles.Mchild=1+handles.rcoeff*(handles.jmax-handles.jmin);
set(handles.editMchild,'String',num2str(handles.Mchild-1));
return
