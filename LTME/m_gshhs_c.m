function m_gshhs_c(varargin);
% M_GSHHS_C Add a coastline to a given map using the 'crude' resolution of
%           the Global Self-consistant Hierarchical High-resolution 
%           Shorelines.
%
%         M_GSHHS_C( (standard line option,...,...) ) draws the coastline
%         as a simple line.
%         M_GSHHS_C('patch' ( ,standard patch options,...,...) ) draws the 
%         coastline as a number of patches. 
%
%         M_GSHHS_C('save',FILENAME) saves the extracted coastline data
%         for the current projection in a file FILENAME. This allows 
%         speedier replotting using M_USERCOAST(FILENAME). 
%    
%         See also M_PROJ, M_GRID, M_COAST, M_GSHHS_L, M_GSHHS_I, M_GSHHS_H 
%         M_USERCOAST    

% Rich Pawlowicz (rich@ocgy.ubc.ca) 15/June/98
%
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.


FILNAME='private/gshhs_c.b';

% Set current projection to geographic
Currentmap=m_coord('set');
m_coord('geographic');

if length(varargin)>1 & varargin{1}(1)=='s',
  [ncst,Area,k]=mu_coast('c',FILNAME);
  eval(['save ' varargin{2} ' ncst k Area']);
else
  mu_coast('c',FILNAME,varargin{:},'tag','m_gshhs_c');
end;

m_coord(Currentmap.name);

