function h=m_psliceuv(long,lat,w,isub,sca,color);
% M_PSLICEUV Makes a arrow plot on a map (PSLICEUV-style)
%    M_PSLICEUV(LONG,LAT,w) plots velocity vectors as arrows with components 
%    (U+Vi) at the points (LONG,LAT) on the currently defined map.  The 
%    matrices LONG,LAT,w must all be the same size. w (U + Vi) contain the 
%    eastward and northward components of velocity. Arrow scaling is automatic.
% 
%
% PSLICEUV plots a horizontal matrix of 
%          velocity from ECOMSI using arrows
%
%  USAGE: h=psliceuv(x,y,w,isub,sca,color)
% x is array of x points
% y is array of y points 
% w is array of velocities
% isub is number to subsample
% sca is scale factor for arrows
% color is color for arrows
% 
% EXAMPLE:  psliceuv(x,y,w,3,20,'white');
%


global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;



[X,Y]=m_ll2xy(long,lat,'clip','point');

[XN,YN]=m_ll2xy(long,lat+.001,'clip','point');
[XE,YE]=m_ll2xy(long+(.001)./cos(lat*pi/180),lat,'clip','point');
u = real(w);
v = imag(w);
mU=u.*(XE-X)*100 + v.*(XN-X)*100;
mV=u.*(YE-Y)*100 + v.*(YN-Y)*100;
W = (mU+sqrt(-1).*mV);
h=psliceuv(X,Y,W,isub,sca,color);
set(h,'tag','m_psliceuv');

if nargout==0,,
 clear h
end;
