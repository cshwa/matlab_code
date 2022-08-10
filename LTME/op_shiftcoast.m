function [Lon, Lat] = op_shiftcoast(long,lat)
% Shifts Matlab's coast.mat data from [-180,180] to [0,360].
%
% coast.mat is part of the mapping toolbox
%
% [Lon,Lat] = op_shiftcoast(long,lat)
%
% INPUT(S):
%
%   long -- longitude vector (obtained from coast.mat)
%   lat -- latitude vector (obtained from coast.mat)
%
% OUTPUT(S):
%
%   Lon -- shifted longitude
%   Lat -- latitude
%
% NOTES:
%
% Matlab's toolbox provides a coastline (load coast). This coastline
% consists of about 10,000 longitude, latitude (long,lat) pairs. Longitude
% ranges from -180 to 180. Fields retrieved with the OPeNDAP Ocean Toolbox 
% cover a variety of ranges, some from -180 to 180, others from 0 to 360 and
% a few neither of these; e.g., global HYCOM fields range from 74 to 434 
% (or something like that, can't remember right now). One can not simply
% add 180 to the Matlab longitudes if one is plotting on a 0 to 360 range
% because of breaks in the outline segments. This routine addresses this 
% problem.
%
%  Returned variables are the longitude and latitude coastline pairs for a
%  0 to 360 map.
%
% OPeNDAP Science Team
% Copyright 2007, 2008
% $Version 2.0.0$

%==========================================================================
% Peter Cornillon
%
% REVISION HISTORY:
% 2007/01/01 0.0.x created, pcornillon
% 2008/03/01 2.0.0 released
% 2008/03/18 2.0.1 modified "help" section
%==========================================================================

j = 0;
imax = length(long);
for i=1:imax-1,

   j = j + 1;
   lonp(j) = long(i);
   latp(j) = lat(i);

   prod = long(i) * long(i+1);
   if prod <= 0
      j = j + 1;
      lonp(j) = nan;
      latp(j) = nan;
   end
   if long(i) == 0 & long(i+1) > 0
      j = j - 1;
      lonp(j) = nan;
      latp(j) = nan;
      j = j + 1;
      lonp(j) = long(i);
      latp(j) = lat(i);
   end
end

Lat = latp;
nn = find(lonp < 0);
Lon = lonp;
Lon(nn) = 360 + lonp(nn);

return

