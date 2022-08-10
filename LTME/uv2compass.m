function [WMAG, WDIR] = uv2compass(WVE, WVN)
% [WMAG, WDIR] = uv2compass(WVE, WVN)
%
% Abstract:  calculate magnitude and compass direction from 
%   velocity components.
%
%   WVE = eastward component of velocity
%   WVN = northward component of velocity
%
%   WMAG = magnitude of velocity vector in same units as WVE & WVN
%   WDIR = compass direction of vector (0 deg North).
% 
WMAG = real( sqrt((WVE.^2) + (WVN.^2)) ); 
  
ratio = WVN ./ WVE;
  
% if ratio less than zero, make positive 
whereLESS = find( ratio < 0.0 ); 
ratio(whereLESS) = -1.0 *ratio(whereLESS);

% convert radian to degree 
d = 180/pi; 

% compute wind direction 
WDIR = (atan(ratio)) .* d;

  % determine which quadrant direction wind is pointing and then set
  % compass direction which
  %  
  %    QUADRANTS
  %        |
  %   4th  |  1st
  % ----------------
  %   3rd  |  2nd
  %        |
  % 
  
  % first quad
  whichones = find( WVE >= 0.0 & WVN >= 0.0);
  WDIR(whichones) = 90.0 - WDIR(whichones);
  
  % fourth quad
  whichones = find( WVE < 0.0 & WVN >= 0.0);
  WDIR(whichones) = 270.0 + WDIR(whichones);

  % second quad
  whichones = find( WVE >= 0.0 & WVN < 0.0);
  WDIR(whichones) = 90.0 + WDIR(whichones);
  
  % third quad
  whichones = find( WVE < 0.0 & WVN < 0.0);
  WDIR(whichones) = 270.0 - WDIR(whichones);
  

