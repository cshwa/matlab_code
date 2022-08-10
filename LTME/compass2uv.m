function [WVE, WVN] = compass2uv(WMAG, WDIR)
% [WVE, WVN] = compass2uv(WMAG, WDIR)
%
% Abstract:  calculate magnitude and compass direction from 
%   velocity components.
%
%   WMAG = magnitude of velocity vector 
%   WDIR = compass direction of vector (0 deg North).
% 
%   WVE = eastward component of velocity in same units as WMAG
%   WVN = northward component of velocity
%


r=pi/180;

WVE = WMAG .*(sin(WDIR*r));
WVN = WMAG .*(cos(WDIR*r));
