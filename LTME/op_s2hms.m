function [hr,min,sec]=op_s2hms(secs)
% Converts seconds to interger hour, minute, and seconds.
%
% [hr,min,sec] = op_s2hms(secs)
%
% See documentation below for more information.
%
% INPUT(S):
%
%   secs -- number of seconds
%
% OUTPUT(S):
%
%   hr -- hour
%   min -- minutes
%   sec -- seconds
%
% OPeNDAP Science Team
% Copyright 2007, 2008
% $Version 2.0.0$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2007/01/01 0.0.x edited, ceb
% 2008/02/28 2.0.0 prepare for distribution, ceb
%==========================================================================
% S2HMS: converts seconds to interger hour, minute, and seconds.
% [hr,min,sec]=S2HMS(secs) converts seconds to integer hour, minute,
% and seconds.

%==========================================================================
% 3/11/96: version 1.0
% 8/5/99: version 2.0
%==========================================================================

sec=round(secs);
hr=floor(sec./3600);
min=floor(rem(sec,3600)./60);
sec=round(rem(sec,60));

return%endoffunction