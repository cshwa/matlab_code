function gtime = op_s2gregorian(sec)
% Converts seconds past 1970 to gregorian.
%
% INPUT(S):
%
%   secs -- seconds since January 1, 1970, 00:00:00 UTC
%
% OUTPUT(S):
%
%   gtime is a six component Gregorian time vector
%                gtime=[year mo da hr mi sec]
%            
% Example: [2003 04 10 10 30 00] = s2gregorian(1049970600)
%
% Note: a pointer-function to op_epoch2greg(...)
%
% OPeNDAP Science Team
% Copyright 2007, 2008
% $Version 2.0.0$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2007/06/01 0.0.x created
% 2008/03/01 2.0.0 released
%==========================================================================

% sec = double(sec);
% days = sec/86400;
% gtime = datevec((datenum([1970,01,01,00,00,00]))+days);

gtime = op_epoch2greg(sec);

return