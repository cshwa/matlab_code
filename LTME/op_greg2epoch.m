function etime = op_greg2epoch(varargin)
% Converts Gregorian calendar time to seconds past 1970/01/01 00:00:00.
%
% etime = op_greg2epoch(gtime)
%
% INPUT(S):
%   gtime -- 6-element vector (or arrays of six-element vectors) 
%            with elements assumed to be gregorian in content
%            (i.e., [year, mo, day, hr, min, sec])
%
% OUTPUT(S):
%   etime -- seconds since January 1, 1970, 00:00:00 UTC (epoch)
%
% If "gtime" is a matrix, such that it is N rows by 6 columns, then
% etime will be a vector.
%
% OPeNDAP Science Team
% Copyright 2007, 2008
% $Version 2.0.0$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2007/06/28 0.0.x created, ceb
% 2008/03/01 2.0.0 released
%==========================================================================

if nargin == 0
    error('no input arguments')
elseif nargin == 1
    gtime = varargin{1};
else
    error('too many input arguments')
end

s = ndims(gtime);
if s > 2
    error('gtime is larger than two dimensions')
end

[ndat,nele] = size(gtime);
if ndat == 6 && nele == 6
    disp('Program assumes input is N rows by 6 elements.');
end
if nele ~= 6
    if nele ~=3
        gtime = reshape(gtime,nele,ndat);
        tmp = ndat;
        ndat = nele;
        nele = tmp;
    end
end

etime = zeros(1,ndat);
gtime0 = [1970,01,01,00,00,00];
etime = (datenum(gtime) - datenum(gtime0)).*86400;

return%endoffunction