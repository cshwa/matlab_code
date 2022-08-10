function SSTImageOut = op_maskland(varargin)
% Sets SST data to user-specified value (Land Mask-related)
%
% This function will set locations in SSTImageIn for which the LandMask
% value is X to Y.
%
% INPUT(S):
%
%   SSTImageIn -- input structure having fieldname "sst" as one of its
%              fields
%   LandMask -- land mask for data
%   x        -- (optional) values of LandMask for which SST values
%              are set to "y", default = 2
%   y        -- (optional) values to set SST to, default = -5
%
% NOTE: number of input arguments must be 2 or 4.
%
% OPeNDAP Science Team
% Copyright 2008
% $Version 2.0.0$

%==========================================================================
% Peter Cornillon
%
% REVISION HISTORY:
% 2008/01/01 1.0.0 (beta) created, pcornillon
% 2008/03/01 2.0.0 released
% 2008/03/18 2.0.1 modified, ceb
%==========================================================================

x = 2;
y = -5;

if nargin == 2
    
    SSTImageOut = varargin{1};
    LandMask = varargin{2};
    
elseif nargin == 4
    
    SSTImageOut = varargin{1};
    LandMask = varargin{2};
    x = varargin{3};
    y = varargin{4};
    
else
    
    error('Incorrect number of input arguments')
    
end

%SSTImageOut = SSTImageIn;

nn = find(LandMask == x);

if ~isfield(SSTImageOut,'sst')
    
    error('"sst" is not a field')
    
end

SSTImageOut.sst(nn) = y;

return