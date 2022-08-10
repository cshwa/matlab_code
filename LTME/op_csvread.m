function [hdr,data] = op_csvread(varargin)
% Reads csv files that contain database information.
%
% data = op_csvread(filename)
%
% INPUT(S):
%
%   filename -- csv file to be read;
%               text files assumed to contain two-line header
%
% OUTPUT(S):
%
%   data -- output matrix (two-dimensional cell array)
%
% EXAMPLE:
%
% OPeNDAP Science Team
% Copyright 2008
% $Version 2.0.0$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2008/01/23 1.0.0 (beta) created, ceb
% 2008/03/01 2.0.0 released, ceb
%==========================================================================

if nargin == 1
    filename = varargin{1};
else
    error('Number of inputs must be 1.')
end

hdr = [];
data = [];

status = exist(filename,'file');
if status == 0
    return
end%endofprogram

%==========================================================================
% READ THE FILE.

if strcmpi(filename(end-2:end),'csv') && status==2
    
    fid = fopen(filename);
    data = textscan(fid,'%s','delimiter',',','treatAsEmpty','-');
    fclose(fid);
    
    %nvar = 106;
    %nrow = nvar+1;
    %ncol = 20;
    
    data = data{1};
    ndata = length(data);
    
    nvar = data{1}; %ignore data{3},data{4},..., data{ncol}
    nvar = str2double(nvar);
    ncol = data{2};
    ncol = str2double(ncol);
    nrow = nvar + 1;
    
    data = data(ncol+1:ndata);
    
    hdr = data(1:ncol); %interesting -- one does not use hdg = data{1:ncol}
    data = data(ncol+1:ncol*nrow);
    data = reshape(data,ncol,nvar); data = data';
       
end%csv-file
%==========================================================================

return %endoffunction