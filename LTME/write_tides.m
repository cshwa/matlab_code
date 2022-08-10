function write_tides(Tide, Gname, Tname, base_date, Dname)

% WRITE_TIDES:  Creates and writes tidal forcing ROMS NetCDF file
%
% write_tides(Tide, Gname, Tname, base_date, Dname)
%
% This creates ROMS tidal forcing NetCDF file and writes data extracted
% from either OTPS or ADCIRC and processed with the "t_tide" utility.
%
% On Input:
%
%    Tide        Extracted tide data (struct array):
%
%                  Tide.period      Tide angular period (hour)
%                  Tide.names       Tide constituents name (string)
%                  Tide.Ephase      Tide elevation phase angle (degree)
%                  Tide.Eamp        Tide elevation amplitude (m)
%                  Tide.Cmax        Maximum tidal current (m/s)
%                  Tide.Cmin        Minimum tidal current (m/s)
%                  Tide.Cangle      Tide current inclination angle (degree)
%                  Tide.Cphase      Tide current phase angle (degree)
%
%    Gname       ROMS grid NetCDF file name (string)
%            or, an existing grid structure generated by 'get_roms_grid'
%
%    Tname       ROMS tidal forcing NetCDF file name (string)
%
%    base_date   Tidal data reference day (serial date number), for
%                  example:
%
%                  base_date = datenum('01-Jan-2000')
%
%    Dname       Tidal dataset name (string)
%

% svn $Id: write_tides.m 796 2016-05-11 01:45:21Z arango $
%=========================================================================%
%  Copyright (c) 2002-2016 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
  
% Initialize.
  
nctype    = 'nc_double';           % double precision NetCDF variables
Unlimited = true;                  % use unlimited dimension for period

mode = netcdf.getConstant('CLOBBER');                    % overwrite!!!
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

% Get Grid structure.

if (~isstruct(Gname)),
  G = get_roms_grid(Gname);
else
  G = Gname;
end

% Inquire about ROMS grid dimension for Grid NetCDF file.

if (~isstruct(Gname)),
  I = nc_inq(Gname);
else
  I = nc_inq(char(G.Filename))
end

Lr = I.Dimensions(strcmp({I.Dimensions.Name},'xi_rho' )).Length;
Mr = I.Dimensions(strcmp({I.Dimensions.Name},'eta_rho')).Length;

% Check input tide structure.

vargrid = {'spherical','lon_rho','lat_rho'};

varlist = {'period','Ephase','Eamp',                                    ...
           'Cphase','Cangle','Cmin','Cmax'};

roms_vars = [vargrid,                                                   ...
	     'tide_period', 'tide_Ephase', 'tide_Eamp',                 ...
             'tide_Cphase', 'tide_Cangle', 'tide_Cmin', 'tide_Cmax'];

for var = varlist,
  field = char(var);
  if (isfield(Tide,field)),
    D.(field) = size(Tide.(field));
  else
    error(['WRITE_TIDES: unable to find field ''Ephase''',              ...
         ' in structure: Tide']);
  end
end

Nperiod = length(Tide.period);

%--------------------------------------------------------------------------
% Set tidal forcing Metadata structure.
%--------------------------------------------------------------------------

S.Filename = Tname;

% Set global attributes.

S.Attributes(1).Name      = 'type';
S.Attributes(1).Value     = 'ROMS Forcing file';

S.Attributes(2).Name      = 'title';
S.Attributes(2).Value     = ['ROMS Tidal Forcing from ', Dname];

S.Attributes(3).Name      = 'grid_file';
if (isstruct(Gname)),
  S.Attributes(3).Value   = G.Filename;
else
  S.Attributes(3).Value   = Gname;
end

S.Attributes(4).Name      = 'base_date';
S.Attributes(4).Value     = ['days since ' datestr(base_date,31)];

S.Attributes(5).Name      = 'components';
if (iscell(Tide.names)),
  Avalue                  = upper(sprintf('%s, ',Tide.names{:}));
  S.Attributes(5).Value   = Avalue(1:end-2);
else
  S.Attributes(5).Value   = Tide.names;
end

S.Attributes(6).Name      = 'source';
S.Attributes(6).Value     = 'OTPS';

S.Attributes(7).Name      = 'history';
S.Attributes(7).Value     = sprintf('%s:  Created by %s with %s.m',     ...
                                    datestr(now), getenv('USER'),       ...
                                    mfilename);

% Set file dimensions.

S.Dimensions(1).Name      = 'xi_rho';
S.Dimensions(1).Length    = Lr;
S.Dimensions(1).Unlimited = false;

S.Dimensions(2).Name      = 'eta_rho';
S.Dimensions(2).Length    = Mr;
S.Dimensions(2).Unlimited = false;

S.Dimensions(3).Name      = 'tide_period';
S.Dimensions(3).Length    = nc_constant('nc_unlimited');
S.Dimensions(3).Unlimited = true;

% Set file variables.
%
% Notice that it is a good idea to save the tidal period in another
% variable. ROMS only uses those periods that are not zero. Therefore,
% the strategy is to extract several tidal components and zero out the
% undesired components periods for a particular run.  Here, the original
% periods are also stored in another variable so the zero out periods
% can be restored by rewriting the 'tide_period' variable.

ic = 0;

for var = roms_vars,
  vname = char(var);
  ic = ic + 1;
  S.Variables(ic) = roms_metadata(vname, G.spherical, nctype, Unlimited);
  if (strcmp(vname, 'tide_period')),
    ic = ic + 1;
    S.Variables(ic) = roms_metadata(vname, G.spherical, nctype, Unlimited);
    S.Variables(ic).Name = 'period';
    S.Variables(ic).Attributes(1).Value = 'saved tide angular period';
  end
end

% Check ROMS metadata structure.  Fill unassigned fields.

S = check_metadata(S);
  
% Create NetCDF file.

ncid = nc_create(S.Filename, mode, S);

%--------------------------------------------------------------------------
% Write out tide data into NetCDF file.
%--------------------------------------------------------------------------

% Grid variables.

for var = vargrid
  field = char(var);
  status = nc_write(Tname, field, G.(field));
  if (status ~= 0),
    error(['WRITE_TIDES: error while writing ', field, 'into file: ',   ...
           Tname]);
  end
end

% Tide data.

for Rec = 1: Nperiod,
  for var = varlist,
    field = char(var);
    vname = strcat('tide_', field);
    if (strcmp(field, 'period')),
      f = squeeze(Tide.(field)(Rec));
    else
      f = squeeze(Tide.(field)(Rec,:,:));
    end      
    status = nc_write(Tname, vname, f, Rec);
    if (status ~= 0),
      error(['WRITE_TIDES: error while writing ', vname, 'into file: ', ...
             Tname]);
    end
  end
end

return
