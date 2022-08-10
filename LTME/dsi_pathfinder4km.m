function dsi_pathfinder4km(varargin)
% DataSetInventory and Interface for the Pathfinder4km
%
% Meri Sheremet
% July 2009

% DataSetBranch=function of (Temporal,TAvg,Pass)
List.Temporal={
'Climatology',                'Climatology';
'TimeSeries',                 'Time Series';
};
List.TAvg={
'daily',               '';
'5day',                '';
'7day',                '';
'8day',                '';
'monthly',             '';
%'Seasonal',            '';
%'Annual',              '';
'yearly',              '';
};
List.Pass={
'Day', ''
'Night',''
'Combined',''
};

if nargin < 1
List.DataSetBranch={
'Pathfinder4km_Climatology'
%'Pathfinder4km_Climatology_Daily_Night'
%'Pathfinder4km_Climatology_Daily_Day'
%'Pathfinder4km_Climatology_Daily_Combined'
%'Pathfinder4km_Climatology_5day_Night'
%'Pathfinder4km_Climatology_5day_Day'
%'Pathfinder4km_Climatology_5day_Combined'
%'Pathfinder4km_Climatology_7day_Night'
%'Pathfinder4km_Climatology_7day_Day'
%'Pathfinder4km_Climatology_7day_Combined'
%'Pathfinder4km_Climatology_8day_Night'
%'Pathfinder4km_Climatology_8day_Day'
%'Pathfinder4km_Climatology_8day_Combined'
%'Pathfinder4km_Climatology_Monthly_Night'
%'Pathfinder4km_Climatology_Monthly_Day'
%'Pathfinder4km_Climatology_Monthly_Combined'
%'Pathfinder4km_Climatology_Seasonal_Night'
%'Pathfinder4km_Climatology_Seasonal_Day'
%'Pathfinder4km_Climatology_Seasonal_Combined'
%'Pathfinder4km_Climatology_Annual_Night'
%'Pathfinder4km_Climatology_Annual_Day'
%'Pathfinder4km_Climatology_Annual_Combined'
'Pathfinder4km_daily'
'Pathfinder4km_5day'
'Pathfinder4km_7day'
'Pathfinder4km_8day'
'Pathfinder4km_monthly'
'Pathfinder4km_yearly'
};

    for k=1:length(List.DataSetBranch)
    disp(['DSI_' List.DataSetBranch{k} ])
    eval(['DSI_' List.DataSetBranch{k} ';'])
    end

    dsi_filename = 'dsi_pathfinder4km.mat';
    %if (~isdeployed)
        pname = mfilename('fullpath');
        pname = pname(1:length(pname) - length(mfilename));
        dsi_filename = [pname,dsi_filename];
    %end

    if ~exist(dsi_filename,'file')
        save(dsi_filename,'List')
    else
        save(dsi_filename,'List','-append')
    end

elseif nargin > 1
    disp('Incorrect number of arguments.'); return
else
% with one argument for the specified DataSetBranch only
eval(varargin{1});
end


function DSI=DSI_Pathfinder4km_daily
DSI.DataSetName='Pathfinder4km';DSI.TAvg='daily';  DSI.DataSetBranch=['Pathfinder4km_' DSI.TAvg]; DSI=DSI_common2(DSI);
function DSI=DSI_Pathfinder4km_5day
DSI.DataSetName='Pathfinder4km';DSI.TAvg='5day';   DSI.DataSetBranch=['Pathfinder4km_' DSI.TAvg]; DSI=DSI_common2(DSI);
function DSI=DSI_Pathfinder4km_7day
DSI.DataSetName='Pathfinder4km';DSI.TAvg='7day';   DSI.DataSetBranch=['Pathfinder4km_' DSI.TAvg]; DSI=DSI_common2(DSI);
function DSI=DSI_Pathfinder4km_8day
DSI.DataSetName='Pathfinder4km';DSI.TAvg='8day';   DSI.DataSetBranch=['Pathfinder4km_' DSI.TAvg]; DSI=DSI_common2(DSI);
function DSI=DSI_Pathfinder4km_monthly
DSI.DataSetName='Pathfinder4km';DSI.TAvg='monthly';DSI.DataSetBranch=['Pathfinder4km_' DSI.TAvg]; DSI=DSI_common2(DSI);
function DSI=DSI_Pathfinder4km_yearly
DSI.DataSetName='Pathfinder4km';DSI.TAvg='yearly'; DSI.DataSetBranch=['Pathfinder4km_' DSI.TAvg]; DSI=DSI_common2(DSI);

function DSI=DSI_common2(DSI)
LL={'latitude','longitude'};
LaLo={'lat','lon'};
% qual must be 1st in the Table - it is used for deriving sst_masked
 FieldsTable={
%    opendap name       long name                            name in     returned name    menu name                   opendap dataset CBNum Units
%                                                            dataset                                                  coord   coord
    'qual',            'Quality',                               'qual' , 'qual.qual',     'Quality (qual)',              LL,    LaLo,    9, '';
    'sst',             'Sea Surface Temperature',               'sst'  , 'sst.sst',       'SST',                         LL,    LaLo,    5, 'deg C' ;
    'sst_masked',      'Sea Surface Temperature Masked',        'sst'  , '',              'Land mask',                   LL,    LaLo,    6, '';
    'ssta',            'Sea Surface Temperature Anomaly',       'sst'  , '',              'SST anomaly',                 LL,    LaLo,    7, 'deg C';
    'ssta_masked',     'Sea Surface Temperature Anomaly Masked','sst'  , '',              'Quality mask',                LL,    LaLo,    8, '';
    'bsst',            'First guess SST',                       'bsst' , 'bsst.bsst',     'First guess SST (BSST)',      LL,    LaLo,   10, 'deg C';
    'sdev',            'Standard Deviation',                    'sst' ,  'sst.sst',       'Standard Deviation (sdev)',   LL,    LaLo,   11, 'deg C';
    'num',             'Number of Observations',                'num'  , 'num.num',       'Number of Observations (num)',LL,    LaLo,   13, '';
    'mask1',           'Mask 1',                                'mask1', 'mask1.mask1',   'Mask 1',                      LL,    LaLo,   12, '';
    'mask2',           'Mask 2',                                'mask2', 'mask2.mask2',   'Mask 2',                      LL,    LaLo,   14, '';
    'landmask',        'LandMask',                       'dsp_band_1'  , 'dsp_band_1.dsp_band_1', 'Land mask',           LL,    LaLo,    2, '';
     };

DSI.Fields         =FieldsTable(:,1);
DSI.Fields_NameLong=FieldsTable(:,2);
DSI.Fields_NameDS  =FieldsTable(:,3);
DSI.Fields_NameReturned=FieldsTable(:,4);
DSI.Fields_NameMenu=FieldsTable(:,5);
DSI.Fields_Coord   =FieldsTable(:,6);
DSI.Fields_CoordDS =FieldsTable(:,7);
DSI.Fields_CBNum   =FieldsTable(:,8);
DSI.Fields_Units   =FieldsTable(:,9);

CoordinatesTable={
'latitude','latitude','degrees_north';
'longitude','longitude','degrees_east';
%'depth','depth','m';
'time','time','days since 0000-01-01 00:00:00';
};

% get available lat,lon,time
R.TAvg=DSI.TAvg;R.TimeRange='y';
[DSI.Time,DSI.URLSITE,DSI.URLPATH,DSI.URLFILE]=cat_pathfinder4km(R);
% catalog of files for Pathfinder4km
% based on the file naming conventions described in PFV50_UserGuide.pdf
% input R.TAVG,R.Passes, optional R.DATE1,R.DATE2,R.DATEINCR

%DSI.URLLandMask='http://satdat1.gso.uri.edu/opendap/Pathfinder/landmask/pfv50_land.m04.hdf';
%DSI.URLClimSite='http://satdat1.gso.uri.edu/opendap/sea_surface_temperature/climatology/pathfinder/Version5.0';
DSI.URLLandMask='http://data.nodc.noaa.gov/opendap/pathfinder/Version5.0/pfv50_land.m04.hdf';
DSI.URLClimSite='http://data.nodc.noaa.gov/opendap/pathfinder/Version5.0_Climatologies';

    if 0 % disabled use theoretical values
URL=[DSI.URLSITE DSI.URLPATH.sst{1} DSI.URLFILE.sst{1} '?lat,lon']
loaddap(['+v','-e'],URL);
%dods_err = 0 means no error.
%dods_err = 1 means error.
%dods_err_msg [variable where error message is stored]
if dods_err
    disp(dods_err_msg)
    msgbox('Could not access dataset site. Sorry, try again later.');
    return
end
DSI.Latitude=sst.lat;
DSI.Longitude=sst.lon;
    end
% theoretical values 
d=180/4096;
DSI.Latitude=(90-d/2:-d:-90+d/2)';
DSI.Longitude=(-180+d/2:d:180-d/2)';

DSI.Coordinates=CoordinatesTable(:,1);
DSI.Dimensions.latitude=length(DSI.Latitude);
DSI.Dimensions.longitude=length(DSI.Longitude);
DSI.Dimensions.time=length(DSI.Time);
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
eval(['DSI.Variables.' FIELD '=DSI.Fields_Coord{k};'])
end
DSI.Variables.latitude={'latitude'};
DSI.Variables.longitude={'longitude'};
%DSI.Variables.depth={'depth'};
DSI.Variables.time={'time'};

for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
eval(['DSI.Attributes.' FIELD '.long_name=DSI.Fields_NameLong{k};']);
eval(['DSI.Attributes.' FIELD '.units    =DSI.Fields_Units{k};']);
end
for k=1:length(DSI.Coordinates)
    COORD=DSI.Coordinates{k};
eval(['DSI.Attributes.' COORD '.long_name=CoordinatesTable{k,2};']);
eval(['DSI.Attributes.' COORD '.units    =CoordinatesTable{k,3};']);
end

% Rename returned variables and change dimensions order
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
    FIELDDS=DSI.Fields_NameDS{k}; % sdev is called sst in DataSet

    FREN=['f=squeeze(' DSI.Fields_NameReturned{k} ');'];
eval(['DSI.FormulaRen.' FIELD '=''' FREN ''';']);
    FCNV=['f=A.' FIELDDS '.scale_factor*f+A.' FIELDDS '.add_off;'];
eval(['DSI.FormulaCnv.' FIELD '=''' FCNV ''';']);
end
% bad values 0,1,2,3
%DSI.FormulaNaN='f(find(f<4))=NaN;';
%DSI.FormulaOrderCoord='f=f'';';

DSI.Range.Time1=min(DSI.Time);
DSI.Range.Time2=max(DSI.Time);
DSI.Range.Latitude1=min(DSI.Latitude);
DSI.Range.Latitude2=max(DSI.Latitude);
DSI.Range.Longitude1=min(DSI.Longitude);
DSI.Range.Longitude2=max(DSI.Longitude);

DSI=rmfield(DSI,'Time');DSI=rmfield(DSI,'URLSITE');DSI=rmfield(DSI,'URLPATH');DSI=rmfield(DSI,'URLFILE');
% don't need these; use cat in get

%DSI.Attributes.latitude.units='degrees_north';
%DSI.Attributes.latitude.valid_range=[min(DSI.Latitude) max(DSI.Latitude)];
%DSI.Attributes.longitude.units='degrees_east';
%DSI.Attributes.longitude.valid_range=[min(DSI.Longitude) max(DSI.Longitude)];
%DSI.Attributes.depth.units='m';
%DSI.Attributes.depth.valid_range=[min(DSI.Depth) max(DSI.Depth)];

S1=['time in seconds since 1970-01-01 00:00:00'];
S2=['missing values replaced with NaN'];
S3=['variables converted to physical units;'];
DSI.Readme=[S1 '/' S2 '/' S3];

eval(['DSI_' DSI.DataSetBranch '=DSI;'])

dsi_filename = 'dsi_pathfinder4km.mat';
%if (~isdeployed)
    pname = mfilename('fullpath');
    pname = pname(1:length(pname) - length(mfilename));
    dsi_filename = [pname,dsi_filename];
%end

%if ~exist(dsi_filename,'file')
%    eval(['save ' dsi_filename ' DSI_' DSI.DataSetBranch ';'])
%end
%eval(['save ' dsi_filename ' -append DSI_' DSI.DataSetBranch ';'])

if ~exist(dsi_filename,'file')
save(dsi_filename,['DSI_' DSI.DataSetBranch])
end
save(dsi_filename,['DSI_' DSI.DataSetBranch],'-append')



 % Time Series
%http://satdat1.gso.uri.edu/opendap/Pathfinder/landmask/pfv50_land.m04.hdf?dsp_band_1[1934:2:2047][2958:2:3071]
%http://satdat1.gso.uri.edu/opendap/sea_surface_temperature/climatology/pathfinder/Version5.0//Annual/Day/annual_day.hdf
%http://satdat1.gso.uri.edu/opendap/sea_surface_temperature/climatology/pathfinder/Version5.0//Annual/Day/annual_day.hdf?Clim_SST[1934:2:2047][2958:2:3071]
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/qual/1985005.m04d3pfv50-qual.hdf?qual[1934:2:2047][2958:2:3071]
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/bsst/1985005.s04d3pfv50-bsst-16b.hdf?bsst[1934:2:2047][2958:2:3071]
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/sdev/1985005.s04d3pfv50-sdev-16b.hdf?sst[1934:2:2047][2958:2:3071]
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/msk1/1985005.m04d3pfv50-msk1.hdf?mask1[1934:2:2047][2958:2:3071]
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/msk2/1985005.m04d3pfv50-msk2.hdf?mask2[1934:2:2047][2958:2:3071]
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/nobs/1985005.m04d3pfv50-num.hdf?num[1934:2:2047][2958:2:3071]
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/sst/1985005.s04d3pfv50-sst-16b.hdf?sst[1934:2:2047][2958:2:3071]

%A=loaddap('-A','http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/qual/1985005.m04d3pfv50-qual.hdf?qual[1934:2:2047][2958:2:3071]')
%A =
%                 qual: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.qual
%ans =
%           dsp_PixelType: 1
%           dsp_PixelSize: 1
%                dsp_Flag: 0
%               dsp_nBits: 8
%            dsp_LineSize: 0
%            dsp_cal_name: '"N/A"'
%                   units: '"N/A"'
%       dsp_cal_eqnNumber: 2
%    dsp_cal_CoeffsLength: 8
%          dsp_cal_coeffs: [2x1 double]
%            scale_factor: 1
%                 add_off: 1
%       DODS_ML_Real_Name: 'qual'
%                    qual: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
%A=loaddap('-A','http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/bsst/1985005.s04d3pfv50-bsst-16b.hdf?bsst[1934:2:2047][2958:2:3071]')
%A =
%                 bsst: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.bsst
%ans =
%           dsp_PixelType: 1
%           dsp_PixelSize: 2
%                dsp_Flag: 0
%               dsp_nBits: 16
%            dsp_LineSize: 0
%            dsp_cal_name: '"Temperature"'
%                   units: '"Temp"'
%       dsp_cal_eqnNumber: 2
%    dsp_cal_CoeffsLength: 8
%          dsp_cal_coeffs: [2x1 double]
%            scale_factor: 0.0750
%                 add_off: -3
%       DODS_ML_Real_Name: 'bsst'
%                    bsst: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
%A=loaddap('-A','http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/sdev/1985005.s04d3pfv50-sdev-16b.hdf?sst[1934:2:2047][2958:2:3071]')
%A =
%                  sst: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.sst
%ans =
%           dsp_PixelType: 1
%           dsp_PixelSize: 2
%                dsp_Flag: 0
%               dsp_nBits: 16
%            dsp_LineSize: 0
%            dsp_cal_name: '"Temperature"'
%                   units: '"Temp"'
%       dsp_cal_eqnNumber: 2
%    dsp_cal_CoeffsLength: 8
%          dsp_cal_coeffs: [2x1 double]
%            scale_factor: 0.1500
%                 add_off: 0
%       DODS_ML_Real_Name: 'sst'
%                     sst: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
%A=loaddap('-A','http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/msk1/1985005.m04d3pfv50-msk1.hdf?mask1[1934:2:2047][2958:2:3071]')
%A =
%                mask1: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.mask1
%ans =
%           dsp_PixelType: 1
%           dsp_PixelSize: 1
%                dsp_Flag: 0
%               dsp_nBits: 8
%            dsp_LineSize: 0
%            dsp_cal_name: '"N/A"'
%                   units: '"N/A"'
%       dsp_cal_eqnNumber: 2
%    dsp_cal_CoeffsLength: 8
%          dsp_cal_coeffs: [2x1 double]
%            scale_factor: 1
%                 add_off: 0
%       DODS_ML_Real_Name: 'mask1'
%                   mask1: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
%A=loaddap('-A','http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/msk2/1985005.m04d3pfv50-msk2.hdf?mask2[1934:2:2047][2958:2:3071]')
%A =
%                mask2: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.mask2
%ans =
%           dsp_PixelType: 1
%           dsp_PixelSize: 1
%                dsp_Flag: 0
%               dsp_nBits: 8
%            dsp_LineSize: 0
%            dsp_cal_name: '"N/A"'
%                   units: '"N/A"'
%       dsp_cal_eqnNumber: 2
%    dsp_cal_CoeffsLength: 8
%          dsp_cal_coeffs: [2x1 double]
%            scale_factor: 1
%                 add_off: 0
%       DODS_ML_Real_Name: 'mask2'
%                   mask2: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
%A=loaddap('-A','http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/nobs/1985005.m04d3pfv50-num.hdf?num[1934:2:2047][2958:2:3071]')
%A =
%                  num: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.num
%ans =
%           dsp_PixelType: 1
%           dsp_PixelSize: 1
%                dsp_Flag: 0
%               dsp_nBits: 8
%            dsp_LineSize: 0
%            dsp_cal_name: '"N/A"'
%                   units: '"N/A"'
%       dsp_cal_eqnNumber: 2
%    dsp_cal_CoeffsLength: 8
%          dsp_cal_coeffs: [2x1 double]
%            scale_factor: 1
%                 add_off: 0
%       DODS_ML_Real_Name: 'num'
%                     num: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
%A=loaddap('-A','http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/day/04km/1985/sst/1985005.s04d3pfv50-sst-16b.hdf?sst[1934:2:2047][2958:2:3071]')
%A =
%                  sst: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.sst
%ans =
%           dsp_PixelType: 1
%           dsp_PixelSize: 2
%                dsp_Flag: 0
%               dsp_nBits: 16
%            dsp_LineSize: 0
%            dsp_cal_name: '"Temperature"'
%                   units: '"Temp"'
%       dsp_cal_eqnNumber: 2
%    dsp_cal_CoeffsLength: 8
%          dsp_cal_coeffs: [2x1 double]
%            scale_factor: 0.0750
%                 add_off: -3
%       DODS_ML_Real_Name: 'sst'
%                     sst: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DSI=DSI_Pathfinder4km_Climatology
DSI.DataSetName='Pathfinder4km';
DSI.DataSetBranch='Pathfinder4km_Climatology';
DSI=DSI_clim_common2(DSI);

% Climatology
%http://satdat1.gso.uri.edu/opendap/Pathfinder/landmask/pfv50_land.m04.hdf?dsp_band_1[1934:2:2047][2958:2:3071]
%http://satdat1.gso.uri.edu/opendap/sea_surface_temperature/climatology/pathfinder/Version5.0//Daily/Night/day004_night.hdf
%http://satdat1.gso.uri.edu/opendap/sea_surface_temperature/climatology/pathfinder/Version5.0//Daily/Night/day004_night.hdf?
%Clim_SST[1934:2:2047][2958:2:3071],?Clim_StandardDeviation[1934:2:2047][2958:2:3071],?Clim_Counts[1934:2:2047][2958:2:3071
%Land Mask
%A=loaddap('-A','http://satdat1.gso.uri.edu/opendap/Pathfinder/landmask/pfv50_land.m04.hdf?dsp_band_1[1934:2:2047][2958:2:3071]')


%A =
%           dsp_band_1: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.dsp_band_1
%ans =
%           dsp_PixelType: 1
%           dsp_PixelSize: 1
%                dsp_Flag: 0
%               dsp_nBits: 8
%            dsp_LineSize: 0
%            dsp_cal_name: '"Counts"'
%                   units: '"<no units>"'
%       dsp_cal_eqnNumber: 2
%    dsp_cal_CoeffsLength: 8
%          dsp_cal_coeffs: [2x1 double]
%            scale_factor: 1
%                 add_off: 0
%       DODS_ML_Real_Name: 'dsp_band_1'
%              dsp_band_1: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]

%A=loaddap('-A','http://satdat1.gso.uri.edu/opendap/sea_surface_temperature/climatology/pathfinder/Version5.0//Daily/Day/day005_day.hdf?Clim_SST[1934:2:2047][2958:2:3071],?Clim_StandardDeviation[1934:2:2047][2958:2:3071],?Clim_Counts[1934:2:2047][2950:2:3071]')
%A =
%                  Clim_SST: [1x1 struct]
%    Clim_StandardDeviation: [1x1 struct]
%               Clim_Counts: [1x1 struct]
%         Global_Attributes: [1x1 struct]
%A. Clim_SST
%ans =
%           solar_corr: '"N/A"'
%            limb_corr: '"N/A"'
%       nonlinear_corr: '"N/A"'
%         sst_equation: '"Pathfinder Version 5.0 NLSST"'
%         percent_good: '"N/A"'
%         scale_factor: 0.0750
%     scale_factor_err: 0
%           add_offset: -3
%       add_offset_err: 0
%        calibrated_nt: 23
%            long_name: '"Climatological Mean Sea Surface Temperature"'
%                units: '"Degrees C"'
%               format: '"F5.2"'
%             coordsys: '"geographic"'
%        ml__FillValue: 0
%          valid_range: [2x1 double]
%             C_format: '"%5.2f"'
%        missing_value: 0
%    DODS_ML_Real_Name: 'Clim_SST'
%             Clim_SST: [1x1 struct]
%             Latitude: [1x1 struct]
%            Longitude: [1x1 struct]
%A.Clim_StandardDeviation
%ans =
%                solar_corr: '"N/A"'
%                 limb_corr: '"N/A"'
%            nonlinear_corr: '"N/A"'
%              sst_equation: '"Pathfinder Version 5.0 NLSST"'
%              percent_good: '"N/A"'
%              scale_factor: 0.0750
%          scale_factor_err: 0
%                add_offset: 0
%            add_offset_err: 0
%             calibrated_nt: 23
%                 long_name: '"Standard Deviation of Mean"'
%                     units: '"Degrees C"'
%                    format: '"F5.2"'
%                  coordsys: '"geographic"'
%             ml__FillValue: 0
%               valid_range: [2x1 double]
%                  C_format: '"%5.2f"'
%             missing_value: 0
%         DODS_ML_Real_Name: 'Clim_StandardDeviation'
%    Clim_StandardDeviation: [1x1 struct]
%                  Latitude: [1x1 struct]
%                 Longitude: [1x1 struct]
%A.Clim_Counts
%ans =
%           solar_corr: '"N/A"'
%            limb_corr: '"N/A"'
%       nonlinear_corr: '"N/A"'
%         sst_equation: '"Pathfinder Version 5.0 NLSST"'
%         percent_good: '"N/A"'
%         scale_factor: 1
%     scale_factor_err: 0
%           add_offset: 0
%       add_offset_err: 0
%        calibrated_nt: 21
%            long_name: '"Number of Observations in Mean"'
%                units: '"Counts"'
%               format: '"D3"'
%             coordsys: '"geographic"'
%        ml__FillValue: 0
%          valid_range: [2x1 double]
%             C_format: '"%3d"'
%        missing_value: 0
%    DODS_ML_Real_Name: 'Clim_Counts'
%          Clim_Counts: [1x1 struct]
%             Latitude: [1x1 struct]
%            Longitude: [1x1 struct]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DSI=DSI_clim_common2(DSI)
LL={'latitude','longitude'};
LaLo={'lat','lon'}; 
FieldsTable={
%    opendap                  long name                 name in                  returned name                                    menu name             opendap dataset CBNum
%    name                                               dataset                                                                                         coord   coord
    'Clim_SST',              'Clim_SST',               'Clim_SST',              'Clim_SST.Clim_SST',                             'Clim_SST',              LL,    LaLo,    1,'deg C';
    'Clim_SST_masked',       'Clim_SST_masked',        'N/A',                   'N/A',                                           'Land mask',             LL,    LaLo,    2,'deg C';
    'Clim_StandardDeviation','Clim_StandardDeviation', 'Clim_StandardDeviation','Clim_StandardDeviation.Clim_StandardDeviation', 'Clim_StandardDeviation',LL,    LaLo,    3,'deg C';
    'Clim_Counts',           'Clim_Counts',            'Clim_Counts',           'Clim_Counts.Clim_Counts',                       'Clim_Counts',           LL,    LaLo,    4,'none';
    'landmask',              'LandMask',               'dsp_band_1'  ,          'dsp_band_1.dsp_band_1',                         'Land mask',             LL,    LaLo,    2,'none';
    };

DSI.Fields         =FieldsTable(:,1);
DSI.Fields_NameLong=FieldsTable(:,2);
DSI.Fields_NameDS  =FieldsTable(:,3);
DSI.Fields_NameReturned=FieldsTable(:,4);
DSI.Fields_NameMenu=FieldsTable(:,5);
DSI.Fields_Coord   =FieldsTable(:,6);
DSI.Fields_CoordDS =FieldsTable(:,7);
DSI.Fields_CBNum   =FieldsTable(:,8);
DSI.Fields_Units   =FieldsTable(:,9);

CoordinatesTable={
'latitude','latitude','degrees_north';
'longitude','longitude','degrees_east';
%'depth','depth','m';
'time','time','days since 0000-01-01 00:00:00';
};


% get available lat,lon,time
%R.TAVG=DSI.TAvg;
%[DSI.Time,DSI.URLSITE,DSI.URLPATH,DSI.URLFILE]=cat_pathfinder4km(R);
% catalog of files for Pathfinder4km
% based on the file naming conventions described in PFV50_UserGuide.pdf
% input R.TAVG,R.Passes, optional R.DATE1,R.DATE2,R.DATEINCR

%DSI.URLLandMask='http://satdat1.gso.uri.edu/opendap/Pathfinder/landmask/pfv50_land.m04.hdf';
%DSI.URLClimSite='http://satdat1.gso.uri.edu/opendap/sea_surface_temperature/climatology/pathfinder/Version5.0';
DSI.URLLandMask='http://data.nodc.noaa.gov/opendap/pathfinder/Version5.0/pfv50_land.m04.hdf';
DSI.URLClimSite='http://data.nodc.noaa.gov/opendap/pathfinder/Version5.0_Climatologies';


    if 0 % disabled use theoretical values
URL=[DSI.URLSITE DSI.URLPATH.sst{1} DSI.URLFILE.sst{1} '?lat,lon']
loaddap(['+v','-e'],URL);
%dods_err = 0 means no error.
%dods_err = 1 means error.
%dods_err_msg [variable where error message is stored]
if dods_err
    disp(dods_err_msg)
    msgbox('Could not access dataset site. Sorry, try again later.');
    return
end
DSI.Latitude=sst.lat;
DSI.Longitude=sst.lon;
    end
% theoretical values 
d=180/4096;
DSI.Latitude=(90-d/2:-d:-90+d/2)';
DSI.Longitude=(-180+d/2:d:180-d/2)';

DSI.Coordinates=CoordinatesTable(:,1);
DSI.Dimensions.latitude=length(DSI.Latitude);
DSI.Dimensions.longitude=length(DSI.Longitude);
%DSI.Dimensions.time=length(DSI.Time);
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
eval(['DSI.Variables.' FIELD '=DSI.Fields_Coord{k};'])
end
DSI.Variables.latitude={'latitude'};
DSI.Variables.longitude={'longitude'};
%DSI.Variables.depth={'depth'};
DSI.Variables.time={'time'};

for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
eval(['DSI.Attributes.' FIELD '.long_name=DSI.Fields_NameLong{k};']);
eval(['DSI.Attributes.' FIELD '.units    =DSI.Fields_Units{k};']);
end
for k=1:length(DSI.Coordinates)
    COORD=DSI.Coordinates{k};
eval(['DSI.Attributes.' COORD '.long_name=CoordinatesTable{k,2};']);
eval(['DSI.Attributes.' COORD '.units    =CoordinatesTable{k,3};']);
end

% Rename returned variables and change dimensions order
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
    FIELDDS=DSI.Fields_NameDS{k}; % sdev is called sst in DataSet

    FREN=['f=squeeze(' DSI.Fields_NameReturned{k} ');'];
eval(['DSI.FormulaRen.' FIELD '=''' FREN ''';']);
end
DSI.FormulaCnv.Clim_SST='A.Clim_SST.scale_factor*f+A.Clim_SST.add_offset;';
DSI.FormulaCnv.Clim_StandardDeviation='A.Clim_StandardDeviation.scale_factor*f+A.Clim_StandardDeviation.add_offset;';
DSI.FormulaCnv.Clim_Counts='';
DSI.FormulaCnv.Clim_SST_masked='';

% bad values 0,1,2,3
%DSI.FormulaNaN='f(find(f<4))=NaN;';
%DSI.FormulaOrderCoord='f=f'';';

%DSI.Range.Time1=min(DSI.TimeA);
%DSI.Range.Time2=max(DSI.TimeB);
DSI.Range.Time1=datenum('00000101','yyyymmdd');%min(DSI.Time);
DSI.Range.Time2=datenum('00001231','yyyymmdd');%max(DSI.Time);
DSI.Range.Latitude1=min(DSI.Latitude);
DSI.Range.Latitude2=max(DSI.Latitude);
DSI.Range.Longitude1=min(DSI.Longitude);
DSI.Range.Longitude2=max(DSI.Longitude);

%DSI.Attributes.latitude.units='degrees_north';
%DSI.Attributes.latitude.valid_range=[min(DSI.Latitude) max(DSI.Latitude)];
%DSI.Attributes.longitude.units='degrees_east';
%DSI.Attributes.longitude.valid_range=[min(DSI.Longitude) max(DSI.Longitude)];
%DSI.Attributes.depth.units='m';
%DSI.Attributes.depth.valid_range=[min(DSI.Depth) max(DSI.Depth)];

S1=['time in seconds since 1970-01-01 00:00:00'];
S2=['missing values replaced with NaN'];
S3=['variables converted to physical units;'];
DSI.Readme=[S1 '/' S2 '/' S3];

eval(['DSI_' DSI.DataSetBranch '=DSI;'])

dsi_filename = 'dsi_pathfinder4km.mat';
%if (~isdeployed)
    pname = mfilename('fullpath')
    pname = pname(1:length(pname) - length(mfilename))
    dsi_filename = [pname,dsi_filename]
%end

if ~exist(dsi_filename,'file')
save(dsi_filename,['DSI_' DSI.DataSetBranch])
end
save(dsi_filename,['DSI_' DSI.DataSetBranch],'-append')



