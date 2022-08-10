% To do: uncertain about times of ascending and descending passes in
% GHRSST_AMSRE and GHRSST_TMI:
% according to file info time variable 0 and 24h - this must be a mistake - assumed 0 and 12h.

function dsi_ghrsst(varargin)
% DataSetInventory and Interface for GHRSST 
% DataSetBranch=function of (Level,Product,Region)
List.Level={
% short name for programming, full name for menus;    
'L2P'                         ,'L2P';
'L2P_G'                       ,'L2P_GRIDDED';
'L4'                          ,'L4';
};

List.Product={
% short name for programming, full name for menus;    
'SEVIRI_SST'    ,'SEVIRI_SST'                             ;%L2P         Eastern Atlantic
%'GOES11'        ,'GOES11'                                 ;%L2P         Eastern Pacific
%'GOES12'        ,'GOES12'                                 ;%L2P         Western Atlantic
'AMSRE'         ,'AMSRE'                                  ;%L2P_GRIDDED GLOB
'TMI'           ,'TMI'                                    ;%L2P_GRIDDED Subtropics
'AUS'           ,'AUS'                                    ;%L4           AUS
'NAVO'          ,'NAVO'                                   ;%L4           GLOB
'ODYSSEA'       ,'ODYSSEA'                                ;%L4           GLOB, MED, GAL, NWE
'AVHRR_AMSR_OI' ,'AVHRR_AMSR_OI'                          ;%L4           GLOB
'AVHRR_OI'      ,'AVHRR_OI'                               ;%L4           GLOB
'mw_ir_OI'      ,'mw_ir_OI'                               ;%L4           GLOB
'OSTIA'         ,'OSTIA'                                  ;%L4           GLOB
'MED'           ,'MED'                                    ;%L4           MED
'DMI_OI'        ,'DMI_OI'                                 ;%L4           NSEABALTIC
};

List.Region={
% short name for programming, full name for menus;    
'EATL'                        ,'Eastern Atlantic';
%'EPAC'                        ,'Eastern Pacific';
%'WATL'                        ,'Western Atlantic';
'TROPICS'                     ,'Tropics';           
'GLOB'                        ,'Global Ocean';
'MED'                         ,'Mediterranean';
'NWE'                         ,'Northwest Europe';
'NSEABALTIC'                  ,'North Sea and Baltic';
'AUS'                         ,'Australia';
'GAL'                         ,'Galapagos Islands';
};

if nargin < 1 
% run through all DataSetBranches
    
List.DataSetBranch={    
% short name for programming, full name for menus,      Level         Product         Region T     DSI.TRes   DSI.HRes
'GHRSST_SEVIRI_SST'                  ,'GHRSST_SEVIRI_SST'             ,'L2P'        ,'SEVIRI_SST'   ,'EATL'     ,'daily'  ,'1/10 deg, 11 km';
%'GHRSST_GOES11'                      ,'GHRSST_GOES11'                 ,'L2P'        ,'GOES11'       ,'EPAC'     ,'daily'  ,'6 km';
%'GHRSST_GOES12'                      ,'GHRSST_GOES12'                 ,'L2P'        ,'GOES12'       ,'WATL'     ,'daily'  ,'6 km';
'GHRSST_AMSRE'                       ,'GHRSST_AMSRE'                  ,'L2P_G'      ,'AMSRE'        ,'GLOB'      ,'daily'  ,'1/4 deg, 25 km'; % 27.75 km
'GHRSST_TMI'                         ,'GHRSST_TMI'                    ,'L2P_G'      ,'TMI'          ,'TROPICS'   ,'daily'  ,'1/4 deg, 25 km'; % 27.75 km
'GHRSST_AUS'                         ,'GHRSST_AUS'                    ,'L4'         ,'AUS'          ,'AUS'       ,'daily'  ,'1/12 deg, 9 km'; % 9.25 km
'GHRSST_GAL'                         ,'GHRSST_GAL'                    ,'L4'         ,'ODYSSEA'      ,'GAL'       ,'daily'  ,'1/50 deg, 2 km'; % 2.22 km
'GHRSST_GLOB_EUR'                    ,'GHRSST_GLOB_EUR'               ,'L4'         ,'ODYSSEA'      ,'GLOB'      ,'daily'  ,'1/10 deg, 11 km'; % 11.1 km
'GHRSST_GLOB_NAVO'                   ,'GHRSST_GLOB_NAVO'              ,'L4'         ,'NAVO'         ,'GLOB'      ,'daily'  ,'1/10 deg, 11 km'; % 11.1 km Latitude runs N to S
'GHRSST_GLOB_NCDC_AVHRR_AMSR_OI'     ,'GHRSST_GLOB_NCDC_AVHRR_AMSR_OI','L4'         ,'AVHRR_AMSR_OI','GLOB'      ,'daily'  ,'1/4 deg, 25 km'; % 27.75 km; 
'GHRSST_GLOB_NCDC_AVHRR_OI'          ,'GHRSST_GLOB_NCDC_AVHRR_OI'     ,'L4'         ,'AVHRR_OI'     ,'GLOB'      ,'daily'  ,'1/4 deg, 25 km'; % 27.75 km
'GHRSST_GLOB_REMSS'                  ,'GHRSST_GLOB_REMSS'             ,'L4'         ,'mw_ir_OI'     ,'GLOB'      ,'daily'  ,'1/11 deg, 10 km';
'GHRSST_GLOB_UKMO'                   ,'GHRSST_GLOB_UKMO'              ,'L4'         ,'OSTIA'        ,'GLOB'      ,'daily'  ,'1/20 deg, 6 km'; % 5.55 km
'GHRSST_MED'                         ,'GHRSST_MED'                    ,'L4'         ,'MED'          ,'MED'       ,'daily'  ,'1/50 deg, 2 km'; % 2.22 km
'GHRSST_MED_ODYSSEA'                 ,'GHRSST_MED_ODYSSEA'            ,'L4'         ,'ODYSSEA'      ,'MED'       ,'daily'  ,'1/50 deg, 2 km'; % 2.22 km
'GHRSST_NSEABALTIC'                  ,'GHRSST_NSEABALTIC'             ,'L4'         ,'DMI_OI'       ,'NSEABALTIC','daily'  ,'1/33 deg, 3 km';
'GHRSST_NWE'                         ,'GHRSST_NWE'                    ,'L4'         ,'ODYSSEA'      ,'NWE'       ,'daily'  ,'1/50 deg, 2 km'; % 2.22 km
};

    for k=1:length(List.DataSetBranch)
    disp(['DSI_' List.DataSetBranch{k} ])
    eval(['DSI_' List.DataSetBranch{k} ';'])
    end

    if ~exist('dsi_ghrsst.mat','file')
    save dsi_ghrsst.mat List
    else
    save dsi_ghrsst.mat -append List
    end
 
elseif nargin > 1
disp('Incorrect number of arguments.'); return  
else
% with one argument for the specified DataSetBranch only     
eval(varargin{1}); 
end

% a general request employed by get_DataSet program has the following pattern
% loaddap([URLSITE URLPATH URLFILE '?' URLCVAR URLCTIME URLCLAT URLCLON URLCDEPTH])
% constraints URLCTIME URLCLAT URLCLON URLCDEPTH can be in different order
% constraints URLCLAT URLCLON URLCDEPTH will be generated by program get

%L2P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DSI=DSI_GHRSST_SEVIRI_SST
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_SEVIRI_SST';
%old
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P/SEVIRI_SST/EUR/
%2008/045/20080214-SEVIRI_SST-EUR-L2P-sst3mlml_20080214_0700-v01.nc.bz2?

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P/SEVIRI_SST/EUR/
%2008/092/20080401-SEVIRI_SST-EUR-L2P-sst3mlml_20080401_1900-v01.nc.bz2?
%time[0:1:0],lat[0:1:1200],lon[0:1:1450],
%sst_dtime[0:1:0][0:1:1200][0:1:1450],
%sea_surface_temperature[0:1:0][0:1:1200][0:1:1450],
%sea_ice_fraction[0:1:0][0:1:1200][0:1:1450],
%sources_of_sea_ice_fraction[0:1:0][0:1:1200][0:1:1450],
%DT_analysis[0:1:0][0:1:1200][0:1:1450],
%rejection_flag[0:1:0][0:1:1200][0:1:1450],
%proximity_confidence[0:1:0][0:1:1200][0:1:1450],
%confidence_flag[0:1:0][0:1:1200][0:1:1450],
%SSES_bias_error[0:1:0][0:1:1200][0:1:1450],
%SSES_standard_deviation_error[0:1:0][0:1:1200][0:1:1450],
%aerosol_optical_depth[0:1:0][0:1:1200][0:1:1450],
%aod_dtime_from_sst[0:1:0][0:1:1200][0:1:1450],
%sources_of_aod[0:1:0][0:1:1200][0:1:1450],
%wind_speed[0:1:0][0:1:1200][0:1:1450],
%wind_speed_dtime_from_sst[0:1:0][0:1:1200][0:1:1450],
%sources_of_wind_speed[0:1:0][0:1:1200][0:1:1450],
%surface_solar_irradiance[0:1:0][0:1:1200][0:1:1450],
%ssi_dtime_from_sst[0:1:0][0:1:1200][0:1:1450],
%sources_of_ssi[0:1:0][0:1:1200][0:1:1450]

% define list of fields or variables available in the dataset
LaLo={'latitude','longitude'};
TLaLo={'time','Latitude','Longitude'};
LaLoT={'Latitude','Longitude','time'};

DSI.FieldsTable={
%    opendap    file name           name in                         returned name                                                   menu name                           opendap dataset returned   
%    name                           dataset                                                                                                                              coord  coord   coord    CBNum 
    'sst'            ,''            ,'sea_surface_temperature'      ,'sea_surface_temperature.sea_surface_temperature'             ,'sea surface temperature'            ,LaLo   ,TLaLo  ,LaLoT   ,1 ,'Kelvin';
    'sst_dtime'      ,''            ,'sst_dtime'                    ,'sst_dtime.sst_dtime'                                         ,'time difference from reference time',LaLo   ,TLaLo  ,LaLoT   ,2 ,'second';
    'SSES_bias_error',''            ,'SSES_bias_error'              ,'SSES_bias_error.SSES_bias_error'                             ,'SSES bias error'                    ,LaLo   ,TLaLo  ,LaLoT   ,3 ,'Kelvin';
    'SSES_standard_deviation',''    ,'SSES_standard_deviation_error','SSES_standard_deviation_error.SSES_standard_deviation_error' ,'SSES standard deviation error'      ,LaLo   ,TLaLo  ,LaLoT   ,4 ,'Kelvin';
    'wind_speed'     ,''            ,'wind_speed'                   ,'wind_speed.wind_speed'                                       ,'wind speed'                         ,LaLo   ,TLaLo  ,LaLoT   ,5 ,'m s-1';
    'rejection_flag' ,''            ,'rejection_flag'               ,'rejection_flag.rejection_flag'                               ,'rejection flag'                     ,LaLo   ,TLaLo  ,LaLoT   ,6 ,'';
    'confidence_flag',''            ,'confidence_flag'              ,'confidence_flag.confidence_flag'                             ,'confidence flag'                    ,LaLo   ,TLaLo  ,LaLoT   ,7 ,'';

    'proximity_confidence',''       ,'proximity_confidence'         ,'proximity_confidence.proximity_confidence'                   ,'proximity confidence'               ,LaLo   ,TLaLo  ,LaLoT   ,8 ,'';         
    'sea_ice_fraction',''           ,'sea_ice_fraction'             ,'sea_ice_fraction.sea_ice_fraction'                           ,'sea ice fraction'                   ,LaLo   ,TLaLo  ,LaLoT   ,13,'percent';    
    'sources_of_sea_ice_fraction','','sources_of_sea_ice_fraction'  ,'sources_of_sea_ice_fraction.sources_of_sea_ice_fraction'     ,'sources of sea ice fraction'        ,LaLo   ,TLaLo  ,LaLoT   ,12,'';
    'DT_analysis',''                ,'DT_analysis'                  ,'DT_analysis.DT_analysis'                                     ,'deviation from sst reference'       ,LaLo   ,TLaLo  ,LaLoT   ,9 ,'Kelvin';
    'aerosol_optical_depth',''      ,'aerosol_optical_depth'        ,'aerosol_optical_depth.aerosol_optical_depth'                 ,'aerosol optical depth'              ,LaLo   ,TLaLo  ,LaLoT   ,14,'count';
    'aod_dtime_from_sst',''         ,'aod_dtime_from_sst'           ,'aod_dtime_from_sst.aod_dtime_from_sst'                       ,'aod dtime from sst'                 ,LaLo   ,TLaLo  ,LaLoT   ,15,'hour';
    'sources_of_aod',''             ,'sources_of_aod'               ,'sources_of_aod.sources_of_aod'                               ,'sources of aerosol optical depth'   ,LaLo   ,TLaLo  ,LaLoT   ,16,'';
    'wind_speed_dtime_from_sst',''  ,'wind_speed_dtime_from_sst'    ,'wind_speed_dtime_from_sst.wind_speed_dtime_from_sst'         ,'wind speed dtime from sst'          ,LaLo   ,TLaLo  ,LaLoT   ,17,'hour';
    'sources_of_wind_speed',''      ,'sources_of_wind_speed'        ,'sources_of_wind_speed.sources_of_wind_speed'                 ,'sources of wind speed'              ,LaLo   ,TLaLo  ,LaLoT   ,19,'';
    'surface_solar_irradiance',''   ,'surface_solar_irradiance'     ,'surface_solar_irradiance.surface_solar_irradiance'           ,'surface solar irradiance'           ,LaLo   ,TLaLo  ,LaLoT   ,20,'watt m-2';
    'ssi_dtime_from_sst',''         ,'ssi_dtime_from_sst'           ,'ssi_dtime_from_sst.ssi_dtime_from_sst'                       ,'ssi dtime from sst'                 ,LaLo   ,TLaLo  ,LaLoT   ,21,'hour';
    'sources_of_ssi',''             ,'sources_of_ssi'               ,'sources_of_ssi.sources_of_ssi'                               ,'sources of ssi'                     ,LaLo   ,TLaLo  ,LaLoT   ,22,'';
    };

DSI.Fields=DSI.FieldsTable(:,1);
DSI.Fields_NameLong=DSI.FieldsTable(:,3);
DSI.Fields_NameDS  =DSI.FieldsTable(:,3);
DSI.Fields_NameReturned=DSI.FieldsTable(:,4);
DSI.Fields_NameMenu=DSI.FieldsTable(:,5);
DSI.Fields_Coord   =DSI.FieldsTable(:,6);
DSI.Fields_CoordDS =DSI.FieldsTable(:,7);
DSI.Fields_CoordReturned =DSI.FieldsTable(:,8);
DSI.Fields_CBNum   =DSI.FieldsTable(:,9);
DSI.Fields_Units   =DSI.FieldsTable(:,10);

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P/SEVIRI_SST/EUR/2008/092/
%'20090410-SEVIRI_SST-EUR-L2P-sst3mlml_20090410_0400-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_0400-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_0700-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_1000-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_1300-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_1600-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_1900-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_2200-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080403_0100-v01.nc.bz2';

DSI.Range.Time1_NODC=datenum(2005,01,031);
%DSI.Range.Time1_NODC=datenum(2005,01,032);
%DSI.Range.Time2_NODC=datenum(2009,01,238); % update
DSI.Range.Time2_JPL=now;
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L2P/SEVIRI_SST/EUR';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P/SEVIRI_SST/EUR';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P/SEVIRI_SST/EUR/2008/092/
%'20080401-SEVIRI_SST-EUR-L2P-sst3mlml_20080401_1900-v01.nc.bz2';
%YYYY=datestr(Time-4/24,'yyyy');
%yearday=(Time-4/24-datenum(YYYY,'yyyy')+1);
%DDD=sprintf('%03d',floor(yearday));

% need to evaluate YYYY DDD specially for this DataSetBranch
DSI.URLPATH='YYYY=datestr(Time-4/24,''yyyy'');yearday=(Time-4/24-datenum(YYYY,''yyyy'')+1);DDD=sprintf(''%03d'',floor(yearday));URLPATH=[''/'' datestr(Time-4/24,''yyyy'') ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[datestr(Time-4/24,''yyyymmdd'') ''-SEVIRI_SST-EUR-L2P-sst3mlml_'' datestr(Time,''yyyymmdd_HHMM'') ''-v01.nc.bz2''];';

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 


DSI.TimeofFrame=[4 7 10 13 16 19 22 25]/24;

Nt=length(DSI.Range.Time1:DSI.Range.Time2)*length(DSI.TimeofFrame);
% 1 frame per file
kt=0;
DSI.URLCTIME='[0]';
%DSI.TIME=cell(Nt,1);
DSI.Time=zeros(Nt,1);

for Time=DSI.Range.Time1:DSI.Range.Time2   
    for h=DSI.TimeofFrame
        kt=kt+1;
DSI.Time(kt)=Time+h;
%TIME=datestr(Time,'yyyy-mm-dd HH:MM:SS');
%HH=sprintf('%02d',round(h*24));
%TIME(12:13)=HH;        % to make hour=25
%DSI.TIME{kt}=TIME;      
    end
end

%Time=DSI.Time(21);
%TIME=DSI.TIME{21};
%YYYY=datestr(Time,'yyyy');
%yearday=(Time-datenum(YYYY,'yyyy')+1);
%DDD=sprintf('%03d',floor(yearday));
%YYYYMMDD=datestr(Time,'yyyymmdd');
%YYYYMMDD_HHMM=[YYYYMMDD '_' TIME(12:13) TIME(15:16)];
%URLSITE=DSI.URLSITE; if strfind(URLSITE,'URLSITE=') eval(URLSITE); end
%URLPATH=DSI.URLPATH; if strfind(URLPATH,'URLPATH=') eval(URLPATH); end  
%URLFILE=DSI.URLFILE; if strfind(URLFILE,'URLFILE=') eval(URLFILE); end 

%URL1=[URLSITE URLPATH URLFILE]
%loaddap([URL1 '?lat,lon'])
%NC_GLOBAL.spatial_resolution: 0.1 degree
%NC_GLOBAL.southernmost_latitude: -60.000000000000000
%NC_GLOBAL.northernmost_latitude: 60.000000000000000
%NC_GLOBAL.westernmost_longitude: -100.00000000000000
%NC_GLOBAL.easternmost_longitude: 45.000000000000000
lat=-60:0.1:60;
lon=-100:0.1:45;
DSI.Latitude=lat;
DSI.Longitude=lon;
DSI.Range.Latitude1=min(DSI.Latitude);
DSI.Range.Latitude2=max(DSI.Latitude);
DSI.Range.Longitude1=min(DSI.Longitude);
DSI.Range.Longitude2=max(DSI.Longitude);

% set FormulaNaN and FormulaCnv to standard form
for j=1:length(DSI.Fields);
    FIELD=DSI.Fields{j};
    FIELDDS=DSI.Fields_NameDS{j};
    eval(['DSI.FormulaNaN.' FIELD '=''f(find(f==A.' FIELDDS '.ml__FillValue))=NaN;'';']);
    eval(['DSI.FormulaCnv.' FIELD '=''f=A.' FIELDDS '.scale_factor*f+A.' FIELDDS '.add_offset;'';']);
end
% if they are non standard change them in the calling function

%A=loaddap('-A','http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P/SEVIRI_SST/EUR/2008/180/20080628-SEVIRI_SST-EUR-L2P-sst3mlml_20080628_0400-v01.nc.bz2')
%A = 
%                             time: [1x1 struct]
%                              lat: [1x1 struct]
%                              lon: [1x1 struct]
%                        sst_dtime: [1x1 struct]
%          sea_surface_temperature: [1x1 struct]
%                 sea_ice_fraction: [1x1 struct]
%      sources_of_sea_ice_fraction: [1x1 struct]
%                      DT_analysis: [1x1 struct]
%                   rejection_flag: [1x1 struct]
%             proximity_confidence: [1x1 struct]
%                  confidence_flag: [1x1 struct]
%                  SSES_bias_error: [1x1 struct]
%    SSES_standard_deviation_error: [1x1 struct]
%            aerosol_optical_depth: [1x1 struct]
%               aod_dtime_from_sst: [1x1 struct]
%                   sources_of_aod: [1x1 struct]
%                       wind_speed: [1x1 struct]
%        wind_speed_dtime_from_sst: [1x1 struct]
%            sources_of_wind_speed: [1x1 struct]
%         surface_solar_irradiance: [1x1 struct]
%               ssi_dtime_from_sst: [1x1 struct]
%                   sources_of_ssi: [1x1 struct]
%                Global_Attributes: [1x1 struct]
%A.sst_dtime
%ans = 
%            long_name: '"time difference from reference time"'
%                units: '"second"'
%        ml__FillValue: -32768
%           add_offset: 0
%         scale_factor: 1
%            valid_min: -32767
%            valid_max: 32767
%    DODS_ML_Real_Name: 'sst_dtime'
%            sst_dtime: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]
%A.sea_surface_temperature
%ans = 
%                  long_name: '"sea surface temperature"'
%                      units: '"kelvin"'
%              ml__FillValue: -32768
%                 add_offset: 273.1500
%               scale_factor: 0.0100
%                  valid_min: -500
%                  valid_max: 4000
%                     source: '"EUMETSAT SAF O&SI"'
%          DODS_ML_Real_Name: 'sea_surface_temperature'
%    sea_surface_temperature: [1x1 struct]
%                       time: [1x1 struct]
%                        lat: [1x1 struct]
%                        lon: [1x1 struct]

%            long_name: '"sea ice fraction"'
%                units: '"percent"'
%        ml__FillValue: 128
%           add_offset: 0
%         scale_factor: 0.0100
%            valid_min: 0
%            valid_max: 100
%    DODS_ML_Real_Name: 'sea_ice_fraction'
%     sea_ice_fraction: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]

%                      long_name: '"sources of sea ice fraction"'
%                  ml__FillValue: 128
%                        comment: [1x123 char]
%              DODS_ML_Real_Name: 'sources_of_sea_ice_fraction'
%    sources_of_sea_ice_fraction: [1x1 struct]
%                           time: [1x1 struct]
%                            lat: [1x1 struct]
%                            lon: [1x1 struct]

%            long_name: '"deviation from sst reference"'
%                units: '"kelvin"'
%        ml__FillValue: 128
%           add_offset: 0
%         scale_factor: 0.1000
%            valid_min: 129
%            valid_max: 127
%            reference: '"GHRSST Analysis"'
%    DODS_ML_Real_Name: 'DT_analysis'
%          DT_analysis: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]

%            long_name: '"rejection_flag"'
%              comment: [1x158 char]
%    DODS_ML_Real_Name: 'rejection_flag'
%       rejection_flag: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]
DSI.FormulaNaN.rejection_flag='';
DSI.FormulaCnv.rejection_flag='';

%               long_name: '"proximity confidence value"'
%           ml__FillValue: 128
%       DODS_ML_Real_Name: 'proximity_confidence'
%    proximity_confidence: [1x1 struct]
%                    time: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
DSI.FormulaCnv.proximity_confidence='';

%            long_name: '"confidence flag"'
%              comment: [1x325 char]
%    DODS_ML_Real_Name: 'confidence_flag'
%      confidence_flag: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]
DSI.FormulaNaN.confidence_flag='';
DSI.FormulaCnv.confidence_flag='';

%            long_name: '"SSES bias error based on confidence flags"'
%                units: '"kelvin"'
%        ml__FillValue: 128
%           add_offset: 0
%         scale_factor: 0.0100
%            valid_min: 129
%            valid_max: 127
%    DODS_ML_Real_Name: 'SSES_bias_error'
%      SSES_bias_error: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]

%                        long_name: [1x57 char]
%                            units: '"kelvin"'
%                    ml__FillValue: 128
%                       add_offset: 1
%                     scale_factor: 0.0100
%                        valid_min: 129
%                        valid_max: 127
%                DODS_ML_Real_Name: 'SSES_standard_deviation_error'
%    SSES_standard_deviation_error: [1x1 struct]
%                             time: [1x1 struct]
%                              lat: [1x1 struct]
%                              lon: [1x1 struct]

%                long_name: '"aerosol optical depth"'
%                    units: '"count"'
%            ml__FillValue: 128
%               add_offset: 0
%             scale_factor: 0.1000
%                valid_min: 129
%                valid_max: 127
%        DODS_ML_Real_Name: 'aerosol_optical_depth'
%    aerosol_optical_depth: [1x1 struct]
%                     time: [1x1 struct]
%                      lat: [1x1 struct]
%                      lon: [1x1 struct]

%             long_name: [1x57 char]
%                 units: '"hour"'
%         ml__FillValue: 128
%            add_offset: 0
%          scale_factor: 0.4000
%             valid_min: 129
%             valid_max: 127
%     DODS_ML_Real_Name: 'aod_dtime_from_sst'
%    aod_dtime_from_sst: [1x1 struct]
%                  time: [1x1 struct]
%                   lat: [1x1 struct]
%                   lon: [1x1 struct]

%            long_name: '"sources of aerosol optical depth"'
%        ml__FillValue: 128
%              comment: [1x99 char]
%    DODS_ML_Real_Name: 'sources_of_aod'
%       sources_of_aod: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]
DSI.FormulaCnv.sources_of_aod='';


%            long_name: '"wind speed"'
%                units: '"m s-1"'
%        ml__FillValue: 128
%           add_offset: 0
%         scale_factor: 1
%            valid_min: 129
%            valid_max: 127
%    DODS_ML_Real_Name: 'wind_speed'
%           wind_speed: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]

%                    long_name: [1x64 char]
%                        units: '"hour"'
%                ml__FillValue: 128
%                   add_offset: 0
%                 scale_factor: 0.2000
%                    valid_min: 129
%                    valid_max: 127
%            DODS_ML_Real_Name: 'wind_speed_dtime_from_sst'
%    wind_speed_dtime_from_sst: [1x1 struct]
%                         time: [1x1 struct]
%                          lat: [1x1 struct]
%                          lon: [1x1 struct]

%                long_name: '"sources of wind speed"'
%            ml__FillValue: 128
%                  comment: [1x173 char]
%        DODS_ML_Real_Name: 'sources_of_wind_speed'
%    sources_of_wind_speed: [1x1 struct]
%                     time: [1x1 struct]
%                      lat: [1x1 struct]
%                      lon: [1x1 struct]
DSI.FormulaCnv.sources_of_wind_speed='';


%                   long_name: '"surface solar irradiance"'
%                       units: '"watt m-2"'
%               ml__FillValue: 128
%                  add_offset: 600
%                scale_factor: 5
%                   valid_min: 129
%                   valid_max: 127
%           DODS_ML_Real_Name: 'surface_solar_irradiance'
%    surface_solar_irradiance: [1x1 struct]
%                        time: [1x1 struct]
%                         lat: [1x1 struct]
%                         lon: [1x1 struct]

%             long_name: '"time difference of surface solar irradiance
%                          measurement from sst measurement"'
%                 units: '"hour"'
%         ml__FillValue: 128
%            add_offset: 0
%          scale_factor: 0.2000
%             valid_min: 129
%             valid_max: 127
%     DODS_ML_Real_Name: 'ssi_dtime_from_sst'
%    ssi_dtime_from_sst: [1x1 struct]
%                  time: [1x1 struct]
%                   lat: [1x1 struct]
%                   lon: [1x1 struct]

%            long_name: '"sources of surface solar irradiance"'
%        ml__FillValue: 128
%              comment: [1x89 char]
%    DODS_ML_Real_Name: 'sources_of_ssi'
%       sources_of_ssi: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]
DSI.FormulaCnv.sources_of_ssi='';


S1=['time in seconds since reference time;'];
S2=['missing values replaced with NaN;'];
S3=['variables converted to physical units.'];
%S4='';%['time is MATLAB datenum.'];
DSI.Readme=[S1 '/' S2 '/' S3];
DSI=DSI_common3(DSI);
SaveDSI(DSI)
        
function DSI=DSI_GHRSST_GOES11
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_GOES11';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P/GOES11/OSDPD/
%2008/043/20080212-GOES11-OSDPD-L2P-GOES11_North_1000Z-v01.nc.gz?
%lat[0:1:1326][0:1:3399],lon[0:1:1326][0:1:3399],time[0:1:0],
%sea_surface_temperature[0:1:0][0:1:1326][0:1:3399],
%sst_dtime[0:1:0][0:1:1326][0:1:3399],
%SSES_bias_error[0:1:0][0:1:1326][0:1:3399],
%SSES_standard_deviation_error[0:1:0][0:1:1326][0:1:3399],
%DT_analysis[0:1:0][0:1:1326][0:1:3399],
%ssi_dtime_from_sst[0:1:0][0:1:1326][0:1:3399],
%wind_speed[0:1:0][0:1:1326][0:1:3399],
%wind_speed_dtime_from_sst[0:1:0][0:1:1326][0:1:3399],
%sea_ice_fraction[0:1:0][0:1:1326][0:1:3399],
%aerosol_optical_depth[0:1:0][0:1:1326][0:1:3399],
%aod_dtime_from_sst[0:1:0][0:1:1326][0:1:3399],
%sources_of_wind_speed[0:1:0][0:1:1326][0:1:3399],
%sources_of_ssi[0:1:0][0:1:1326][0:1:3399],
%sources_of_sea_ice_fraction[0:1:0][0:1:1326][0:1:3399],
%sources_of_aod[0:1:0][0:1:1326][0:1:3399],
%satellite_zenith_angle[0:1:0][0:1:1326][0:1:3399],
%rejection_flag[0:1:0][0:1:1326][0:1:3399],
%confidence_flag[0:1:0][0:1:1326][0:1:3399],
%proximity_confidence[0:1:0][0:1:1326][0:1:3399],
%probability_of_clear_sky[0:1:0][0:1:1326][0:1:3399]

DSI.URLSITE='';
DSI.URLPATH='';

function DSI=DSI_GHRSST_GOES12
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_GOES12';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P/GOES12/OSDPD/
%2008/045/20080214-GOES12-OSDPD-L2P-GOES12_North_0715Z-v01.nc.gz?
%lat[0:1:1826][0:1:3459],lon[0:1:1826][0:1:3459],time[0:1:0],
%sea_surface_temperature[0:1:0][0:1:1826][0:1:3459],
%sst_dtime[0:1:0][0:1:1826][0:1:3459],
%SSES_bias_error[0:1:0][0:1:1826][0:1:3459],
%SSES_standard_deviation_error[0:1:0][0:1:1826][0:1:3459],
%DT_analysis[0:1:0][0:1:1826][0:1:3459],
%surface_solar_irradiance[0:1:0][0:1:1826][0:1:3459],
%ssi_dtime_from_sst[0:1:0][0:1:1826][0:1:3459],
%wind_speed[0:1:0][0:1:1826][0:1:3459],
%wind_speed_dtime_from_sst[0:1:0][0:1:1826][0:1:3459],
%sea_ice_fraction[0:1:0][0:1:1826][0:1:3459],
%aerosol_optical_depth[0:1:0][0:1:1826][0:1:3459],
%aod_dtime_from_sst[0:1:0][0:1:1826][0:1:3459],
%sources_of_wind_speed[0:1:0][0:1:1826][0:1:3459],
%sources_of_ssi[0:1:0][0:1:1826][0:1:3459],
%sources_of_sea_ice_fraction[0:1:0][0:1:1826][0:1:3459],
%sources_of_aod[0:1:0][0:1:1826][0:1:3459],
%satellite_zenith_angle[0:1:0][0:1:1826][0:1:3459],
%rejection_flag[0:1:0][0:1:1826][0:1:3459],
%confidence_flag[0:1:0][0:1:1826][0:1:3459],
%proximity_confidence[0:1:0][0:1:1826][0:1:3459],
%probability_of_clear_sky[0:1:0][0:1:1826][0:1:3459]

DSI.URLSITE='';
DSI.URLPATH='';

%L2P_GRIDDED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DSI=DSI_GHRSST_AMSRE
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_AMSRE';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/AMSRE/REMSS/
%2002/152/20020601-AMSRE-REMSS-L2P_GRIDDED_25-amsre_20020601v5-v01.nc.gz?
%lat[0:1:719],lon[0:1:1439],time[0:1:1],
%sea_surface_temperature[0:1:1][0:1:719][0:1:1439],
%sst_dtime[0:1:1][0:1:719][0:1:1439],
%SSES_bias_error[0:1:1][0:1:719][0:1:1439],
%SSES_standard_deviation_error[0:1:1][0:1:719][0:1:1439],
%wind_speed[0:1:1][0:1:719][0:1:1439],
%rejection_flag[0:1:1][0:1:719][0:1:1439],
%confidence_flag[0:1:1][0:1:719][0:1:1439],
%proximity_confidence[0:1:1][0:1:719][0:1:1439],
%diurnal_amplitude[0:1:1][0:1:719][0:1:1439],
%cool_skin[0:1:1][0:1:719][0:1:1439]

% define list of fields or variables available in the dataset
LaLo={'latitude','longitude'};
TLaLo={'time','Latitude','Longitude'};
LaLoT={'Latitude','Longitude','time'};

DSI.FieldsTable={
%    opendap    file name        name in                         returned name                                                   menu name                           opendap dataset returned 
%    name                        dataset                                                                                                                              coord  coord   coord   CBNum
    'sst'            ,''        ,'sea_surface_temperature'      ,'sea_surface_temperature.sea_surface_temperature'             ,'sea surface temperature'            ,LaLo   ,TLaLo  ,LaLoT   ,1 ,'Kelvin';
    'sst_dtime'      ,''        ,'sst_dtime'                    ,'sst_dtime.sst_dtime'                                         ,'time difference from reference time',LaLo   ,TLaLo  ,LaLoT   ,2 ,'second';
    'SSES_bias_error',''        ,'SSES_bias_error'              ,'SSES_bias_error.SSES_bias_error'                             ,'SSES bias error'                    ,LaLo   ,TLaLo  ,LaLoT   ,3 ,'Kelvin';
    'SSES_standard_deviation','','SSES_standard_deviation_error','SSES_standard_deviation_error.SSES_standard_deviation_error' ,'SSES standard deviation error'      ,LaLo   ,TLaLo  ,LaLoT   ,4 ,'Kelvin';
    'wind_speed'     ,''        ,'wind_speed'                   ,'wind_speed.wind_speed'                                       ,'wind speed'                         ,LaLo   ,TLaLo  ,LaLoT   ,5 ,'m s-1';
    'rejection_flag' ,''        ,'rejection_flag'               ,'rejection_flag.rejection_flag'                               ,'rejection flag'                     ,LaLo   ,TLaLo  ,LaLoT   ,6 ,'';
    'confidence_flag',''        ,'confidence_flag'              ,'confidence_flag.confidence_flag'                             ,'confidence flag'                    ,LaLo   ,TLaLo  ,LaLoT   ,7 ,'';
    'proximity_confidence',''   ,'proximity_confidence'         ,'proximity_confidence.proximity_confidence'                   ,'proximity confidence'               ,LaLo   ,TLaLo  ,LaLoT   ,8 ,'';         
    'diurnal_amplitude'   ,''   ,'diurnal_amplitude'            ,'diurnal_amplitude.diurnal_amplitude'                         ,'diurnal amplitude'                  ,LaLo   ,TLaLo  ,LaLoT   ,9 ,'Kelvin';
    'cool_skin'           ,''   ,'cool_skin'                    ,'cool_skin.cool_skin'                                         ,'cool skin'                          ,LaLo   ,TLaLo  ,LaLoT   ,10,'Kelvin';
     };
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/AMSRE/REMSS/
%2007/025/20070125-AMSRE-REMSS-L2P_GRIDDED_25-amsre_20070125rt-v01.nc.gz?
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/AMSRE/REMSS/
%2007/025/20070125-AMSRE-REMSS-L2P_GRIDDED_25-amsre_20070125v5-v01.nc.gz?
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/AMSRE/REMSS/
%2008/010/20080110-AMSRE-REMSS-L2P_GRIDDED_25-amsre_20080110rt-v01.nc.gz
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/AMSRE/REMSS/
%2008/010/20080110-AMSRE-REMSS-L2P_GRIDDED_25-amsre_20080110v5-v01.nc.gz

DSI.Range.Time1_NODC=datenum(2002,01,152);
%DSI.Range.Time2_NODC=datenum(2009,01,239); % update
DSI.Range.Time2_JPL=now;
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L2P_GRIDDED/AMSRE/REMSS';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/AMSRE/REMSS';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';
%DSI.URLPATH='/2002/152/';
%DSI.URLFILE='20020601-AMSRE-REMSS-L2P_GRIDDED_25-amsre_20020601v5-v01.nc.gz';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-AMSRE-REMSS-L2P_GRIDDED_25-amsre_'' YYYYMMDD ''v5-v01.nc.gz''];';

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI.TimeofFrame=[0 0.5];
DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    

%                  long_name: '"sea surface temperature"'
%                      units: '"kelvin"'
%              ml__FillValue: -32768
%                 add_offset: 273.1500
%               scale_factor: 0.0100
%                  valid_min: -5000
%                  valid_max: 5000

  % long_name: '"time difference from reference time"'
  %              units: '"second"'
  %      ml__FillValue: -32768
  %         add_offset: 0
  %       scale_factor: 10
  %  DODS_ML_Real_Name: 'sst_dtime'
 
  % long_name: '"SSES bias error"'
  %              units: '"kelvin"'
  %      ml__FillValue: 128
  %         add_offset: 0
  %       scale_factor: 0.0100
  %          valid_min: 129
  %          valid_max: 127
  % DODS_ML_Real_Name: 'SSES_bias_error'
  
  %        long_name: '"SSES standard deviation "'
  %                          units: '"kelvin"'
  %                  ml__FillValue: 128
  %                     add_offset: 0.7500
  %                   scale_factor: 0.0100
  %                      valid_min: 129
  %                      valid_max: 127
  %              DODS_ML_Real_Name: 'SSES_standard_deviation_error'
  
  %    long_name: '"wind speed"'
  %                      units: '"m s-1"'
  %              ml__FillValue: 128
  %                 add_offset: 25
  %               scale_factor: 0.2000
  %                  valid_min: 129
  %                  valid_max: 127
  %                     source: '"native_AMSRE_wind"'
  %  dtime_from_sst_in_seconds: '"0"'
  %          DODS_ML_Real_Name: 'wind_speed'
            
% 
%            long_name: '"rejection flag"'
%              comment: [1x87 char]
%    DODS_ML_Real_Name: 'rejection_flag'
%       rejection_flag: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]
%A.rejection_flag.comment
%ans =
%"b0:1=rain;b1:1=high wind/sunglint;b2:1=ice;b3:1=near land;b4:1=no data;b5:1=land;    "
DSI.FormulaNaN.rejection_flag='';
DSI.FormulaCnv.rejection_flag='';

%  A.confidence_flag
%ans = 
%            long_name: '"confidence flag"'
%              comment: [1x310 char]
%    DODS_ML_Real_Name: 'confidence_flag'
%      confidence_flag: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]
%A.confidence_flag.comment
%ans =
%"b0:1=within 50km rain, 0.6 dif mwoisst;b1:1=within 100km rain, 0.8 dif mwoisst;b2:1=within 150km ice, 0.6 dif mwoisst;b3:1=more than 5deg dif mwoisst;b4:1=3-sigma test;b5:1=(tmi only) within 150 km of land and 0.6 warmer than mwoisst;b6:1=diurnal estimate > 0.6 warming;b7:1=diurnal estimate > 0.3 warming   "
DSI.FormulaNaN.confidence_flag='';
DSI.FormulaCnv.confidence_flag='';

%A.proximity_confidence
%ans = 
%               long_name: '"proximity confidence "'
%                 comment: [1x164 char]
%       DODS_ML_Real_Name: 'proximity_confidence'
%    proximity_confidence: [1x1 struct]
%                    time: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
%A.proximity_confidence.comment
%ans =
%"1=Bad, data rejected;2=Suspected Bad, data that has any confidence flags bit 0-5 thrown;3=Unprocessed proximity confidence flag, should be Good data;4=Good data  "
DSI.FormulaNaN.proximity_confidence='';
DSI.FormulaCnv.proximity_confidence='';

  %          long_name: '"Diurnal warming amplitude"'
  %              units: '"kelvin"'
  %      ml__FillValue: 128
  %         add_offset: 1
  %       scale_factor: 0.0200
  %          valid_min: 129
  %          valid_max: 127
  %            comment: '"non L2P_GRIDDED core field"'
  %  DODS_ML_Real_Name: 'diurnal_amplitude'
  
   % long_name: '"cool skin"'
   %             units: '"kelvin"'
   %     ml__FillValue: 128
   %        add_offset: -1
   %      scale_factor: 0.0100
   %         valid_min: 129
   %         valid_max: 127
   %           comment: '"non L2P_GRIDDED core field"'
   % DODS_ML_Real_Name: 'cool_skin'
  
DSI=DSI_common3(DSI);    
SaveDSI(DSI)


function DSI=DSI_GHRSST_TMI
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_TMI';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/TMI/REMSS/
%1998/001/19980101-TMI-REMSS-L2P_GRIDDED_25-tmi_19980101v4-v01.nc.gz?
%lat[0:1:319],lon[0:1:1439],time[0:1:1],
%sea_surface_temperature[0:1:1][0:1:319][0:1:1439],
%sst_dtime[0:1:1][0:1:319][0:1:1439],
%SSES_bias_error[0:1:1][0:1:319][0:1:1439],
%SSES_standard_deviation_error[0:1:1][0:1:319][0:1:1439],
%wind_speed[0:1:1][0:1:319][0:1:1439],
%rejection_flag[0:1:1][0:1:319][0:1:1439],
%confidence_flag[0:1:1][0:1:319][0:1:1439],
%proximity_confidence[0:1:1][0:1:319][0:1:1439],
%diurnal_amplitude[0:1:1][0:1:319][0:1:1439],
%cool_skin[0:1:1][0:1:319][0:1:1439]

% define list of fields or variables available in the dataset
LaLo={'latitude','longitude'};
TLaLo={'time','Latitude','Longitude'};
LaLoT={'Latitude','Longitude','time'};

DSI.FieldsTable={
%    opendap    file name        name in                         returned name                                                   menu name                           opendap dataset returned 
%    name                        dataset                                                                                                                              coord  coord   coord   CBNum
    'sst'      ,''              ,'sea_surface_temperature'      ,'sea_surface_temperature.sea_surface_temperature'             ,'sea surface temperature'            ,LaLo   ,TLaLo  ,LaLoT   ,1 ,'Kelvin';
    'sst_dtime',''              ,'sst_dtime'                    ,'sst_dtime.sst_dtime'                                         ,'time difference from reference time',LaLo   ,TLaLo  ,LaLoT   ,2 ,'second';
    'SSES_bias_error',''        ,'SSES_bias_error'              ,'SSES_bias_error.SSES_bias_error'                             ,'SSES bias error'                    ,LaLo   ,TLaLo  ,LaLoT   ,3 ,'Kelvin';
    'SSES_standard_deviation','','SSES_standard_deviation_error','SSES_standard_deviation_error.SSES_standard_deviation_error' ,'SSES standard deviation error'      ,LaLo   ,TLaLo  ,LaLoT   ,4 ,'Kelvin';
    'wind_speed'  ,''           ,'wind_speed'                   ,'wind_speed.wind_speed'                                       ,'wind speed'                         ,LaLo   ,TLaLo  ,LaLoT   ,5 ,'m s-1';
    'rejection_flag',''         ,'rejection_flag'               ,'rejection_flag.rejection_flag'                               ,'rejection flag'                     ,LaLo   ,TLaLo  ,LaLoT   ,6 ,'';
    'confidence_flag',''        ,'confidence_flag'              ,'confidence_flag.confidence_flag'                             ,'confidence flag'                    ,LaLo   ,TLaLo  ,LaLoT   ,7 ,'';
    'proximity_confidence',''   ,'proximity_confidence'         ,'proximity_confidence.proximity_confidence'                   ,'proximity confidence'               ,LaLo   ,TLaLo  ,LaLoT   ,8 ,'';         
    'diurnal_amplitude'   ,''   ,'diurnal_amplitude'            ,'diurnal_amplitude.diurnal_amplitude'                         ,'diurnal amplitude'                  ,LaLo   ,TLaLo  ,LaLoT   ,9 ,'Kelvin';
    'cool_skin'           ,''   ,'cool_skin'                    ,'cool_skin.cool_skin'                                         ,'cool skin'                          ,LaLo   ,TLaLo  ,LaLoT   ,10,'Kelvin';
     };

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/TMI/REMSS/
%1998/001/19980101-TMI-REMSS-L2P_GRIDDED_25-tmi_19980101v4-v01.nc.gz?
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/TMI/REMSS/
%2006/001/20060101-TMI-REMSS-L2P_GRIDDED_25-20060101tm.dat-v01.nc.gz
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/TMI/REMSS/
%2007/022/20070122-TMI-REMSS-L2P_GRIDDED_25-tmi_20070122rt-v01.nc.gz

DSI.Range.Time1_NODC=datenum(1998,01,01);
DSI.Range.Time2_NODC=datenum(2009,01,238); % update
DSI.Range.Time2_JPL=now;
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L2P_GRIDDED/TMI/REMSS';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P_GRIDDED/TMI/REMSS';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';
%DSI.URLPATH='/1998/001/';
%DSI.URLFILE='19980101-TMI-REMSS-L2P_GRIDDED_25-tmi_19980101v4-v01.nc.gz';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-TMI-REMSS-L2P_GRIDDED_25-tmi_'' YYYYMMDD ''v4-v01.nc.gz''];';
DSI.TimeofFrame=[0 0.5];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty 
% fill values are the same as in AMSRE 
DSI.FormulaNaN.rejection_flag='';
DSI.FormulaCnv.rejection_flag='';
DSI.FormulaNaN.confidence_flag='';
DSI.FormulaCnv.confidence_flag='';
DSI.FormulaNaN.proximity_confidence='';
DSI.FormulaCnv.proximity_confidence='';

DSI=DSI_common3(DSI);
SaveDSI(DSI)

%L4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DSI=DSI_GHRSST_AUS
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_AUS';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/AUS/ABOM/RAMSSA_09km/
%2008/092/20080401-ABOM-L4HRfnd-AUS-v01-fv01_0-RAMSSA_09km.nc.bz2?
%lon[0:1:1560],lat[0:1:1080],time[0:1:0],
%analysed_sst[0:1:0][0:1:1080][0:1:1560],
%analysis_error[0:1:0][0:1:1080][0:1:1560],
%sea_ice_fraction[0:1:0][0:1:1080][0:1:1560],
%mask[0:1:0][0:1:1080][0:1:1560]

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/AUS/ABOM/RAMSSA_09km/
%2008/092/20080505-ABOM-L4HRfnd-AUS-v01-fv02_0-RAMSSA_09km.nc.bz2?

% define list of fields or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name          opendap dataset returned 
%    name                              dataset                                                                     coord   coord   coord    CBNum
    'analysed_sst'      ,''          ,'analysed_sst'      ,'analysed_sst.analysed_sst'        ,'analysed sst'      ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'    ,''          ,'analysis_error'    ,'analysis_error.analysis_error'    ,'analysis error'    ,LaLo   ,LaLo   ,LaLo   ,12,'Kelvin';      
    'sea_ice_fraction'  ,''          ,'sea_ice_fraction'  ,'sea_ice_fraction.sea_ice_fraction','sea ice fraction'  ,LaLo   ,LaLo   ,LaLo   ,13,'';    
    'mask'              ,''          ,'mask'              ,'mask.mask'                        ,'mask'              ,LaLo   ,LaLo   ,LaLo   ,14,'';
    };

DSI.Range.Time1_NODC=datenum(2008,01,092);
DSI.Range.Time2_NODC=datenum(2009,01,239); % update
DSI.Range.Time2_JPL=now;
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/AUS/ABOM/RAMSSA_09km';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/AUS/ABOM/RAMSSA_09km';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2008/126/';
%DSI.URLFILE='20080505-ABOM-L4HRfnd-AUS-v01-fv01_0-RAMSSA_09km.nc.bz2';
%http://data.nodc.noaa.gov/opendap/ghrsst/L4/AUS/ABOM/RAMSSA_09km/2008/092/20080401-ABOM-L4HRfnd-AUS-v01-fv02_0-RAMSSA_09km.nc.bz2
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-ABOM-L4HRfnd-AUS-v01-'' ''fv02_0-RAMSSA_09km.nc.bz2''];';
DSI.TimeofFrame=[0.5];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    

 %  long_name: '"analysed sea surface temperature"'
 %       standard_name: '"sea_surface_temperature"'
 %                type: '"foundation"'
 %               units: '"kelvin"'
 %          add_offset: 273.1500
 %        scale_factor: 0.0100
 %       ml__FillValue: -32768
 %           valid_min: -300
 %           valid_max: 4500
 %   DODS_ML_Real_Name: 'analysed_sst'
 
  %          long_name: '"estimated error standard deviation of analysed_sst"'
  %              units: '"kelvin"'
  %         add_offset: 0
  %       scale_factor: 0.0100
  %      ml__FillValue: -32768
  %          valid_min: 0
  %          valid_max: 32767
  %  DODS_ML_Real_Name: 'analysis_error'
 
    %        long_name: '"sea ice area fraction"'
    %    standard_name: '"sea_ice_area_fraction"'
    %            units: '"none (fraction)"'
    %       add_offset: 0
    %     scale_factor: 0.0100
    %    ml__FillValue: -128
    %        valid_min: 0
    %        valid_max: 100
    %           source: '"NCEP_ICE"'
    %DODS_ML_Real_Name: 'sea_ice_fraction'

     %       long_name: '"sea/land/lake/ice field composite mask"'
     %   ml__FillValue: -128
     %         comment: [1x213 char]
    %DODS_ML_Real_Name: 'mask'
    %             mask: [1x1 struct]
    %             time: [1x1 struct]
    %              lat: [1x1 struct]
    %              lon: [1x1 struct]
%a.mask.comment
%ans =
%"b0: 1=grid cell is open sea water. b1: 1=land is present in this grid cell. b2: 1=lake surface is present in this grid cell. b3: 1=sea ice is present in this grid cell. b4-b7: reserved for future grid mask data."
DSI.FormulaCnv.mask='';
DSI=DSI_common3(DSI);
SaveDSI(DSI)


function DSI=DSI_GHRSST_GAL
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_GAL';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GAL/EUR/ODYSSEA/
%2008/023/20080123-EUR-L4UHRfnd-GAL-v01-fv01-ODYSSEA.nc.bz2?
%time[0:1:0],lat[0:1:1999],lon[0:1:2549],
%analysed_sst[0:1:0][0:1:1999][0:1:2549],
%analysis_error[0:1:0][0:1:1999][0:1:2549],
%mask[0:1:0][0:1:1999][0:1:2549],
%sea_ice_fraction[0:1:0][0:1:1999][0:1:2549]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name          opendap dataset returned 
%    name                              dataset                                                                     coord   coord   coord    CBNum
    'analysed_sst'      ,''          ,'analysed_sst'      ,'analysed_sst.analysed_sst'        ,'analysed sst'      ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'    ,''          ,'analysis_error'    ,'analysis_error.analysis_error'    ,'analysis error'    ,LaLo   ,LaLo   ,LaLo   ,12,'percent';      
    'sea_ice_fraction'  ,''          ,'sea_ice_fraction'  ,'sea_ice_fraction.sea_ice_fraction','sea ice fraction'  ,LaLo   ,LaLo   ,LaLo   ,13,'';    
    'mask'              ,''          ,'mask'              ,'mask.mask'                        ,'mask'              ,LaLo   ,LaLo   ,LaLo   ,14,'';
    };

DSI.Range.Time1_NODC=datenum(2008,01,023);
DSI.Range.Time2_NODC=datenum(2008,01,350); % update
DSI.Range.Time2_JPL=DSI.Range.Time2_NODC; %now; JPL has same data
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/GAL/EUR/ODYSSEA';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GAL/EUR/ODYSSEA';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2008/023/';
%DSI.URLFILE='20080123-EUR-L4UHRfnd-GAL-v01-fv01-ODYSSEA.nc.bz2';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-EUR-L4UHRfnd-GAL-v01-'' ''fv01-ODYSSEA.nc.bz2''];';
DSI.TimeofFrame=[0];

YYYY='2008' %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    

  %    long_name: '"analysed sea surface temperature"'
  %      standard_name: '"sea_surface_temperature"'
  %               type: '"foundation"'
  %              units: '"kelvin"'
  %      ml__FillValue: -32768
  %         add_offset: 273.1500
  %       scale_factor: 0.0100
  %          valid_min: -300
  %          valid_max: 4500
  %  DODS_ML_Real_Name: 'analysed_sst'
    
   %         long_name: '"estimated error percentage"'
   %             units: '"percent"'
   %     ml__FillValue: -32768
   %        add_offset: 0
   %      scale_factor: 0.0100
   %         valid_min: 0
   %         valid_max: 10000
   % DODS_ML_Real_Name: 'analysis_error'
   
    %        long_name: '"sea/land/lake/ice field composite mask"'
    %      flag_values: '"b0, 1=sea, b1: 1=land, b2: 1=lake, b3: 1=ice"'
    %    flag_meanings: '"sea land lake ice"'
    %DODS_ML_Real_Name: 'mask'
DSI.FormulaNaN.mask='';
DSI.FormulaCnv.mask='';

   
    %        long_name: '"sea ice area fraction"'
    %            units: '" "'
    %    ml__FillValue: 128
    %       add_offset: 0
    %     scale_factor: 0.0100
    %        valid_min: 0
    %        valid_max: 100
    %      institution: '"SAF O&SI"'
    %DODS_ML_Real_Name: 'sea_ice_fraction'
   
DSI=DSI_common3(DSI);
SaveDSI(DSI)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DSI=DSI_GHRSST_GLOB_EUR
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_GLOB_EUR';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/EUR/ODYSSEA/
%2007/274/20071001-EUR-L4HRfnd-GLOB-v01-fv01-ODYSSEA.nc.bz2?
%time[0:1:0],lat[0:1:1599],lon[0:1:3599],
%analysed_sst[0:1:0][0:1:1599][0:1:3599],
%analysis_error[0:1:0][0:1:1599][0:1:3599],
%mask[0:1:0][0:1:1599][0:1:3599],
%sea_ice_fraction[0:1:0][0:1:1599][0:1:3599]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name          opendap dataset returned 
%    name                              dataset                                                                     coord   coord   coord    CBNum
    'analysed_sst'      ,''          ,'analysed_sst'      ,'analysed_sst.analysed_sst'        ,'analysed sst'      ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'    ,''          ,'analysis_error'    ,'analysis_error.analysis_error'    ,'analysis error'    ,LaLo   ,LaLo   ,LaLo   ,12,'percent';      
    'sea_ice_fraction'  ,''          ,'sea_ice_fraction'  ,'sea_ice_fraction.sea_ice_fraction','sea ice fraction'  ,LaLo   ,LaLo   ,LaLo   ,13,'';    
    'mask'              ,''          ,'mask'              ,'mask.mask'                        ,'mask'              ,LaLo   ,LaLo   ,LaLo   ,14,'';
    };

DSI.Range.Time1_NODC=datenum(2007,01,274);
DSI.Range.Time2_NODC=datenum(2009,01,239); % update
DSI.Range.Time2_JPL=now; 
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/GLOB/EUR/ODYSSEA';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/EUR/ODYSSEA';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2007/274/';
%DSI.URLFILE='20071001-EUR-L4HRfnd-GLOB-v01-fv01-ODYSSEA.nc.bz2';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-EUR-L4HRfnd-GLOB-v01-'' ''fv01-ODYSSEA.nc.bz2''];';
DSI.TimeofFrame=[0];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    

 
  %          long_name: '"analysed sea surface temperature"'
  %      standard_name: '"sea_surface_temperature"'
  %               type: '"foundation"'
  %              units: '"kelvin"'
  %      ml__FillValue: -32768
  %         add_offset: 273.1500
  %       scale_factor: 0.0100
  %          valid_min: -300
  %          valid_max: 4500
  %  DODS_ML_Real_Name: 'analysed_sst'

%            long_name: '"estimated error percentage"'
%                units: '"percent"'
%        ml__FillValue: -32768
%           add_offset: 0
%         scale_factor: 0.0100
%            valid_min: 0
%            valid_max: 10000
%    DODS_ML_Real_Name: 'analysis_error'

%            long_name: '"sea/land/lake/ice field composite mask"'
%          flag_values: '"b0, 1=sea, b1: 1=land, b2: 1=lake, b3: 1=ice"'
%        flag_meanings: '"sea land lake ice"'
%    DODS_ML_Real_Name: 'mask'
DSI.FormulaNaN.mask='';
DSI.FormulaCnv.mask='';

%          long_name: '"sea ice area fraction"'
%                units: '" "'
%        ml__FillValue: 128
%           add_offset: 0
%         scale_factor: 0.0100
%            valid_min: 0
%            valid_max: 100
%          institution: '"SAF O&SI"'
%    DODS_ML_Real_Name: 'sea_ice_fraction'

DSI=DSI_common3(DSI);
SaveDSI(DSI)

function DSI=DSI_GHRSST_GLOB_NAVO
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_GLOB_NAVO';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/NAVO/K10_SST/
%2008/092/20080401-NAVO-L4HR1m-GLOB-v01-fv01_0-K10_SST.nc.bz2?
%lat[0:1:1800],lon[0:1:3599],time[0:1:0],
%analysed_sst[0:1:0][0:1:1800][0:1:3599],
%analysis_error[0:1:0][0:1:1800][0:1:3599],
%mask[0:1:0][0:1:1800][0:1:3599]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name          opendap dataset returned 
%    name                              dataset                                                                     coord   coord   coord    CBNum
    'analysed_sst'      ,''          ,'analysed_sst'      ,'analysed_sst.analysed_sst'        ,'analysed sst'      ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'    ,''          ,'analysis_error'    ,'analysis_error.analysis_error'    ,'analysis error'    ,LaLo   ,LaLo   ,LaLo   ,12,'Kelvin';      
    'mask'              ,''          ,'mask'              ,'mask.mask'                        ,'mask'              ,LaLo   ,LaLo   ,LaLo   ,14,'';
    };


DSI.Range.Time1_NODC=datenum(2008,01,092);
DSI.Range.Time2_NODC=datenum(2009,01,238); % update
DSI.Range.Time2_JPL=now; 
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/GLOB/NAVO/K10_SST';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/NAVO/K10_SST';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2008/092/';
%DSI.URLFILE='20080401-NAVO-L4HR1m-GLOB-v01-fv01_0-K10_SST.nc.bz2';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-NAVO-L4HR1m-GLOB-v01-'' ''fv01_0-K10_SST.nc.bz2''];';
DSI.TimeofFrame=[0];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    
 
%         long_name: '"Analyzed Sea Surface Temperature"'
%        standard_name: '"sea_surface_temperature"'
%                 type: '"depth 1m"'
%                units: '"kelvin"'
%        ml__FillValue: -32768
%           add_offset: 273.1500
%         scale_factor: 0.1000
%            valid_min: -20
%            valid_max: 350
%    DODS_ML_Real_Name: 'analysed_sst'
  
 %           long_name: '"Analysis Error"'
 %               units: '"kelvin"'
 %       ml__FillValue: -32768
 %        scale_factor: 0.0100
 %           valid_min: 0.4000
 %           valid_max: 1.5000
 %   DODS_ML_Real_Name: 'analysis_error'
 FIELD='analysis_error';
 FIELDDS='analysis_error';
 eval(['DSI.FormulaCnv.' FIELD '=''f=A.' FIELDDS '.scale_factor*f;'';']);

  %          long_name: '"Land Sea Mask"'
  %        flag_values: [2x1 double]
  %      flag_meanings: '"water, land"'
  %            comment: '"bit 0 = 1 if grid cell contains water; "'
  %  DODS_ML_Real_Name: 'mask'
DSI.FormulaNaN.mask='';
DSI.FormulaCnv.mask='';

DSI=DSI_common3(DSI);
SaveDSI(DSI)

function DSI=DSI_GHRSST_GLOB_NCDC_AVHRR_AMSR_OI
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_GLOB_NCDC_AVHRR_AMSR_OI';
 
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/NCDC/AVHRR_AMSR_OI/
%2002/152/20020601-NCDC-L4LRblend-GLOB-v01-fv01_0-AVHRR_AMSR_OI.nc.bz2?
%lat[0:1:719],lon[0:1:1439],time[0:1:0],
%analysed_sst[0:1:0][0:1:719][0:1:1439],
%analysis_error[0:1:0][0:1:719][0:1:1439],
%mask[0:1:0][0:1:719][0:1:1439],
%sea_ice_fraction[0:1:0][0:1:719][0:1:1439]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name          opendap dataset returned 
%    name                              dataset                                                                     coord   coord   coord    CBNum
    'analysed_sst'      ,''          ,'analysed_sst'      ,'analysed_sst.analysed_sst'        ,'analysed sst'      ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'    ,''          ,'analysis_error'    ,'analysis_error.analysis_error'    ,'analysis error'    ,LaLo   ,LaLo   ,LaLo   ,12,'Kelvin';      
    'sea_ice_fraction'  ,''          ,'sea_ice_fraction'  ,'sea_ice_fraction.sea_ice_fraction','sea ice fraction'  ,LaLo   ,LaLo   ,LaLo   ,13,'percent';    
    'mask'              ,''          ,'mask'              ,'mask.mask'                        ,'mask'              ,LaLo   ,LaLo   ,LaLo   ,14,'';
    };

DSI.Range.Time1_NODC=datenum(2002,01,152);
DSI.Range.Time2_NODC=datenum(2009,01,238); % update
DSI.Range.Time2_JPL=now; 
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/GLOB/NCDC/AVHRR_AMSR_OI';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/NCDC/AVHRR_AMSR_OI';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2002/152/';
%DSI.URLFILE='20020601-NCDC-L4LRblend-GLOB-v01-fv01_0-AVHRR_AMSR_OI.nc.bz2';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-NCDC-L4LRblend-GLOB-v01-'' ''fv01_0-AVHRR_AMSR_OI.nc.bz2''];';
DSI.TimeofFrame=[0];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    
 
%        long_name: '"analysed sea surface temperature"'
%        standard_name: '"sea_surface_temperature"'
%                 type: '"foundation"'
%                units: '"kelvin"'
%        ml__FillValue: -32768
%           add_offset: 273.1500
%         scale_factor: 0.0100
%            valid_min: -300
%            valid_max: 4500
%    DODS_ML_Real_Name: 'analysed_sst'
 
 %           long_name: '"estimated error standard deviation of analysed_sst"'
 %               units: '"kelvin"'
 %       ml__FillValue: -32768
 %          add_offset: 0
 %        scale_factor: 0.0100
 %           valid_min: 0
 %           valid_max: 127
 %   DODS_ML_Real_Name: 'analysis_error'
    
  %          long_name: '"sea/land field composite mask"'
  %      ml__FillValue: 128
  %        flag_values: 1
  %      flag_meanings: '"sea land lake ice"'
  %            comment: [1x202 char]
  %  DODS_ML_Real_Name: 'mask'
DSI.FormulaCnv.mask='';

   %         long_name: '"sea ice area fraction"'
   %     standard_name: '"sea ice area fraction"'
   %             units: '"percent"'
   %     ml__FillValue: 128
   %        add_offset: 0
   %      scale_factor: 0.0100
   %         valid_min: 0
   %         valid_max: 100
   % DODS_ML_Real_Name: 'sea_ice_fraction'
DSI=DSI_common3(DSI);
SaveDSI(DSI)

function DSI=DSI_GHRSST_GLOB_NCDC_AVHRR_OI
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_GLOB_NCDC_AVHRR_OI';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/NCDC/AVHRR_OI/
%1985/004/19850104-NCDC-L4LRblend-GLOB-v01-fv01_0-AVHRR_OI.nc.bz2?
%lat[0:1:719],lon[0:1:1439],time[0:1:0],
%analysed_sst[0:1:0][0:1:719][0:1:1439],
%analysis_error[0:1:0][0:1:719][0:1:1439],
%mask[0:1:0][0:1:719][0:1:1439],
%sea_ice_fraction[0:1:0][0:1:719][0:1:1439]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name          opendap dataset returned 
%    name                              dataset                                                                     coord   coord   coord    CBNum
    'analysed_sst'      ,''          ,'analysed_sst'      ,'analysed_sst.analysed_sst'        ,'analysed sst'      ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'    ,''          ,'analysis_error'    ,'analysis_error.analysis_error'    ,'analysis error'    ,LaLo   ,LaLo   ,LaLo   ,12,'Kelvin';      
    'sea_ice_fraction'  ,''          ,'sea_ice_fraction'  ,'sea_ice_fraction.sea_ice_fraction','sea ice fraction'  ,LaLo   ,LaLo   ,LaLo   ,13,'percent';    
    'mask'              ,''          ,'mask'              ,'mask.mask'                        ,'mask'              ,LaLo   ,LaLo   ,LaLo   ,14,'';
    };

DSI.Range.Time1_NODC=datenum(1981,01,244);
DSI.Range.Time2_NODC=datenum(2009,01,238); % update
DSI.Range.Time2_JPL=now; 
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/GLOB/NCDC/AVHRR_OI';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/NCDC/AVHRR_OI';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/1985/004/';
%DSI.URLFILE='19850104-NCDC-L4LRblend-GLOB-v01-fv01_0-AVHRR_OI.nc.bz2';
%             19820102-NCDC-L4LRblend-GLOB-v01-fv02_0-AVHRR_OI.nc.bz2
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-NCDC-L4LRblend-GLOB-v01-'' ''fv02_0-AVHRR_OI.nc.bz2''];';
DSI.TimeofFrame=[0];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    
 
 %           long_name: '"analysed sea surface temperature"'
 %       standard_name: '"sea_surface_temperature"'
 %                type: '"foundation"'
 %               units: '"kelvin"'
 %       ml__FillValue: -32768
 %          add_offset: 273.1500
 %        scale_factor: 0.0100
 %           valid_min: -300
 %           valid_max: 4500
 %   DODS_ML_Real_Name: 'analysed_sst'

  %          long_name: '"estimated error standard deviation of analysed_sst"'
  %              units: '"kelvin"'
  %      ml__FillValue: -32768
  %         add_offset: 0
  %       scale_factor: 0.0100
  %          valid_min: 0
  %          valid_max: 127
  %  DODS_ML_Real_Name: 'analysis_error'

 %           long_name: '"sea/land field composite mask"'
 %       ml__FillValue: 128
 %         flag_values: 1
 %       flag_meanings: '"sea land lake ice"'
 %             comment: [1x202 char]
 %   DODS_ML_Real_Name: 'mask'
DSI.FormulaCnv.mask='';

  %  long_name: '"sea ice area fraction"'
  %      standard_name: '"sea ice area fraction"'
  %              units: '"percent"'
  %      ml__FillValue: 128
  %         add_offset: 0
  %       scale_factor: 0.0100
  %          valid_min: 0
  %          valid_max: 100
  %  DODS_ML_Real_Name: 'sea_ice_fraction'
DSI=DSI_common3(DSI);
SaveDSI(DSI)

function DSI=DSI_GHRSST_GLOB_REMSS
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_GLOB_REMSS';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/REMSS/mw_ir_OI/
%2005/233/20050821-REMSS-L4HRfnd-GLOB-v01-fv01-mw_ir_OI.nc.gz?
%lat[0:1:2047],lon[0:1:4095],time[0:1:0],
%analysed_sst[0:1:0][0:1:2047][0:1:4095],
%analysis_error[0:1:0][0:1:2047][0:1:4095],
%mask[0:1:0][0:1:2047][0:1:4095]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name          opendap dataset returned 
%    name                              dataset                                                                     coord   coord   coord    CBNum
    'analysed_sst'      ,''          ,'analysed_sst'      ,'analysed_sst.analysed_sst'        ,'analysed sst'      ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'    ,''          ,'analysis_error'    ,'analysis_error.analysis_error'    ,'analysis error'    ,LaLo   ,LaLo   ,LaLo   ,12,'Kelvin';      
    'mask'              ,''          ,'mask'              ,'mask.mask'                        ,'mask'              ,LaLo   ,LaLo   ,LaLo   ,14,'';
    };


DSI.Range.Time1_NODC=datenum(2005,01,233);
DSI.Range.Time2_NODC=datenum(2009,01,238); % update
DSI.Range.Time2_JPL=now; 
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/GLOB/REMSS/mw_ir_OI';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/REMSS/mw_ir_OI';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2005/233/';
%DSI.URLFILE='20050821-REMSS-L4HRfnd-GLOB-v01-fv01-mw_ir_OI.nc.gz';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-REMSS-L4HRfnd-GLOB-v01-'' ''fv01-mw_ir_OI.nc.gz''];';
DSI.TimeofFrame=[0];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    
 
 %           long_name: '"analysed sea surface temperature"'
 %       standard_name: '"sea_surface_temperature"'
 %                type: '"foundation"'
 %               units: '"kelvin"'
 %       ml__FillValue: -32768
 %          add_offset: 273.1500
 %        scale_factor: 0.0100
 %           valid_min: -300
 %           valid_max: 4500
 %              source: '"REMSS"'
 %   DODS_ML_Real_Name: 'analysed_sst'

 %   long_name: '"estimated error standard deviation of analysed_sst"'
 %               units: '"kelvin"'
 %       ml__FillValue: -32768
 %          add_offset: 0
 %        scale_factor: 0.0100
 %           valid_min: 0
 %           valid_max: 32767
 %   DODS_ML_Real_Name: 'analysis_error'

 %   long_name: '"sea/land/lake/ice field composite mask"'
 %       ml__FillValue: 128
 %         flag_values: [4x1 double]
 %       flag_meanings: '"sea land lake ice"'
 %             comment: [1x209 char]
 %   DODS_ML_Real_Name: 'mask'
DSI.FormulaCnv.mask='';
DSI=DSI_common3(DSI);
SaveDSI(DSI)

function DSI=DSI_GHRSST_GLOB_UKMO
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_GLOB_UKMO';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/UKMO/OSTIA/
%2006/091/20060401-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc.bz2?
%time[0:1:0],lat[0:1:3599],lon[0:1:7199],
%analysed_sst[0:1:0][0:1:3599][0:1:7199],
%analysis_error[0:1:0][0:1:3599][0:1:7199],
%sea_ice_fraction[0:1:0][0:1:3599][0:1:7199],
%mask[0:1:0][0:1:3599][0:1:7199]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name          opendap dataset returned 
%    name                              dataset                                                                     coord   coord   coord    CBNum
    'analysed_sst'      ,''          ,'analysed_sst'      ,'analysed_sst.analysed_sst'        ,'analysed sst'      ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'    ,''          ,'analysis_error'    ,'analysis_error.analysis_error'    ,'analysis error'    ,LaLo   ,LaLo   ,LaLo   ,12,'Kelvin';      
    'sea_ice_fraction'  ,''          ,'sea_ice_fraction'  ,'sea_ice_fraction.sea_ice_fraction','sea ice fraction'  ,LaLo   ,LaLo   ,LaLo   ,13,'';    
    'mask'              ,''          ,'mask'              ,'mask.mask'                        ,'mask'              ,LaLo   ,LaLo   ,LaLo   ,14,'';
    };

DSI.Range.Time1_NODC=datenum(2006,01,91);
DSI.Range.Time2_NODC=datenum(2009,01,238); % update
DSI.Range.Time2_JPL=now; 
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/GLOB/UKMO/OSTIA';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/UKMO/OSTIA';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2006/091/';
%DSI.URLFILE='20060401-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc.bz2';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-UKMO-L4HRfnd-GLOB-v01-'' ''fv02-OSTIA.nc.bz2''];';
DSI.TimeofFrame=[0.5];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    
 
  %          long_name: '"analysed sea surface temperature"'
  %      standard_name: '"sea_surface_temperature"'
  %               type: '"foundation"'
  %              units: '"kelvin"'
  %      ml__FillValue: -32768
  %         add_offset: 273.1500
  %       scale_factor: 0.0100
  %          valid_min: -300
  %          valid_max: 4500
  %  DODS_ML_Real_Name: 'analysed_sst'
 
 %           long_name: '"estimated error standard deviation of analysed_sst"'
 %               units: '"kelvin"'
 %       ml__FillValue: -32768
 %          add_offset: 0
 %        scale_factor: 0.0100
 %           valid_min: 0
 %           valid_max: 32767
 %   DODS_ML_Real_Name: 'analysis_error'

  %  long_name: '"sea ice fraction"'
  %      standard_name: '"sea_ice_area_fraction"'
  %              units: '"1"'
  %      ml__FillValue: 128
  %         add_offset: 0
  %       scale_factor: 0.0100
  %          valid_min: 0
  %          valid_max: 100
  %             source: '"EUMETSAT OSI-SAF"'
  %  DODS_ML_Real_Name: 'sea_ice_fraction'

   % long_name: '"sea/land/lake/ice field composite mask"'
   %     ml__FillValue: 128
   %       flag_values: [4x1 double]
   %     flag_meanings: '"sea land lake ice"'
   %           comment: [1x224 char]
   % DODS_ML_Real_Name: 'mask'
DSI.FormulaCnv.mask='';
DSI=DSI_common3(DSI);
SaveDSI(DSI)

function DSI=DSI_GHRSST_MED
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_MED';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/MED/EUR/
%2006/355/20061221-EUR-L4UHFnd-MED-v01.nc.bz2?
%time[0:1:0],skin_time[0:1:11],lon[0:1:2749],lat[0:1:824],
%sst_foundation[0:1:0][0:1:824][0:1:2749],
%normalised_analysis_error[0:1:0][0:1:824][0:1:2749],
%bias[0:1:0][0:1:824][0:1:2749],
%sea_ice_fraction[0:1:0][0:1:824][0:1:2749],
%mask[0:1:0][0:1:824][0:1:2749],
%DT_sst_skin[0:1:11][0:1:824][0:1:2749],
%sst_skin_quality_flag[0:1:824][0:1:2749]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};
TLaLo={'time','Latitude','Longitude'};
LaLoT={'Latitude','Longitude','time'};

DSI.FieldsTable={
%    opendap                  file         name in                    returned name                                         menu name                    opendap dataset returned 
%    name                     name         dataset                                                                                                       coord   coord   coord    CBNum
    'sst_foundation'           ,''        ,'sst_foundation'           ,'sst_foundation.sst_foundation'                       ,'sst foundation'           ,LaLo   ,LaLo   ,LaLo   ,15,'Kelvin';
    'normalised_analysis_error',''        ,'normalised_analysis_error','normalised_analysis_error.normalised_analysis_error' ,'normalised analysis error',LaLo   ,LaLo   ,LaLo   ,16,'percent';   
    'bias'                     ,''        ,'bias'                     ,'bias.bias'                                           ,'analysis error bias'      ,LaLo   ,LaLo   ,LaLo   ,17,'Kelvin';         
    'sea_ice_fraction'         ,''        ,'sea_ice_fraction'         ,'sea_ice_fraction.sea_ice_fraction'                   ,'sea ice fraction'         ,LaLo   ,LaLo   ,LaLo   ,13,'percent';    
    'mask'                     ,''        ,'mask'                     ,'mask.mask'                                           ,'mask'                     ,LaLo   ,LaLo   ,LaLo   ,14,'';
%???    'DT_sst_skin'              ,''        ,'DT_sst_skin'              ,'DT_sst_skin.DT_sst_skin'                             ,'skin sst'                 ,LaLo   ,TLaLo  ,LaLoT  ,18,'Kelvin';   
%    'sst_skin_quality_flag'    ,''        ,'sst_skin_quality_flag'    ,'sst_skin_quality_flag.sst_skin_quality_flag'         ,'sst skin quality flag'    ,LaLo   ,LaLo   ,LaLo  ,19,'';
    };

% warning: disabled
% DT_sst_skin[0:1:11][0:1:824][0:1:2749]  strange: one more dimension with 12 values ??? 
% sst_skin_quality_flag[0:1:824][0:1:2749] does not allow CTIME='[0]'; 
% 

DSI.Range.Time1_NODC=datenum(2005,01,116);
DSI.Range.Time2_NODC=datenum(2008,01,000); % 2008-01-13 ?
DSI.Range.Time2_JPL=DSI.Range.Time2_NODC; %now; % followed by ODYSSEA
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/MED/EUR';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/MED/EUR';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2006/355/';
%DSI.URLFILE='20061221-EUR-L4UHFnd-MED-v01.nc.bz2';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-EUR-L4UHFnd-MED-'' ''v01.nc.bz2''];';
DSI.TimeofFrame=[0.];

YYYY='2008' %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    
 
 %           long_name: '"foundation sea surface temperature"'
 %               units: '"kelvin"'
 %       ml__FillValue: -32768
 %          add_offset: 273.1500
 %        scale_factor: 0.0100
 %           valid_min: -32767
 %           valid_max: 32767
 %   DODS_ML_Real_Name: 'sst_foundation'

  %                  long_name: '"error estimate from the analysis"'
  %                      units: '"percent"'
  %              ml__FillValue: 128
  %                 add_offset: 0
  %               scale_factor: 0.0100
  %                  valid_min: 129
  %                  valid_max: 127
  %          DODS_ML_Real_Name: 'normalised_analysis_error'
   
   %         long_name: '"analysis error bias"'
   %             units: '"kelvin"'
   %     ml__FillValue: 128
   %        add_offset: 0
   %      scale_factor: 0.0100
   %         valid_min: 129
   %         valid_max: 127
   % DODS_ML_Real_Name: 'bias'
   
  %          long_name: '"sea ice fraction"'
  %              units: '"percent"'
  %      ml__FillValue: 128
  %         add_offset: 0
  %       scale_factor: 0.0100
  %          valid_min: 129
  %          valid_max: 127
  %  DODS_ML_Real_Name: 'sea_ice_fraction'
 
  %          long_name: '"Land/Sea/Ice mask"'
  %      ml__FillValue: 128
  %            comment: '"b0:1=sea; b1:1=land; b2:1=lake; b3:1=ice"'
  %  DODS_ML_Real_Name: 'mask'
DSI.FormulaCnv.mask='';

   % long_name: '"skin sea surface temperature"'
   %             units: '"kelvin"'
   %     ml__FillValue: 128
   %        add_offset: 0
   %      scale_factor: 0.1000
   %         valid_min: -127
   %         valid_max: 127
   %            source: '"Stuart-Menteth model"'
   % DODS_ML_Real_Name: 'DT_sst_skin'
   
  %              long_name: '"quality control indicator for sst skin"'
  %          ml__FillValue: -2.1475e+009
  %                comment: [1x185 char]
  %      DODS_ML_Real_Name: 'sst_skin_quality_flag'
FIELD='sst_skin_quality_flag';
FIELDDS='sst_skin_quality_flag';
eval(['DSI.FormulaNaN.' FIELD '=''f(find(f<0.95*A.' FIELDDS '.ml__FillValue))=NaN;'';']);
DSI.FormulaCnv.sst_skin_quality_flag='';
% a.sst_skin_quality_flag.comment
%ans =
%"each pair of bits is related to one of the 12  DT_sst_skin values in chronological order, with the following code convention:bit(i+1,i):00 = good, 01 = fair, 10 = uncertain, 11 = poor"
DSI=DSI_common3(DSI);
SaveDSI(DSI)

function DSI=DSI_GHRSST_MED_ODYSSEA
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_MED_ODYSSEA';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/MED/EUR/ODYSSEA/
%2008/001/20080101-EUR-L4UHRfnd-MED-v01-fv01-ODYSSEA.nc.bz2?
%time[0:1:0],lat[0:1:824],lon[0:1:2749],
%analysed_sst[0:1:0][0:1:824][0:1:2749],
%analysis_error[0:1:0][0:1:824][0:1:2749],
%mask[0:1:0][0:1:824][0:1:2749],
%sea_ice_fraction[0:1:0][0:1:824][0:1:2749]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name                 opendap dataset returned 
%    name                              dataset                                                                             coord   coord   coord    CBNum
    'analysed_sst'         ,''        ,'analysed_sst'         ,'analysed_sst.analysed_sst'       ,'analysed sst'           ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'       ,''        ,'analysis_error'       ,'analysis_error.analysis_error'   ,'analysis error'         ,LaLo   ,LaLo   ,LaLo   ,12,'percent';
    'sea_ice_fraction'     ,''        ,'sea_ice_fraction'    ,'sea_ice_fraction.sea_ice_fraction','sea ice fraction'       ,LaLo   ,LaLo   ,LaLo   ,13,'';    
    'mask'                 ,''        ,'mask'                ,'mask.mask'                        ,'mask'                   ,LaLo   ,LaLo   ,LaLo   ,14,'';
};

DSI.Range.Time1_NODC=datenum(2008,01,001);
DSI.Range.Time2_NODC=datenum(2009,01,239); % update
DSI.Range.Time2_JPL=now; 
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/MED/EUR/ODYSSEA';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/MED/EUR/ODYSSEA';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2008/001/';
%DSI.URLFILE='20080101-EUR-L4UHRfnd-MED-v01-fv01-ODYSSEA.nc.bz2';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-EUR-L4UHRfnd-MED-v01-'' ''fv01-ODYSSEA.nc.bz2''];';
DSI.TimeofFrame=[0.];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    
 
 %          long_name: '"analysed sea surface temperature"'
 %       standard_name: '"sea_surface_temperature"'
 %                type: '"foundation"'
 %               units: '"kelvin"'
 %       ml__FillValue: -32768
 %          add_offset: 273.1500
 %        scale_factor: 0.0100
 %           valid_min: -300
 %           valid_max: 4500
 %   DODS_ML_Real_Name: 'analysed_sst'

  %  long_name: '"estimated error percentage"'
  %              units: '"percent"'
  %      ml__FillValue: -32768
  %         add_offset: 0
  %       scale_factor: 0.0100
  %          valid_min: 0
  %          valid_max: 10000
  %  DODS_ML_Real_Name: 'analysis_error'

  %          long_name: '"sea/land/lake/ice field composite mask"'
  %        flag_values: '"b0, 1=sea, b1: 1=land, b2: 1=lake, b3: 1=ice"'
  %      flag_meanings: '"sea land lake ice"'
  %  DODS_ML_Real_Name: 'mask'
DSI.FormulaNaN.mask='';
DSI.FormulaCnv.mask='';
 
  %          long_name: '"sea ice area fraction"'
  %              units: '" "'
  %      ml__FillValue: 128
  %         add_offset: 0
  %       scale_factor: 0.0100
  %          valid_min: 0
  %          valid_max: 100
  %        institution: '"SAF O&SI"'
  %  DODS_ML_Real_Name: 'sea_ice_fraction'
DSI=DSI_common3(DSI);
SaveDSI(DSI)

function DSI=DSI_GHRSST_NSEABALTIC
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_NSEABALTIC';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/NSEABALTIC/DMI/DMI_OI/
%2007/155/20070604-DMI-L4UHfnd-NSEABALTIC-v01-fv01-DMI_OI.nc.gz?
%lon[0:1:1333],lat[0:1:599],time[0:1:0],
%analysed_sst[0:1:0][0:1:599][0:1:1333],
%analysis_error[0:1:0][0:1:599][0:1:1333],
%mask[0:1:0][0:1:599][0:1:1333]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name          opendap dataset returned 
%    name                              dataset                                                                     coord   coord   coord    CBNum
    'analysed_sst'      ,''          ,'analysed_sst'      ,'analysed_sst.analysed_sst'        ,'analysed sst'      ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'    ,''          ,'analysis_error'    ,'analysis_error.analysis_error'    ,'analysis error'    ,LaLo   ,LaLo   ,LaLo   ,12,'Kelvin';      
    'mask'              ,''          ,'mask'              ,'mask.mask'                        ,'mask'              ,LaLo   ,LaLo   ,LaLo   ,14,'';
    };

DSI.Range.Time1_NODC=datenum(2007,01,155);
DSI.Range.Time2_NODC=datenum(2009,01,239); % update
DSI.Range.Time2_JPL=now; 
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/NSEABALTIC/DMI/DMI_OI';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/NSEABALTIC/DMI/DMI_OI';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2007/155/';
%DSI.URLFILE='20070604-DMI-L4UHfnd-NSEABALTIC-v01-fv01-DMI_OI.nc.gz';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-DMI-L4UHfnd-NSEABALTIC-v01-'' ''fv01-DMI_OI.nc.gz''];';
DSI.TimeofFrame=[0.];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    

  %     standard_name: '"sea_surface_temperature"'
  %          long_name: '"analysed sea surface temperature"'
  %              units: '"kelvin"'
  %          valid_min: -300
  %          valid_max: 4500
  %       scale_factor: 0.0100
  %         add_offset: 273.1500
  %      ml__FillValue: -32768
  %               type: '"foundation"'
  %  DODS_ML_Real_Name: 'analysed_sst'

   %         long_name: '"estimated error standard deviation of analysed_sst"'
   %             units: '"kelvin"'
   %         valid_min: 0
   %         valid_max: 32767
   %      scale_factor: 0.0100
   %        add_offset: 0
   %     ml__FillValue: -32768
   % DODS_ML_Real_Name: 'analysis_error'

    %        long_name: '"sea/land/lake/ice field composite mask"'
    %    ml__FillValue: 128
    %          comment: [1x208 char]
    %      flag_values: '"1b,2b,4b,8b"'
    %    flag_meanings: '"sea land lake ice"'
   % DODS_ML_Real_Name: 'mask'
DSI.FormulaCnv.mask='';

DSI=DSI_common3(DSI);
SaveDSI(DSI)

function DSI=DSI_GHRSST_NWE
DSI.DataSetName='GHRSST';
DSI.DataSetBranch='GHRSST_NWE';

%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/NWE/EUR/ODYSSEA/
%2008/023/20080123-EUR-L4UHRfnd-NWE-v01-fv01-ODYSSEA.nc.bz2?
%time[0:1:0],lat[0:1:849],lon[0:1:1099],
%analysed_sst[0:1:0][0:1:849][0:1:1099],
%analysis_error[0:1:0][0:1:849][0:1:1099],
%mask[0:1:0][0:1:849][0:1:1099],
%sea_ice_fraction[0:1:0][0:1:849][0:1:1099]

% define list of FieldsTable or variables available in the dataset
LaLo={'latitude','longitude'};

DSI.FieldsTable={
%    opendap          file name        name in                returned name                     menu name          opendap dataset returned 
%    name                              dataset                                                                     coord   coord   coord    CBNum
    'analysed_sst'      ,''          ,'analysed_sst'      ,'analysed_sst.analysed_sst'        ,'analysed sst'      ,LaLo   ,LaLo   ,LaLo   ,11,'Kelvin';
    'analysis_error'    ,''          ,'analysis_error'    ,'analysis_error.analysis_error'    ,'analysis error'    ,LaLo   ,LaLo   ,LaLo   ,12,'percent';      
    'sea_ice_fraction'  ,''          ,'sea_ice_fraction'  ,'sea_ice_fraction.sea_ice_fraction','sea ice fraction'  ,LaLo   ,LaLo   ,LaLo   ,13,'';    
    'mask'              ,''          ,'mask'              ,'mask.mask'                        ,'mask'              ,LaLo   ,LaLo   ,LaLo   ,14,'';
    };

DSI.Range.Time1_NODC=datenum(2008,01,023);
DSI.Range.Time2_NODC=datenum(2009,01,239); % update
DSI.Range.Time2_JPL=now; 
DSI.Range.Time1=DSI.Range.Time1_NODC;
DSI.Range.Time2=DSI.Range.Time2_JPL;
DSI.URLSITE_NODC='http://data.nodc.noaa.gov/opendap/ghrsst/L4/NWE/EUR/ODYSSEA';
DSI.URLSITE_JPL='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/NWE/EUR/ODYSSEA';
DSI.URLSITE='URLSITE=DSI.URLSITE_NODC; if Time > DSI.Range.Time2_NODC URLSITE=DSI.URLSITE_JPL; end';

%DSI.URLPATH='/2008/023/';
%DSI.URLFILE='20080123-EUR-L4UHRfnd-NWE-v01-fv01-ODYSSEA.nc.bz2';
DSI.URLPATH='URLPATH=[''/'' YYYY ''/'' DDD ''/''];';
DSI.URLFILE='URLFILE=[YYYYMMDD ''-EUR-L4UHRfnd-NWE-v01-'' ''fv01-ODYSSEA.nc.bz2''];';
DSI.TimeofFrame=[0.];

YYYY='2009'; %update once a year
URL=[DSI.URLSITE_NODC '/' YYYY '/']
ddd=LastNODC(URL)
DSI.Range.Time2_NODC=datenum(str2num(YYYY),01,ddd-1); 

DSI=DSICommon(DSI);
% FormulaNaN and FormulaCnv are standard except for flags: set them empty    

 
  %        long_name: '"analysed sea surface temperature"'
  %      standard_name: '"sea_surface_temperature"'
  %               type: '"foundation"'
  %              units: '"kelvin"'
  %      ml__FillValue: -32768
  %         add_offset: 273.1500
  %       scale_factor: 0.0100
  %          valid_min: -300
  %          valid_max: 4500
  %  DODS_ML_Real_Name: 'analysed_sst'

   % long_name: '"estimated error percentage"'
   %             units: '"percent"'
   %     ml__FillValue: -32768
   %        add_offset: 0
   %      scale_factor: 0.0100
   %         valid_min: 0
   %         valid_max: 10000
   % DODS_ML_Real_Name: 'analysis_error'

  %          long_name: '"sea/land/lake/ice field composite mask"'
  %        flag_values: '"b0, 1=sea, b1: 1=land, b2: 1=lake, b3: 1=ice"'
  %      flag_meanings: '"sea land lake ice"'
  %  DODS_ML_Real_Name: 'mask'
DSI.FormulaNaN.mask='';
DSI.FormulaCnv.mask='';

   % long_name: '"sea ice area fraction"'
   %             units: '" "'
   %     ml__FillValue: 128
   %        add_offset: 0
   %      scale_factor: 0.0100
   %         valid_min: 0
   %         valid_max: 100
   %       institution: '"SAF O&SI"'
   % DODS_ML_Real_Name: 'sea_ice_fraction'
DSI=DSI_common3(DSI);
SaveDSI(DSI)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F=DSICommon(DSI)
% do things common to all DSI

DSI.Fields=DSI.FieldsTable(:,1);
DSI.Fields_NameLong=DSI.FieldsTable(:,3);
DSI.Fields_NameDS  =DSI.FieldsTable(:,3);
DSI.Fields_NameReturned=DSI.FieldsTable(:,4);
DSI.Fields_NameMenu=DSI.FieldsTable(:,5);
DSI.Fields_Coord   =DSI.FieldsTable(:,6);
DSI.Fields_CoordDS =DSI.FieldsTable(:,7);
DSI.Fields_CoordReturned =DSI.FieldsTable(:,8);
DSI.Fields_CBNum   =DSI.FieldsTable(:,9);
DSI.Fields_Units   =DSI.FieldsTable(:,10);

Nt=length(DSI.Range.Time1:DSI.Range.Time2)*length(DSI.TimeofFrame); % at most 2 frames per day
kt=0;
%DSI.URLPATH=cell(Nt,1);
%DSI.URLFILE=cell(Nt,1);
DSI.URLCTIME=cell(Nt,1);
%DSI.TIME=cell(Nt,1);
DSI.Time=zeros(Nt,1);

for Time=DSI.Range.Time1:DSI.Range.Time2
%YYYY=datestr(jd,'yyyy');
%yearday=(jd-datenum(YYYY,'yyyy')+1);
%DDD=sprintf('%03d',floor(yearday));
%YYYYMMDD=datestr(jd,'yyyymmdd');
%URLPATH=eval(DSI.FormulaURLPATH); % ['/' YYYY '/' DDD '/'];
%URLFILE=eval(DSI.FormulaURLFILE); %[YYYYMMDD '-AMSRE-REMSS-L2P_GRIDDED_25-amsre_' YYYYMMDD 'v5-v01.nc.gz'];
        for jt=1:length(DSI.TimeofFrame)
        kt=kt+1;
%        DSI.URLPATH{kt}=URLPATH;
%        DSI.URLFILE{kt}=URLFILE;
        DSI.URLCTIME{kt}=['[' num2str(jt-1) ']'];                
        t=Time+DSI.TimeofFrame(jt);
        T=datestr(t,'yyyy-mm-dd HH:MM:SS');
        DSI.Time(kt)=t;  
%        DSI.TIME{kt}=T;  
        end
end
if length(DSI.TimeofFrame)==1
    DSI.CTIME='[0]';
end
%DSI.URLPATH=DSI.URLPATH(1:kt);
%DSI.URLFILE=DSI.URLFILE(1:kt);
%DSI.URLCTIME=DSI.URLCTIME(1:kt);
%DSI.Time=DSI.Time(1:kt);  
%DSI.TIME=DSI.TIME(1:kt);  

%DSI.Range.Time1=DSI.Time(1);
%DSI.Range.Time2=DSI.Time(end);

Time=DSI.Range.Time1;
YYYY=datestr(Time,'yyyy');
yearday=(Time-datenum(YYYY,'yyyy')+1);
DDD=sprintf('%03d',floor(yearday));
YYYYMMDD=datestr(Time,'yyyymmdd');
URLSITE=DSI.URLSITE;if strfind(URLSITE,'URLSITE=') eval(URLSITE); end
URLPATH=DSI.URLPATH;if strfind(URLPATH,'URLPATH=') eval(URLPATH); end % ['/' YYYY '/' DDD '/'];
URLFILE=DSI.URLFILE;if strfind(URLFILE,'URLFILE=') eval(URLFILE); end %[YYYYMMDD '-AMSRE-REMSS-L2P_GRIDDED_25-amsre_' YYYYMMDD 'v5-v01.nc.gz'];

if strcmp(DSI.DataSetBranch,'XXX_AMSRE')
%NC_GLOBAL.southernmost_latitude: -89.87500000
%NC_GLOBAL.northernmost_latitude: 89.87500000
%NC_GLOBAL.westernmost_longitude: 179.8750000
%NC_GLOBAL.easternmost_longitude: -179.8750000
DSI.Latitude=(-89.875:0.25:89.875)';
DSI.Longitude=(-179.875:0.25:179.875)';
else
URL1=[URLSITE URLPATH URLFILE]
loaddap([URL1 '?lat,lon'])
    if exist('lat')
    DSI.Latitude=lat;
    DSI.Longitude=lon;
    elseif exist('sea_surface_temperature')
    DSI.Latitude=sea_surface_temperature.lat;
    DSI.Longitude=sea_surface_temperature.lon;
    elseif exist('analysed_sst')
    DSI.Latitude=analysed_sst.lat;
    DSI.Longitude=analysed_sst.lon;
    elseif exist('sst_foundation')
    DSI.Latitude=sst_foundation.lat;
    DSI.Longitude=sst_foundation.lon;
    end
end
DSI.Range.Latitude1=min(DSI.Latitude);
DSI.Range.Latitude2=max(DSI.Latitude);
DSI.Range.Longitude1=min(DSI.Longitude);
DSI.Range.Longitude2=max(DSI.Longitude);

% set FormulaNaN and FormulaCnv to standard form
for j=1:length(DSI.Fields)
    FIELD=DSI.Fields{j};
    FIELDDS=DSI.Fields_NameDS{j};
    eval(['DSI.FormulaNaN.' FIELD '=''f(find(f==A.' FIELDDS '.ml__FillValue))=NaN;'';']);
    eval(['DSI.FormulaCnv.' FIELD '=''f=A.' FIELDDS '.scale_factor*f+A.' FIELDDS '.add_offset;'';']);
end
% if they are non standard change them in the calling function

S1=['time in seconds since reference time;'];
S2=['missing values replaced with NaN;'];
S3=['variables converted to physical units.'];
%S4='';%['time is MATLAB datenum.'];
DSI.Readme=[S1 '/' S2 '/' S3];

F=DSI;

function SaveDSI(DSI)
eval(['DSI_' DSI.DataSetBranch '=DSI;'])
dsi_filename = 'dsi_ghrsst.mat';
%if (~isdeployed)
    pname = mfilename('fullpath');
    pname = pname(1:length(pname) - length(mfilename));
    dsi_filename = [pname,dsi_filename];
%end
if ~exist(dsi_filename,'file')
save(dsi_filename,['DSI_' DSI.DataSetBranch])
end
save(dsi_filename,['DSI_' DSI.DataSetBranch],'-append')


function DSI=DSI_common3(DSI)
CoordinatesTable={
'latitude','latitude','degrees_north';
'longitude','longitude','degrees_east';
%'depth','depth','m';
'time','time','days since 0000-01-01 00:00:00';
};

DSI.Coordinates=CoordinatesTable(:,1);
DSI.Dimensions.latitude=length(DSI.Latitude);
DSI.Dimensions.longitude=length(DSI.Longitude);
%DSI.Dimensions.depth=length(DSI.Depth);
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

function ddd=LastNODC(URL)
%URL='http://data.nodc.noaa.gov/opendap/ghrsst/L2P/AMSRE/REMSS/2009/';
S=urlread(URL);
%S=regexp(S,[FN ''.*?'' EXT],''match'')
S1='<td align="left">'; 
S2='</td>';

S1='<a href="';
S2='/contents.html">';
S=regexp(S,[S1 '.*?' S2],'match');

%cut the the leading and trailing strings 
L1=length(S1);
L2=length(S2);
SS=S;
for k=1:length(S)
S1=SS{k};    
S{k}=S1(L1+1:end-L2);
end
%S(:);
ddd=str2num(S{end});



