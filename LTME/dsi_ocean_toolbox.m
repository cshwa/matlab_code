% to be placed in directory Common
 function dsi_ocean_toolbox(varargin)
% DataSetInventory and Interface

DSI.DataSetsTable={
% DataSetName,       PBNum,            tooltip,                                     DocumentationURL    
    'AIRS',             1,'Parameters from Atmospheric Infrared Sounder (AIRS)',    'http://disc.gsfc.nasa.gov/AIRS/airsL3_STD.shtml'
    'AVISO_Altimetry',  2,'AVISO Altimetry from TOPEX, Jason, ERS and Merged',      ''
    'GHRSST',           3,'Global High Resolution Sea Surface Temperature (GHRSST)',''
    'GOES',             4,'Sea Surface Temperature from Geostationary Operational Environmental Satellite (GOES) (Regional)','http://coastwatch.pfeg.noaa.gov/infog/GA_ssta_las.html'
    'HYCOM',            5,'HYbrid Coordinate Ocean Model (HYCOM)',                  'http://hycom.coaps.fsu.edu/data/atlantic_info.html'
    'MODIS',            6,'Sea Surface Temperature from Moderate Resolution Imaging Spectroradiometer (MODIS)',''
    'OAFlux',           7,'Objectively Analyzed Air-Sea Flux (OAFlux)',             'http://oaflux.whoi.edu'
    'OceanColor',       8,'Ocean Color from CZCS, SeaWIFS, MODIS and Merged SeaWIFS-MODIS','oceancolor.gsfc.nasa.gov/PRODUCTS/sst.html'
    'OceanWind',        9,'10-m Surface Winds from Scatterometer Missions','http://dods.jpl.nasa.gov/opendap/ocean_wind/ccmp/L3.0/doc/contents.html'                
    'Pathfinder1km',   10,'Sea Surface Temperature (1-km) from AVHRR Using Pathfinder Algorithm (Regional)','http://satdat1.gso.uri.edu/opendap/Pathfinder/Pathfinder1km/pathfinder_1km.html'
    'Pathfinder4km',   11,'Sea Surface Temperature (4-km) from AVHRR Using Pathfinder Algorithm','http://www.nodc.noaa.gov/sog/pathfinder4km/userguide.html' 
        };

% OceanWind documentation url also 'http://podaac-www.jpl.nasa.gov/PRODUCTS/p109.html'

DSI.VariablesTable={
% variable,                   CBNum,  datasets containings this variable    
    'wind',                     1,    {'OceanWind','AVISO_Altimetry','GHRSST','OAFlux'}
    'wind stress',              2,    {'OceanWind'}
%    '',              2,    {'OceanWind'}
    'air temperature',          3,    {'AIRS'}
    'air pressure',             4,    {'AIRS'}
    'cloud cover',              5,    {'AIRS'}
    'water vapor',              6,    {'AIRS'}
    'heat flux',                7,    {'HYCOM','OAFlux'}
    'temperature',              8,    {'HYCOM'}
    'salinity',                 9,    {'HYCOM'}
    'density',                  10,   {'HYCOM'}
    'currents',                 11,   {'HYCOM','AVISO_Altimetry'}
    'mixed layer depth',        12,   {'HYCOM'}
    'sea surface temperature',  13,   {'GOES','GHRSST','HYCOM','OceanColor','Pathfinder1km','Pathfinder4km','MODIS'}
    'sea surface salinity',     14,   {''}
    'sea surface height',       15,   {'AVISO_Altimetry','GHRSST','HYCOM'}
    'ocean color',              16,   {'OceanColor'}
    'ice',                      17,   {'GHRSST'}
    'other',                    18,   {'AIRS','AVISO_Altimetry','HYCOM','OceanColor','Pathfinder4km'}
    };
DSI.TemporalTable={
% parameter,                   checkbox,                   datasets containings this parameter    
    'climatology',            'checkboxClimatology',   {'Pathfinder4km'}
    'time series',            'checkboxTimeSeries',    {'AIRS','AVISO_Altimetry','GHRSST','GOES','HYCOM','MODIS','OAFlux','OceanColor','Pathfinder1km','Pathfinder4km','OceanWind'}
};
DSI.SourceTable={
% parameter,                   checkbox,                  datasets containings this parameter    
    'satellite',              'checkboxSatellite',    {'AIRS','AVISO_Altimetry','GHRSST','GOES','MODIS','OAFlux','OceanColor','Pathfinder1km','Pathfinder4km','OceanWind'}
    'model',                  'checkboxModel',        {'HYCOM','OAFlux'}
};


% get range of time,lat,lon from individual dataset dsi
    for k=1:length(DSI.DataSetsTable(:,1))
        DSName=DSI.DataSetsTable{k,1};
        disp(DSName)
        dsname=lower(DSName);
% fore each:  load('../Pathfinder4km/dsi_pathfinder4km.mat')
        load(['../' DSName '/dsi_' dsname '.mat'])
    end
    
% extract DSI_$DataSetBranch$.Range information    
C=who('DSI_*');
    for k=1:length(C)
        S=C{k}; % 'DSI_Pathfinder4km_daily'
%        fieldnames(eval(S))
        %eval(S)
        DataSetBranch=S(5:end); % 'Pathfinder4km_daily'
        %eval(['DSI_' DataSetBranch ])
        if isfield(eval(S),'Range')
            eval(['DSI.' DataSetBranch '.Range = DSI_' DataSetBranch '.Range;']) % DSI.Pathfinder4km_daily=DSI_Pathfinder4km_daily.Range;       
        end
        if isfield(eval(S),'Ranges')
            eval(['DSI.' DataSetBranch '.Range = DSI_' DataSetBranch '.Ranges;']) % DSI.Pathfinder4km_daily=DSI_Pathfinder4km_daily.Range;       
        end
    end
    
% combine DataSetBranch.Range into DataSet.Range
       
   for k=1:length(DSI.DataSetsTable(:,1))
        DSName=DSI.DataSetsTable{k,1};
        % initialize to calculate min and max
        TimeMin=datenum(9999,12,31);TimeMax=0;LatMin=180;LatMax=-180;LonMin=720;LonMax=-720;
        %C=who(['DSI_' DSName '_*'); % DataSetBranches
        C=fieldnames(DSI); C=C(strmatch(DSName,C)); % select only Branches corresponding to a particular DataSet
        for kk=1:length(C)
            S=C{kk}; % 'Pathfinder4km_daily'
            %disp(S)
            eval(['Time1=DSI.' S '.Range.Time1;'])
            eval(['Time2=DSI.' S '.Range.Time2;'])
            eval(['Lat1=DSI.' S '.Range.Latitude1;'])
            eval(['Lat2=DSI.' S '.Range.Latitude2;'])
            eval(['Lon1=DSI.' S '.Range.Longitude1;'])
            eval(['Lon2=DSI.' S '.Range.Longitude2;'])
            if isempty(strfind(S,'limatology')) % exclude Climatology
            TimeMin=min(TimeMin,Time1);
            %datestr(TimeMin,'yyyy-mm-dd')
            TimeMax=max(TimeMax,Time2);
            end
            LatMin=min(LatMin,Lat1);
            LatMax=max(LatMax,Lat2);
            LonMin=min(LonMin,Lon1);
            LonMax=max(LonMax,Lon2);        
        end
            eval(['DSI.' DSName '.Range.Time1=TimeMin;'])
            eval(['DSI.' DSName '.Range.Time2=TimeMax;'])
            eval(['DSI.' DSName '.Range.Latitude1=LatMin;'])
            eval(['DSI.' DSName '.Range.Latitude2=LatMax;'])
            eval(['DSI.' DSName '.Range.Longitude1=LonMin;'])
            eval(['DSI.' DSName '.Range.Longitude2=LonMax;'])
   end

   % combine DataSet.Range into Ocean_Toolbox.Range

           % initialize to calculate min and max
   TimeMin=datenum(9999,12,31);TimeMax=0;LatMin=180;LatMax=-180;LonMin=720;LonMax=-720;
   for k=1:length(DSI.DataSetsTable(:,1))
            S=DSI.DataSetsTable{k,1};
            eval(['Time1=DSI.' S '.Range.Time1;'])
            eval(['Time2=DSI.' S '.Range.Time2;'])
            eval(['Lat1=DSI.' S '.Range.Latitude1;'])
            eval(['Lat2=DSI.' S '.Range.Latitude2;'])
            eval(['Lon1=DSI.' S '.Range.Longitude1;'])
            eval(['Lon2=DSI.' S '.Range.Longitude2;'])
            TimeMin=min(TimeMin,Time1);
            TimeMax=max(TimeMax,Time2);
            LatMin=min(LatMin,Lat1);
            LatMax=max(LatMax,Lat2);
            LonMin=min(LonMin,Lon1);
            LonMax=max(LonMax,Lon2);        
   end
            S='Ocean_Toolbox';
            eval(['DSI.' S '.Range.Time1=TimeMin;'])
            eval(['DSI.' S '.Range.Time2=TimeMax;'])
            eval(['DSI.' S '.Range.Latitude1=LatMin;'])
            eval(['DSI.' S '.Range.Latitude2=LatMax;'])
            eval(['DSI.' S '.Range.Longitude1=LonMin;'])
            eval(['DSI.' S '.Range.Longitude2=LonMax;'])   
%DSI
%clear DSI_*    

save('dsi_ocean_toolbox.mat','DSI')
