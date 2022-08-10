function dsi_pathfinder1km(varargin)
% DataSetInventory and Interface for the Pathfinder1km 
%
% Meri Sheremet
% March 2009

% DataSetBranch=function of (SpatialCoverage)
List.SpatialCoverage={
'NWA',                        'Northwest Atlantic';
'NEP',                        'Northeast Pacific';
'NWP',                        'Northwest Pacific';    
};

if nargin < 1 
List.DataSetBranch={    
%Satellite_SpatialCoverage
'Pathfinder1km_NWA';
'Pathfinder1km_NEP';
'Pathfinder1km_NWP';
};
    
    for k=1:length(List.DataSetBranch)
%    disp(['DSI_' List.DataSetBranch{k} ])
    eval(['DSI_' List.DataSetBranch{k} ';'])
    end

    dsi_filename = 'dsi_pathfinder1km.mat';
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


function DSI=DSI_Pathfinder1km_NWA
DSI.DataSetName='Pathfinder1km';
DSI.DataSetBranch='Pathfinder1km_NWA';
DSI.SpatialCoverage='NWA';
DSI.SpatialCoverageFN='NWAtlantic';
DSI=DSI_common2(DSI);

function DSI=DSI_Pathfinder1km_NEP
DSI.DataSetName='Pathfinder1km';
DSI.DataSetBranch='Pathfinder1km_NEP';
DSI.SpatialCoverage='NEP';
DSI.SpatialCoverageFN='NEPacific';
DSI=DSI_common2(DSI);

function DSI=DSI_Pathfinder1km_NWP
DSI.DataSetName='Pathfinder1km';
DSI.DataSetBranch='Pathfinder1km_NWP';
DSI.SpatialCoverage='NWP';
DSI.SpatialCoverageFN='NWPacific';
DSI=DSI_common2(DSI);

function DSI=DSI_common2(DSI)
DSI.URLSITE='http://satdat1.gso.uri.edu';
DSI.URLPATH='/thredds/dodsC/';
SCFN=DSI.SpatialCoverageFN;
LL={'latitude','longitude'};
TLL={'time','lat','lon'};
% define list of fields or variables available in the dataset
FieldsTable={
%    opendap      long            file name    name in        returned name        menu name       opendap   dataset  CBNum   units
%    name         name                         dataset                                             coord     coord     
'sst_dec', 'sst_declouded',[SCFN 'Dec_1km'],'dsp_band_1','dsp_band_1.dsp_band_1','SST (Declouded)' ,LL        ,TLL    ,1      ,'degree_C';
'sst_raw', 'sst_raw',      [SCFN 'Raw_1km'],'dsp_band_1','dsp_band_1.dsp_band_1','SST (Raw)'       ,LL        ,TLL    ,2      ,'degree_C';
     };
DSI.Fields         =FieldsTable(:,1); 
DSI.Fields_NameLong=FieldsTable(:,2);
DSI.Fields_FileName=FieldsTable(:,3);
DSI.Fields_NameDS  =FieldsTable(:,4);
DSI.Fields_NameReturned=FieldsTable(:,5);
DSI.Fields_NameMenu=FieldsTable(:,6);
DSI.Fields_Coord   =FieldsTable(:,7);
DSI.Fields_CoordDS =FieldsTable(:,8);
DSI.Fields_CBNum   =FieldsTable(:,9);
DSI.Fields_Units   =FieldsTable(:,10);

CoordinatesTable={
'latitude','latitude','degrees_north';
'longitude','longitude','degrees_east';
%'depth','depth','m';
'time','time','days since 0000-01-01 00:00:00';
};

for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
eval(['DSI.URLFILE.' FIELD '=DSI.Fields_FileName{k};']);
eval(['DSI.URLCVAR.' FIELD '=DSI.Fields_NameDS{k};']);
end

% get available lat,lon,time 
%loaddap('http://thredds1.pfeg.noaa.gov:8080/thredds/dodsC/satellite/GA/ssta/1day?time')
URL=[DSI.URLSITE DSI.URLPATH DSI.URLFILE.sst_raw '?lat,lon,time'];
loaddap(['+v','-e'],URL);
%dods_err = 0 means no error.
%dods_err = 1 means error.
%dods_err_msg [variable where error message is stored]
if dods_err
    disp(dods_err_msg)
    msgbox('Could not access dataset site. Sorry, try again later.');
    return
end
Time70Raw=round(time);
URL=[DSI.URLSITE DSI.URLPATH DSI.URLFILE.sst_dec '?lat,lon,time'];
loaddap(['+v','-e'],URL);
%dods_err = 0 means no error.
%dods_err = 1 means error.
%dods_err_msg [variable where error message is stored]
if dods_err
    disp(dods_err_msg)
    msgbox('Could not access dataset site. Sorry, try again later.');
    return
end
Time70Dec=round(time);

% merge two time arrays with possible gaps
[Time70,iDec,iRaw]=mergeint(Time70Dec,Time70Raw);

%units: '"seconds since 1970-01-01"'
t0=datenum('1970-01-01 00:00:00','yyyy-mm-dd HH:MM:SS'); 
t=Time70/(24*60*60)+t0; % in datenum days since '0000-01-01 00:00:00'
DSI.Time=t;
DSI.Latitude=lat;
DSI.Longitude=lon;

% Note that Longitude may jump from 180E to -179W in Western Pacific Region
% fix that by adding 360 degrees
if strcmp(DSI.SpatialCoverage,'NWP')
i=find(DSI.Longitude<0);
DSI.Longitude(i)=DSI.Longitude(i)+360;
end

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

% Rename returned variables and change dimensions order
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
    FREN=['f=squeeze(' DSI.Fields_NameReturned{k} ');'];
eval(['DSI.FormulaRen.' FIELD '=''' FREN ''';']);      
end
% bad values 0,1,2,3
DSI.FormulaNaN='f(find(f<4))=NaN;';
%DSI.FormulaOrderCoord='f=f'';';
DSI.FormulaCnv='f=A.dsp_band_1.scale_factor*f + A.dsp_band_1.add_off;';

%DSI.Range.Time1=min(DSI.TimeA);
%DSI.Range.Time2=max(DSI.TimeB);
DSI.Range.Time1=min(DSI.Time);
DSI.Range.Time2=max(DSI.Time);
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

dsi_filename = 'dsi_pathfinder1km.mat';
%if (~isdeployed)
    pname = mfilename('fullpath');
    pname = pname(1:length(pname) - length(mfilename));
    dsi_filename = [pname,dsi_filename];
%end

if ~exist(dsi_filename,'file')
save(dsi_filename,['DSI_' DSI.DataSetBranch])
end
save(dsi_filename,['DSI_' DSI.DataSetBranch],'-append')




function [c,pa,pb]=mergeint(a,b)
% [c,pa,pb]=mergeint(a,b) 
%           merges two integer arrays with gaps
% input:    a, b arrays to be merged
% output:   c merged array containing all possible values in a and b without repetition
%           pa,pb index pointers to original arrays a and b; NaN if the value is not available 
%
% Example:
% a=[010;020;030;050;060;070;080;090;100];
% b=[030;040;060;070;080;090;100;120];
% c=[10;20;30;40;50;60;70;80;90;100;120];
% pa=[     1     2     3   NaN     4     5     6     7     8     9   NaN]';
% pb=[   NaN   NaN     1     2   NaN     3     4     5     6     7     8]';


[NA,N2]=size(a);
[NB,N2]=size(b);
NC=NA+NB;
c=[a;b];pa=c;pb=c;
ka=1;kb=1;kc=0;

while ka<=NA | kb<=NB
if ka<=NA & kb<=NB
    if a(ka)==b(kb)
    kc=kc+1;c(kc)=a(ka);pa(kc)=ka;pb(kc)=kb;ka=ka+1;kb=kb+1;
    elseif a(ka)<b(kb)
    kc=kc+1;c(kc)=a(ka);pa(kc)=ka;pb(kc)=NaN;ka=ka+1;
    elseif a(ka)>b(kb)
    kc=kc+1;c(kc)=b(kb);pa(kc)=NaN;pb(kc)=kb;kb=kb+1;
    end
elseif ka>NA & kb<=NB
    kc=kc+1;c(kc)=b(kb);pa(kc)=NaN;pb(kc)=kb;kb=kb+1;
elseif ka<=NA & kb>NB
    kc=kc+1;c(kc)=a(ka);pa(kc)=ka;pb(kc)=NaN;ka=ka+1;
end
end
NC=kc;
c=c(1:NC);
pa=pa(1:NC);
pb=pb(1:NC);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %examples of URLs for SST (Declouded) and SST (Raw)
%http://satdat1.gso.uri.edu/thredds/dodsC/NWAtlanticDec_1km?dsp_band_1[14497][1498:2:1960][4531:2:4894]
%http://satdat1.gso.uri.edu/thredds/dodsC/NWAtlanticRaw_1km?dsp_band_1[14517][1498:2:1960][884:2:1248]

%A=loaddap('-A','http://satdat1.gso.uri.edu/thredds/dodsC/NWAtlanticDec_1km?dsp_band_1[14497][1498:2:1960][4531:2:4894]')
%A = 
%           dsp_band_1: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.dsp_band_1
%ans = 
%            ml__unsigned: '"true"'
%      ml__CoordinateAxes: '"time lat lon "'
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
%            scale_factor: 0.1250
%                 add_off: -4
%       DODS_ML_Real_Name: 'dsp_band_1'
%              dsp_band_1: [1x1 struct]
%                    time: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
%A=loaddap('-A','http://satdat1.gso.uri.edu/thredds/dodsC/NWAtlanticRaw_1km?dsp_band_1[14517][1498:2:1960][884:2:1248]')
%A = 
%           dsp_band_1: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.dsp_band_1
%ans = 
%            ml__unsigned: '"true"'
%      ml__CoordinateAxes: '"time lat lon "'
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
%            scale_factor: 0.1250
%                 add_off: -4
%       DODS_ML_Real_Name: 'dsp_band_1'
%              dsp_band_1: [1x1 struct]
%                    time: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
             
 %examples of URLs for SST (Declouded) and SST (Raw)
%http://satdat1.gso.uri.edu/thredds/dodsC/NEPacificDec_1km?dsp_band_1[4008][1035:2:1497][1131:2:1518]
%http://satdat1.gso.uri.edu/thredds/dodsC/NEPacificRaw_1km?dsp_band_1[4009][1035:2:1497][5013:2:5400]


%A=loaddap('-A','http://satdat1.gso.uri.edu/thredds/dodsC/NEPacificDec_1km?dsp_band_1[4008][1035:2:1497][1131:2:1518]')
%A = 
%           dsp_band_1: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.dsp_band_1
%ans = 
%            ml__unsigned: '"true"'
%      ml__CoordinateAxes: '"time lat lon "'
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
%            scale_factor: 0.1250
%                 add_off: -4
%       DODS_ML_Real_Name: 'dsp_band_1'
%              dsp_band_1: [1x1 struct]
%                    time: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
%A=loaddap('-A','http://satdat1.gso.uri.edu/thredds/dodsC/NEPacificRaw_1km?dsp_band_1[4009][1035:2:1497][5013:2:5400]')
%A = 
%           dsp_band_1: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.dsp_band_1
%ans = 
%            ml__unsigned: '"true"'
%      ml__CoordinateAxes: '"time lat lon "'
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
%            scale_factor: 0.1250
%                 add_off: -4
%       DODS_ML_Real_Name: 'dsp_band_1'
%              dsp_band_1: [1x1 struct]
%                    time: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%examples of URLs for SST (Declouded) and SST (Raw)
%http://satdat1.gso.uri.edu/thredds/dodsC/NWPacificDec_1km?dsp_band_1[12347][1313:2:1775][0:2:6143]
%http://satdat1.gso.uri.edu/thredds/dodsC/NWPacificRaw_1km?dsp_band_1[12474][1313:2:1775][0:2:6143]

%A=loaddap('-A','http://satdat1.gso.uri.edu/thredds/dodsC/NWPacificDec_1km?dsp_band_1[12347][1313:2:1775][0:2:6143]')
%A = 
%           dsp_band_1: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.dsp_band_1
%ans = 
%            ml__unsigned: '"true"'
%      ml__CoordinateAxes: '"time lat lon "'
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
%            scale_factor: 0.1250
%                 add_off: -4
%       DODS_ML_Real_Name: 'dsp_band_1'
%              dsp_band_1: [1x1 struct]
%                    time: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]

%A=loaddap('-A','http://satdat1.gso.uri.edu/thredds/dodsC/NWPacificRaw_1km?dsp_band_1[12474][1313:2:1775][0:2:6143]')
%A = 
%           dsp_band_1: [1x1 struct]
%    Global_Attributes: [1x1 struct]
%A.dsp_band_1
%ans = 
%            ml__unsigned: '"true"'
%      ml__CoordinateAxes: '"time lat lon "'
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
%            scale_factor: 0.1250
%                 add_off: -4
%       DODS_ML_Real_Name: 'dsp_band_1'
%              dsp_band_1: [1x1 struct]
%                    time: [1x1 struct]
%                     lat: [1x1 struct]
%                     lon: [1x1 struct]
                     
                     
                     
                     