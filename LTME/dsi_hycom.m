function dsi_hycom(varargin)
% updates the DataSetInventory structures
% without arguments for all the following DataSetBranches
% url location of the dataset files

if nargin < 1

List.DataSetBranch={
%'HYCOM_Global_Near_Real_Time_Outputs';
'HYCOM_Global_Hindcast';
'HYCOM_Global_Simulation';
'HYCOM_Atlantic_Ocean_Prediction_System';    
};

    for k=1:length(List.DataSetBranch)
    disp(['DSI_' List.DataSetBranch{k} ])
    eval(['DSI_' List.DataSetBranch{k} ';'])
    end
    
dsi_filename ='dsi_hycom.mat';
%if (~isdeployed)
pname = strrep(mfilename('fullpath'),mfilename,'');
dsi_filename = [pname,dsi_filename];
%end
if ~exist(dsi_filename,'file')
save(dsi_filename,'List')
end
save(dsi_filename,'List','-append')

elseif nargin > 1
disp('Incorrect number of arguments.'); return
else
% with one argument for the specified DataSetBranch only
eval(varargin{1});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DSI=DSI_HYCOM_Global_Near_Real_Time_Outputs
% create DataSetInventory of HYCOM
DSI.DataSetName='HYCOM';
% Example of loaddap argument
% http://hycom.coaps.fsu.edu/thredds/dodsC/glb_nrt_analysis?temperature[0:1:1][0:1:2][0:1:3][0:1:4]
DSI.DataSetBranch='HYCOM_Global_Near_Real_Time_Outputs';
% http://hycom.coaps.fsu.edu/thredds/dodsC/glb_nrt_analysis?temperature[0:1:1][0:1:2][0:1:3][0:1:4]
DSI.URLSITE='http://hycom.coaps.fsu.edu';
DSI.URLPATH='/thredds/dodsC/';
DSI.URLFILE='glb_nrt_analysis';
URL=[DSI.URLSITE DSI.URLPATH DSI.URLFILE];
loaddap('+v',[URL '?MT'])
t0=datenum('1900-12-31 00:00:00','yyyy-mm-dd HH:MM:SS');
t=t0+MT;
DSI.TimeA=t;
DSI.TimeB=t+1;
DSI.Time=t+0.5;
DSI.NY=3298;DSI.NX=4500;
loaddap('+v',[URL '?Latitude[0:1:3297][0]'])
DSI.Latitude=Latitude.Latitude;
loaddap('+v',[URL '?Longitude[0][0:1:4499]'])
DSI.Longitude=Longitude.Longitude';
loaddap('+v',[URL '?Depth'])
DSI.Depth=Depth;

% A=loaddap('-A','http://hycom.coaps.fsu.edu/thredds/dodsC/glb_nrt_analysis')
%A = 
%                                   Y: [1x1 struct]
%                                   X: [1x1 struct]
%                            Latitude: [1x1 struct]
%                           Longitude: [1x1 struct]
%                                  MT: [1x1 struct]
%                                Date: [1x1 struct]
%                                qtot: [1x1 struct]
%                                 emp: [1x1 struct]
%           surface_temperature_trend: [1x1 struct]
%              surface_salinity_trend: [1x1 struct]
%                                 ssh: [1x1 struct]
%    surface_boundary_layer_thickness: [1x1 struct]
%               mixed_layer_thickness: [1x1 struct]
%                               Depth: [1x1 struct]
%                                   u: [1x1 struct]
%                                   v: [1x1 struct]
%                         temperature: [1x1 struct]
%                            salinity: [1x1 struct]
%                   Global_Attributes: [1x1 struct]

% define list of fields or variables available in the dataset
LL={'Y','X'};LLD={'Y','X','depth'};TYX={'MT','Y','X'};TDYX={'MT','Depth','Y','X'};
FieldsTable={
%    long_name                            dataset/opendap                    menu                               opendap   dataset
%    name                               , name                             , name                              ,dim,coord,dim,coord,CBNum,units
    'downward_heat_flux_in_air'         ,'qtot'                            ,'downward_heat_flux_in_air (qtot)' ,2  ,LL   ,3  ,TYX ,16,'wm-2';
    'water_flux_into_ocean'             ,'emp'                             ,'water_flux_into_ocean (emp)'      ,2  ,LL   ,3  ,TYX ,17,'kgm-2s-1';
    'surface_temperature_trend'         ,'surface_temperature_trend'       ,'temperature_trend'                ,2  ,LL   ,3  ,TYX ,13,'degC/day';
    'surface_salinity_trend'            ,'surface_salinity_trend'          ,'salinity_trend'                   ,2  ,LL   ,3  ,TYX ,12,'psu/day';
    'sea_surface_elevation'             ,'ssh'                             ,'sea_surface_elevation (ssh)'      ,2  ,LL   ,3  ,TYX ,14,'m';
%    'mixed_layer_u_velocity'            ,'mixed_layer_u_velocity'          ,'u_velocity'                      ,2  ,LL   ,3  ,TYX ,10,'ms-1';
%    'mixed_layer_v_velocity'            ,'mixed_layer_v_velocity'          ,'v_velocity'                      ,2  ,LL   ,3  ,TYX ,11,'ms-1';
    'surface_boundary_layer_thickness'  ,'surface_boundary_layer_thickness','boundary_layer_thickness'         ,2  ,LL   ,3  ,TYX ,15,'m';
    'mixed_layer_thickness'             ,'mixed_layer_thickness'           ,'thickness'                        ,2  ,LL   ,3  ,TYX , 9,'m';
%    'mixed_layer_temperature'           ,'mixed_layer_temperature'         ,'temperature'                      ,2  ,LL   ,3  ,TYX , 8,'degC';
%    'mixed_layer_salinity'              ,'mixed_layer_salinity'            ,'salinity'                         ,2  ,LL   ,3  ,TYX , 7,'psu';
%    'mixed_layer_density'               ,'mixed_layer_density'             ,'density'                          ,2  ,LL   ,3  ,TYX , 6,'sigma';
    'eastward_velocity'                 ,'u'                               ,'eastward_velocity  (u)'           ,3  ,LLD  ,4  ,TDYX, 4,'ms-1';
    'northward_velocity'                ,'v'                               ,'northward_velocity (v)'           ,3  ,LLD  ,4  ,TDYX, 5,'ms-1';
    'temperature'                       ,'temperature'                     ,'temperature'                      ,3  ,LLD  ,4  ,TDYX, 3,'degC';
    'salinity'                          ,'salinity'                        ,'salinity'                         ,3  ,LLD  ,4  ,TDYX, 2,'psu';
%    'density'                           ,'density'                         ,'density'                          ,3  ,LLD  ,4  ,TDYX, 1,'sigma';
    %'latitudeYX'                       ,'Latitude'                        ,'Latitude(Y,X)'                    ,2  ,LL   ,3  ,TYX ,91,'degN';
    %'longitudeYX'                      ,'Longitude'                       ,'Longitude(Y,X)'                   ,2  ,LL   ,3  ,TYX ,92,'degE';
     };

 DSI=DSI_common2(FieldsTable,DSI); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DSI=DSI_HYCOM_Global_Hindcast
% create DataSetInventory of HYCOM
DSI.DataSetName='HYCOM';
% Example of loaddap argument
% http://hycom.coaps.fsu.edu/thredds/dodsC/glb_analysis?temperature[0:1:1][0:1:2][0:1:3][0:1:4]
DSI.DataSetBranch='HYCOM_Global_Hindcast';
% http://hycom.coaps.fsu.edu/thredds/dodsC/glb_analysis?temperature[0:1:1][0:1:2][0:1:3][0:1:4]
DSI.URLSITE='http://hycom.coaps.fsu.edu';
DSI.URLPATH='/thredds/dodsC/';
DSI.URLFILE='glb_analysis';
URL=[DSI.URLSITE DSI.URLPATH DSI.URLFILE];
loaddap('+v',[URL '?MT'])
t0=datenum('1900-12-31 00:00:00','yyyy-mm-dd HH:MM:SS');
t=t0+MT;
DSI.TimeA=t;
DSI.TimeB=t+1;
DSI.Time=t+0.5;
DSI.NY=3298;DSI.NX=4500;
loaddap('+v',[URL '?Latitude[0:1:3297][0]'])
DSI.Latitude=Latitude.Latitude;
loaddap('+v',[URL '?Longitude[0][0:1:4499]'])
DSI.Longitude=Longitude.Longitude';
loaddap('+v',[URL '?Depth'])
DSI.Depth=Depth;

% A=loaddap('-A','http://hycom.coaps.fsu.edu/thredds/dodsC/glb_analysis')
%A = 
%                                   Y: [1x1 struct]
%                                   X: [1x1 struct]
%                            Latitude: [1x1 struct]
%                           Longitude: [1x1 struct]
%                                  MT: [1x1 struct]
%                                Date: [1x1 struct]
%                                qtot: [1x1 struct]
%                                 emp: [1x1 struct]
%           surface_temperature_trend: [1x1 struct]
%              surface_salinity_trend: [1x1 struct]
%                                 ssh: [1x1 struct]
%              mixed_layer_u_velocity: [1x1 struct]
%              mixed_layer_v_velocity: [1x1 struct]
%    surface_boundary_layer_thickness: [1x1 struct]
%               mixed_layer_thickness: [1x1 struct]
%             mixed_layer_temperature: [1x1 struct]
%                mixed_layer_salinity: [1x1 struct]
%                 mixed_layer_density: [1x1 struct]
%                               Depth: [1x1 struct]
%                                   u: [1x1 struct]
%                                   v: [1x1 struct]
%                         temperature: [1x1 struct]
%                            salinity: [1x1 struct]
%                   Global_Attributes: [1x1 struct]

% define list of fields or variables available in the dataset
LL={'Y','X'};LLD={'Y','X','depth'};TYX={'MT','Y','X'};TDYX={'MT','Depth','Y','X'};
FieldsTable={
%    long_name                            dataset/opendap                    menu                               opendap   dataset
%    name                               , name                             , name                              ,dim,coord,dim,coord,CBNum,units
    'downward_heat_flux_in_air'         ,'qtot'                            ,'downward_heat_flux_in_air (qtot)' ,2  ,LL   ,3  ,TYX ,16,'wm-2';
    'water_flux_into_ocean'             ,'emp'                             ,'water_flux_into_ocean (emp)'      ,2  ,LL   ,3  ,TYX ,17,'kgm-2s-1';
    'surface_temperature_trend'         ,'surface_temperature_trend'       ,'temperature_trend'                ,2  ,LL   ,3  ,TYX ,13,'degC/day';
    'surface_salinity_trend'            ,'surface_salinity_trend'          ,'salinity_trend'                   ,2  ,LL   ,3  ,TYX ,12,'psu/day';
    'sea_surface_elevation'             ,'ssh'                             ,'sea_surface_elevation (ssh)'      ,2  ,LL   ,3  ,TYX ,14,'m';
    'mixed_layer_u_velocity'            ,'mixed_layer_u_velocity'          ,'u_velocity'                       ,2  ,LL   ,3  ,TYX ,10,'ms-1';
    'mixed_layer_v_velocity'            ,'mixed_layer_v_velocity'          ,'v_velocity'                       ,2  ,LL   ,3  ,TYX ,11,'ms-1';
    'surface_boundary_layer_thickness'  ,'surface_boundary_layer_thickness','boundary_layer_thickness'         ,2  ,LL   ,3  ,TYX ,15,'m';
    'mixed_layer_thickness'             ,'mixed_layer_thickness'           ,'thickness'                        ,2  ,LL   ,3  ,TYX , 9,'m';
    'mixed_layer_temperature'           ,'mixed_layer_temperature'         ,'temperature'                      ,2  ,LL   ,3  ,TYX , 8,'degC';
    'mixed_layer_salinity'              ,'mixed_layer_salinity'            ,'salinity'                         ,2  ,LL   ,3  ,TYX , 7,'psu';
    'mixed_layer_density'               ,'mixed_layer_density'             ,'density'                          ,2  ,LL   ,3  ,TYX , 6,'sigma';
    'eastward_velocity'                 ,'u'                               ,'eastward_velocity  (u)'           ,3  ,LLD  ,4  ,TDYX, 4,'ms-1';
    'northward_velocity'                ,'v'                               ,'northward_velocity (v)'           ,3  ,LLD  ,4  ,TDYX, 5,'ms-1';
    'temperature'                       ,'temperature'                     ,'temperature'                      ,3  ,LLD  ,4  ,TDYX, 3,'degC';
    'salinity'                          ,'salinity'                        ,'salinity'                         ,3  ,LLD  ,4  ,TDYX, 2,'psu';
%    'density'                           ,'density'                         ,'density'                          ,3  ,LLD  ,4  ,TDYX, 1,'sigma';
    %'latitudeYX'                       ,'Latitude'                        ,'Latitude(Y,X)'                    ,2  ,LL   ,3  ,TYX ,91,'degN';
    %'longitudeYX'                      ,'Longitude'                       ,'Longitude(Y,X)'                   ,2  ,LL   ,3  ,TYX ,92,'degE';
     };

 DSI=DSI_common2(FieldsTable,DSI);  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DSI=DSI_HYCOM_Global_Simulation
% create DataSetInventory of HYCOM
DSI.DataSetName='HYCOM';
% Example of loaddap argument
% http://hycom.coaps.fsu.edu/thredds/dodsC/glb_analysis?temperature[0:1:1][0:1:2][0:1:3][0:1:4]
DSI.DataSetBranch='HYCOM_Global_Simulation';
% http://hycom.coaps.fsu.edu/thredds/dodsC/glb_analysis?temperature[0:1:1][0:1:2][0:1:3][0:1:4]
DSI.URLSITE='http://hycom.coaps.fsu.edu';
DSI.URLPATH='/thredds/dodsC/';
DSI.URLFILE='glb_simulation';
URL=[DSI.URLSITE DSI.URLPATH DSI.URLFILE];
loaddap('+v',[URL '?MT'])
t0=datenum('1900-12-31 00:00:00','yyyy-mm-dd HH:MM:SS');
t=t0+MT;
DSI.TimeA=t;
DSI.TimeB=t+1;
DSI.Time=t+0.5;
DSI.NY=3298;DSI.NX=4500;
loaddap('+v',[URL '?Latitude[0:1:3297][0]'])
DSI.Latitude=Latitude.Latitude;
loaddap('+v',[URL '?Longitude[0][0:1:4499]'])
DSI.Longitude=Longitude.Longitude';
loaddap('+v',[URL '?Depth'])
DSI.Depth=Depth;

%A=loaddap('-A','http://hycom.coaps.fsu.edu/thredds/dodsC/glb_simulation')
%A = 
%                                   Y: [1x1 struct]
%                                   X: [1x1 struct]
%                            Latitude: [1x1 struct]
%                           Longitude: [1x1 struct]
%                                  MT: [1x1 struct]
%                                Date: [1x1 struct]
%                                qtot: [1x1 struct]
%                                 emp: [1x1 struct]
%           surface_temperature_trend: [1x1 struct]
%             surface_salinity_trend: [1x1 struct]
%                                 ssh: [1x1 struct]
%    surface_boundary_layer_thickness: [1x1 struct]
%               mixed_layer_thickness: [1x1 struct]
%                               Depth: [1x1 struct]
%                                   u: [1x1 struct]
%                                   v: [1x1 struct]
%                         temperature: [1x1 struct]
%                            salinity: [1x1 struct]
%                             density: [1x1 struct]
%                   Global_Attributes: [1x1 struct]

% define list of fields or variables available in the dataset
LL={'Y','X'};LLD={'Y','X','depth'};TYX={'MT','Y','X'};TDYX={'MT','Depth','Y','X'};
FieldsTable={
%    long_name                            dataset/opendap                    menu                               opendap   dataset
%    name                               , name                             , name                              ,dim,coord,dim,coord,CBNum,units
    'downward_heat_flux_in_air'         ,'qtot'                            ,'downward_heat_flux_in_air (qtot)' ,2  ,LL   ,3  ,TYX ,16,'wm-2';
    'water_flux_into_ocean'             ,'emp'                             ,'water_flux_into_ocean (emp)'      ,2  ,LL   ,3  ,TYX ,17,'kgm-2s-1';
    'surface_temperature_trend'         ,'surface_temperature_trend'       ,'temperature_trend'                ,2  ,LL   ,3  ,TYX ,13,'degC/day';
    'surface_salinity_trend'            ,'surface_salinity_trend'          ,'salinity_trend'                   ,2  ,LL   ,3  ,TYX ,12,'psu/day';
    'sea_surface_elevation'             ,'ssh'                             ,'sea_surface_elevation (ssh)'      ,2  ,LL   ,3  ,TYX ,14,'m';
%    'mixed_layer_u_velocity'            ,'mixed_layer_u_velocity'          ,'u_velocity'                      ,2  ,LL   ,3  ,TYX ,10,'ms-1';
%    'mixed_layer_v_velocity'            ,'mixed_layer_v_velocity'          ,'v_velocity'                      ,2  ,LL   ,3  ,TYX ,11,'ms-1';
    'surface_boundary_layer_thickness'  ,'surface_boundary_layer_thickness','boundary_layer_thickness'         ,2  ,LL   ,3  ,TYX ,15,'m';
    'mixed_layer_thickness'             ,'mixed_layer_thickness'           ,'thickness'                        ,2  ,LL   ,3  ,TYX , 9,'m';
%    'mixed_layer_temperature'           ,'mixed_layer_temperature'         ,'temperature'                      ,2  ,LL   ,3  ,TYX , 8,'degC';
%    'mixed_layer_salinity'              ,'mixed_layer_salinity'            ,'salinity'                         ,2  ,LL   ,3  ,TYX , 7,'psu';
%    'mixed_layer_density'               ,'mixed_layer_density'             ,'density'                          ,2  ,LL   ,3  ,TYX , 6,'sigma';
    'eastward_velocity'                 ,'u'                               ,'eastward_velocity  (u)'           ,3  ,LLD  ,4  ,TDYX, 4,'ms-1';
    'northward_velocity'                ,'v'                               ,'northward_velocity (v)'           ,3  ,LLD  ,4  ,TDYX, 5,'ms-1';
    'temperature'                       ,'temperature'                     ,'temperature'                      ,3  ,LLD  ,4  ,TDYX, 3,'degC';
    'salinity'                          ,'salinity'                        ,'salinity'                         ,3  ,LLD  ,4  ,TDYX, 2,'psu';
    'density'                           ,'density'                         ,'density'                          ,3  ,LLD  ,4  ,TDYX, 1,'sigma';
    %'latitudeYX'                       ,'Latitude'                        ,'Latitude(Y,X)'                    ,2  ,LL   ,3  ,TYX ,91,'degN';
    %'longitudeYX'                      ,'Longitude'                       ,'Longitude(Y,X)'                   ,2  ,LL   ,3  ,TYX ,92,'degE';
     };

 DSI=DSI_common2(FieldsTable,DSI); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DSI=DSI_HYCOM_Atlantic_Ocean_Prediction_System
% create DataSetInventory of HYCOM
DSI.DataSetName='HYCOM';
% Example of loaddap argument
% http://hycom.coaps.fsu.edu/thredds/dodsC/glb_analysis?temperature[0:1:1][0:1:2][0:1:3][0:1:4]
DSI.DataSetBranch='HYCOM_Atlantic_Ocean_Prediction_System';
% http://hycom.coaps.fsu.edu/thredds/dodsC/glb_analysis?temperature[0:1:1][0:1:2][0:1:3][0:1:4]
DSI.URLSITE='http://hycom.coaps.fsu.edu';
DSI.URLPATH='/thredds/dodsC/';
DSI.URLFILE='atl_ops';
URL=[DSI.URLSITE DSI.URLPATH DSI.URLFILE];
DSI.NY=3298;DSI.NX=4500;
% in this DataSet Latitude/Longitude is 1D array
loaddap('+v',[URL '?Latitude'])
DSI.Latitude=Latitude;
loaddap('+v',[URL '?Longitude'])
DSI.Longitude=Longitude;
loaddap('+v',[URL '?Depth'])
DSI.Depth=Depth;
loaddap('+v',[URL '?MT'])
t0=datenum('1900-12-31 00:00:00','yyyy-mm-dd HH:MM:SS');
t=t0+MT;
DSI.TimeA=t;
DSI.TimeB=t+1;
DSI.Time=t+0.5;

%A=loaddap('-A','http://hycom.coaps.fsu.edu/thredds/dodsC/atl_ops')
%A = 
%                            Latitude: [1x1 struct]
%                           Longitude: [1x1 struct]
%                                  MT: [1x1 struct]
%                                Date: [1x1 struct]
%                                qtot: [1x1 struct]
%                                 emp: [1x1 struct]
%           surface_temperature_trend: [1x1 struct]
%              surface_salinity_trend: [1x1 struct]
%                                 ssh: [1x1 struct]
%              mixed_layer_u_velocity: [1x1 struct]
%              mixed_layer_v_velocity: [1x1 struct]
%    surface_boundary_layer_thickness: [1x1 struct]
%               mixed_layer_thickness: [1x1 struct]
%             mixed_layer_temperature: [1x1 struct]
%                mixed_layer_salinity: [1x1 struct]
%                 mixed_layer_density: [1x1 struct]
%                               Depth: [1x1 struct]
%                                   u: [1x1 struct]
%                                   v: [1x1 struct]
%                         temperature: [1x1 struct]
%                            salinity: [1x1 struct]
%                             density: [1x1 struct]
%                   Global_Attributes: [1x1 struct]

% define list of fields or variables available in the dataset
LL={'Y','X'};LLD={'Y','X','depth'};TYX={'MT','Y','X'};TDYX={'MT','Depth','Y','X'};
FieldsTable={
%    long_name                            dataset/opendap                    menu                               opendap   dataset
%    name                               , name                             , name                              ,dim,coord,dim,coord,CBNum,units
    'downward_heat_flux_in_air'         ,'qtot'                            ,'downward_heat_flux_in_air (qtot)' ,2  ,LL   ,3  ,TYX ,16,'wm-2';
    'water_flux_into_ocean'             ,'emp'                             ,'water_flux_into_ocean (emp)'      ,2  ,LL   ,3  ,TYX ,17,'kgm-2s-1';
    'surface_temperature_trend'         ,'surface_temperature_trend'       ,'temperature_trend'                ,2  ,LL   ,3  ,TYX ,13,'degC/day';
    'surface_salinity_trend'            ,'surface_salinity_trend'          ,'salinity_trend'                   ,2  ,LL   ,3  ,TYX ,12,'psu/day';
    'sea_surface_elevation'             ,'ssh'                             ,'sea_surface_elevation (ssh)'      ,2  ,LL   ,3  ,TYX ,14,'m';
    'mixed_layer_u_velocity'            ,'mixed_layer_u_velocity'          ,'u_velocity'                       ,2  ,LL   ,3  ,TYX ,10,'ms-1';
    'mixed_layer_v_velocity'            ,'mixed_layer_v_velocity'          ,'v_velocity'                       ,2  ,LL   ,3  ,TYX ,11,'ms-1';
    'surface_boundary_layer_thickness'  ,'surface_boundary_layer_thickness','boundary_layer_thickness'         ,2  ,LL   ,3  ,TYX ,15,'m';
    'mixed_layer_thickness'             ,'mixed_layer_thickness'           ,'thickness'                        ,2  ,LL   ,3  ,TYX , 9,'m';
    'mixed_layer_temperature'           ,'mixed_layer_temperature'         ,'temperature'                      ,2  ,LL   ,3  ,TYX , 8,'degC';
    'mixed_layer_salinity'              ,'mixed_layer_salinity'            ,'salinity'                         ,2  ,LL   ,3  ,TYX , 7,'psu';
    'mixed_layer_density'               ,'mixed_layer_density'             ,'density'                          ,2  ,LL   ,3  ,TYX , 6,'sigma';
    'eastward_velocity'                 ,'u'                               ,'eastward_velocity  (u)'           ,3  ,LLD  ,4  ,TDYX, 4,'ms-1';
    'northward_velocity'                ,'v'                               ,'northward_velocity (v)'           ,3  ,LLD  ,4  ,TDYX, 5,'ms-1';
    'temperature'                       ,'temperature'                     ,'temperature'                      ,3  ,LLD  ,4  ,TDYX, 3,'degC';
    'salinity'                          ,'salinity'                        ,'salinity'                         ,3  ,LLD  ,4  ,TDYX, 2,'psu';
    'density'                           ,'density'                         ,'density'                          ,3  ,LLD  ,4  ,TDYX, 1,'sigma';
    %'latitudeYX'                       ,'Latitude'                        ,'Latitude(Y,X)'                    ,2  ,LL   ,3  ,TYX ,91,'degN';
    %'longitudeYX'                      ,'Longitude'                       ,'Longitude(Y,X)'                   ,2  ,LL   ,3  ,TYX ,92,'degE';
     };

 DSI=DSI_common2(FieldsTable,DSI); 

function DSI=DSI_common2(FieldsTable,DSI); 
DSI.Fields         =FieldsTable(:,2); 
DSI.Fields_NameDS  =FieldsTable(:,2);
DSI.Fields_NameMenu=FieldsTable(:,3);
DSI.Fields_Dim     =FieldsTable(:,4);
DSI.Fields_Coord   =FieldsTable(:,5);
DSI.Fields_DimDS   =FieldsTable(:,6);
DSI.Fields_CoordDS =FieldsTable(:,7);
DSI.Fields_CBNum   =FieldsTable(:,8);

CoordinatesTable={
'latitude','latitude','degrees_north';
'longitude','longitude','degrees_east';
'depth','depth','m';
'time','time','days since 0000-01-01 00:00:00';
};
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
eval(['DSI.URLCVAR.' FIELD '=FieldsTable{k,2};']);
end

DSI.Coordinates=CoordinatesTable(:,1);
DSI.Dimensions.Y=DSI.NY;
DSI.Dimensions.X=DSI.NX;
DSI.Dimensions.depth=length(DSI.Depth);
DSI.Dimensions.time=length(DSI.Time);
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
eval(['DSI.Variables.' FIELD '=FieldsTable{k,5};'])
end
DSI.Variables.latitude={'latitude'};
DSI.Variables.longitude={'longitude'};
DSI.Variables.depth={'depth'};
DSI.Variables.time={'time'};

for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
eval(['DSI.Attributes.' FIELD '.long_name=FieldsTable{k,1};']);
eval(['DSI.Attributes.' FIELD '.units    =FieldsTable{k,9};']);
end
for k=1:length(DSI.Coordinates)
    COORD=DSI.Coordinates{k};
eval(['DSI.Attributes.' COORD '.long_name=CoordinatesTable{k,2};']);
eval(['DSI.Attributes.' COORD '.units    =CoordinatesTable{k,3};']);
end

% Rename returned variables and change dimensions order
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
    FREN=['f=squeeze(' FIELD '.' FIELD ');'];
eval(['DSI.FormulaRen.' FIELD '=''' FREN ''';']);      
end
% ml__FillValue: 1.2677e+030 for all fields
DSI.FormulaNaN='f(find(f>1.e30))=NaN;';
DSI.FormulaCnv='';

DSI.Range.Time1=min(DSI.TimeA);
DSI.Range.Time2=max(DSI.TimeB);
DSI.Range.Latitude1=min(DSI.Latitude);
DSI.Range.Latitude2=max(DSI.Latitude);
DSI.Range.Longitude1=min(DSI.Longitude);
DSI.Range.Longitude2=max(DSI.Longitude);

DSI.Attributes.latitude.units='degrees_north';
DSI.Attributes.latitude.valid_range=[min(DSI.Latitude) max(DSI.Latitude)];
DSI.Attributes.longitude.units='degrees_east';
DSI.Attributes.longitude.valid_range=[min(DSI.Longitude) max(DSI.Longitude)];
DSI.Attributes.depth.units='m';
DSI.Attributes.depth.valid_range=[min(DSI.Depth) max(DSI.Depth)];

S1=['time in seconds since 1970-01-01 00:00:00'];
S2=['missing values replaced with NaN'];
%S3=['variables converted to physical units;'];
DSI.Readme=[S1 '/' S2 ];

dsi_filename ='dsi_hycom.mat';
%if (~isdeployed)
pname = strrep(mfilename('fullpath'),mfilename,'');
dsi_filename = [pname,dsi_filename];
%end

eval(['DSI_' DSI.DataSetBranch '=DSI;'])
if ~exist(dsi_filename,'file')
save(dsi_filename,['DSI_' DSI.DataSetBranch])
end
save(dsi_filename,['DSI_' DSI.DataSetBranch],'-append')


 
