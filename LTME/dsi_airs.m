function dsi_airs(varargin)
% updates the DataSetInventory structures
% without arguments for all the following DataSetBranches

if nargin < 1
%DSI_AIRS
DataSetBranchList={
'AIRS_daily';
'AIRS_weekly';
'AIRS_monthly';
};

dsi_filename ='dsi_airs.mat';
%if (~isdeployed)
pname = strrep(mfilename('fullpath'),mfilename,'');
dsi_filename = [pname,dsi_filename];
%end

if ~exist(dsi_filename,'file')
save(dsi_filename,'DataSetBranchList')
end
save(dsi_filename,'DataSetBranchList','-append')

for k=1:length(DataSetBranchList)    
eval(['DSI=DSI_' DataSetBranchList{k} ';'])
eval(['DSI_' DSI.DataSetBranch '=DSI;'])
if ~exist(dsi_filename,'file')
save(dsi_filename,['DSI_' DSI.DataSetBranch])
end
save(dsi_filename,['DSI_' DSI.DataSetBranch],'-append')
end

elseif nargin > 1
disp('Incorrect number of arguments.'); return
else
% with one argument for the specified DataSetBranch only
eval(varargin{1});
eval(['DSI_' DSI.DataSetBranch '=DSI;'])
if ~exist(dsi_filename,'file')
save(dsi_filename,['DSI_' DSI.DataSetBranch])
end
save(dsi_filename,['DSI_' DSI.DataSetBranch],'-append')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DSI=DSI_AIRS_daily
% create DataSetInventory of AIRS

DSI.DataSetName='AIRS';
DSI.DataSetBranch='AIRS_daily';
DSI.TAVG='daily';

% a general request employed by get_DataSet program has the following pattern
% loaddap([URLSITE URLPATH URLFILE '?' URLCVAR URLCTIME URLCLAT URLCLON])
%   constraints URLCTIME URLCLAT URLCLON can be in different order
%   constraints URLCTIME URLCLAT URLCLON will be generated by program get

% Example of loaddap argument
% daily
%http://acdisc.sci.gsfc.nasa.gov/opendap/Aqua_AIRS_Level3/AIRX3STD.005/2007/AIRS.2007.01.01.L3.RetStd001.v5.0.14.0.G07348100416.hdf?
%TotCldLiqH2O_D[88:1:89][135:1:136],TotCldLiqH2O_D_sdev[88:1:89][135:1:136],TotCldLiqH2O_D_ct[88:1:89][135:1:136],
%TotH2OVap_D[88:1:89][135:1:136],TotH2OVap_D_sdev[88:1:89][135:1:136],TotH2OVap_D_ct[88:1:89][135:1:136],
%TotH2OVap_MW_D[88:1:89][135:1:136],TotH2OVap_MW_D_sdev[88:1:89][135:1:136],TotH2OVap_MW_D_ct[88:1:89][135:1:136],
%SurfAirTemp_D[88:1:89][135:1:136],SurfAirTemp_D_sdev[88:1:89][135:1:136],SurfAirTemp_D_ct[88:1:89][135:1:136],
%SurfSkinTemp_D[88:1:89][135:1:136],SurfSkinTemp_D_sdev[88:1:89][135:1:136],SurfSkinTemp_D_ct[88:1:89][135:1:136],
%SurfPres_D[88:1:89][135:1:136],SurfPres_D_sdev[88:1:89][135:1:136],SurfPres_D_ct[88:1:89][135:1:136],
%TotO3_D[88:1:89][135:1:136],TotO3_D_sdev[88:1:89][135:1:136],TotO3_D_ct[88:1:89][135:1:136],
%CloudFrc_D[88:1:89][135:1:136],CloudFrc_D_sdev[88:1:89][135:1:136],CloudFrc_D_ct[88:1:89][135:1:136],
%CloudTopPres_D[88:1:89][135:1:136],CloudTopPres_D_sdev[88:1:89][135:1:136],CloudTopPres_D_ct[88:1:89][135:1:136],
%CloudTopTemp_D[88:1:89][135:1:136],CloudTopTemp_D_sdev[88:1:89][135:1:136],CloudTopTemp_D_ct[88:1:89][135:1:136],
%OLR_D[88:1:89][135:1:136],OLR_D_sdev[88:1:89][135:1:136],OLR_D_ct[88:1:89][135:1:136],
%ClrOLR_D[88:1:89][135:1:136],ClrOLR_D_sdev[88:1:89][135:1:136],ClrOLR_D_ct[88:1:89][135:1:136],
%Temperature_D[0:1:2][88:1:89][135:1:136],Temperature_D_sdev[0:1:2][88:1:89][135:1:136],Temperature_D_ct[0:1:2][88:1:89][135:1:136],
%Temperature_MW_D[0:1:2][88:1:89][135:1:136],Temperature_MW_D_sdev[0:1:2][88:1:89][135:1:136],Temperature_MW_D_ct[0:1:2][88:1:89][135:1:136],
%GPHeight_D[0:1:2][88:1:89][135:1:136],GPHeight_D_sdev[0:1:2][88:1:89][135:1:136],GPHeight_D_ct[0:1:2][88:1:89][135:1:136],
%GPHeight_MW_D[0:1:2][88:1:89][135:1:136],GPHeight_MW_D_sdev[0:1:2][88:1:89][135:1:136],GPHeight_MW_D_ct[0:1:2][88:1:89][135:1:136],
%RelHumid_D[0:1:2][88:1:89][135:1:136],RelHumid_D_sdev[0:1:2][88:1:89][135:1:136],RelHumid_D_ct[0:1:2][88:1:89][135:1:136],
%H2OVapMMR_D[0:1:2][88:1:89][135:1:136],H2OVapMMR_D_sdev[0:1:2][88:1:89][135:1:136],H2OVapMMR_D_ct[0:1:2][88:1:89][135:1:136],
%EmisIR_D[0:1:2][88:1:89][135:1:136],EmisIR_D_sdev[0:1:2][88:1:89][135:1:136],EmisIR_D_ct[0:1:2][88:1:89][135:1:136],
%EmisMW_MW_D[0:1:2][88:1:89][135:1:136],EmisMW_MW_D_sdev[0:1:2][88:1:89][135:1:136],EmisMW_MW_D_ct[0:1:2][88:1:89][135:1:136],
%LandSeaMask[88:1:89][135:1:136]
%http://acdisc.sci.gsfc.nasa.gov/opendap/Aqua_AIRS_Level3/AIRX3STD.005/2007/AIRS.2007.01.01.L3.RetStd001.v5.0.14.0.G07348100416.hdf?TotCldLiqH2O_A[88:1:89][135:1:136],TotCldLiqH2O_A_sdev[88:1:89][135:1:136],TotCldLiqH2O_A_ct[88:1:89][135:1:136],TotH2OVap_A[88:1:89][135:1:136],TotH2OVap_A_sdev[88:1:89][135:1:136],TotH2OVap_A_ct[88:1:89][135:1:136],TotH2OVap_MW_A[88:1:89][135:1:136],TotH2OVap_MW_A_sdev[88:1:89][135:1:136],TotH2OVap_MW_A_ct[88:1:89][135:1:136],SurfAirTemp_A[88:1:89][135:1:136],SurfAirTemp_A_sdev[88:1:89][135:1:136],SurfAirTemp_A_ct[88:1:89][135:1:136],SurfSkinTemp_A[88:1:89][135:1:136],SurfSkinTemp_A_sdev[88:1:89][135:1:136],SurfSkinTemp_A_ct[88:1:89][135:1:136],SurfPres_A[88:1:89][135:1:136],SurfPres_A_sdev[88:1:89][135:1:136],SurfPres_A_ct[88:1:89][135:1:136],TotO3_A[88:1:89][135:1:136],TotO3_A_sdev[88:1:89][135:1:136],TotO3_A_ct[88:1:89][135:1:136],CloudFrc_A[88:1:89][135:1:136],CloudFrc_A_sdev[88:1:89][135:1:136],CloudFrc_A_ct[88:1:89][135:1:136],CloudTopPres_A[88:1:89][135:1:136],CloudTopPres_A_sdev[88:1:89][135:1:136],CloudTopPres_A_ct[88:1:89][135:1:136],CloudTopTemp_A[88:1:89][135:1:136],CloudTopTemp_A_sdev[88:1:89][135:1:136],CloudTopTemp_A_ct[88:1:89][135:1:136],OLR_A[88:1:89][135:1:136],OLR_A_sdev[88:1:89][135:1:136],OLR_A_ct[88:1:89][135:1:136],ClrOLR_A[88:1:89][135:1:136],ClrOLR_A_sdev[88:1:89][135:1:136],ClrOLR_A_ct[88:1:89][135:1:136],Temperature_A[0:1:2][88:1:89][135:1:136],Temperature_A_sdev[0:1:2][88:1:89][135:1:136],Temperature_A_ct[0:1:2][88:1:89][135:1:136],Temperature_MW_A[0:1:2][88:1:89][135:1:136],Temperature_MW_A_sdev[0:1:2][88:1:89][135:1:136],Temperature_MW_A_ct[0:1:2][88:1:89][135:1:136],GPHeight_A[0:1:2][88:1:89][135:1:136],GPHeight_A_sdev[0:1:2][88:1:89][135:1:136],GPHeight_A_ct[0:1:2][88:1:89][135:1:136],GPHeight_MW_A[0:1:2][88:1:89][135:1:136],GPHeight_MW_A_sdev[0:1:2][88:1:89][135:1:136],GPHeight_MW_A_ct[0:1:2][88:1:89][135:1:136],RelHumid_A[0:1:2][88:1:89][135:1:136],RelHumid_A_sdev[0:1:2][88:1:89][135:1:136],RelHumid_A_ct[0:1:2][88:1:89][135:1:136],H2OVapMMR_A[0:1:2][88:1:89][135:1:136],H2OVapMMR_A_sdev[0:1:2][88:1:89][135:1:136],H2OVapMMR_A_ct[0:1:2][88:1:89][135:1:136],EmisIR_A[0:1:2][88:1:89][135:1:136],EmisIR_A_sdev[0:1:2][88:1:89][135:1:136],EmisIR_A_ct[0:1:2][88:1:89][135:1:136],EmisMW_MW_A[0:1:2][88:1:89][135:1:136],EmisMW_MW_A_sdev[0:1:2][88:1:89][135:1:136],EmisMW_MW_A_ct[0:1:2][88:1:89][135:1:136],LandSeaMask[88:1:89][135:1:136]
DSI.URLSITE='http://acdisc.sci.gsfc.nasa.gov/opendap/Aqua_AIRS_Level3/';
DSI.URLPATH='opendap/Aqua_AIRS_Level3';
DSI.URLFILE='/AIRX3STD.005/';
%URL=[DSI.URLSITE DSI.URLPATH DSI.URLFILE];

% define list of fields or variables available in the dataset
LL={'latitude','longitude'};
FieldsTable={
% opendap name,    ( dataset name,   menu name,)           long name,          units      
   'TotCldLiqH2O'       ,'Mean total integrated column cloud liquid water'     ,'kg m-2';
   'TotCldLiqH2O_sdev'  ,'TotCldLiqH2O_sdev'                                   ,'kg m-2';
   'TotCldLiqH2O_ct'    ,'TotCldLiqH2O_ct'                                     ,'';    
   'TotH2OVap'          ,'Total integrated column water vapor burden'          ,'kg m-2';
   'TotH2OVap_sdev'     ,'TotH2OVap_sdev'                                      ,'kg m-2';
   'TotH2OVap_ct'       ,'TotH2OVap_ct'                                        ,'';
   'TotH2OVap_MW'       ,'Total integrated column water vapor burden'          ,'kg m-2';   
   'TotH2OVap_MW_sdev'  ,'TotH2OVap_MW_sdev'                                   ,'kg m-2';
   'TotH2OVap_MW_ct'    ,'TotH2OVap_MW_ct'                                     ,'';
   'SurfAirTemp'        ,'Temperature of the atmosphere at the Earths surface' ,'K';
   'SurfAirTemp_sdev'   ,'SurfAirTemp_sdev'                                    ,'K';
   'SurfAirTemp_ct'     ,'SurfAirTemp_ct'                                      ,'';
   'SurfSkinTemp'       ,'Surface skin temperature'                            ,'K';
   'SurfSkinTemp_sdev'  ,'SurfSkinTemp_sdev'                                   ,'K'; 
   'SurfSkinTemp_ct'    ,'SurfSkinTemp_ct'                                     ,'';
   'SurfPres'           ,'Mean surface pressure'                               ,'mbar';
   'SurfPres_sdev'      ,'SurfPres_sdev'                                       ,'mbar';
   'SurfPres_ct'        ,'SurfPres_ct'                                         ,'';
   'TotO3'              ,'Total integrated column ozone burden'                ,'';
   'TotO3_sdev'         ,'TotO3_sdev'                                          ,'';
   'TotO3_ct'           ,'TotO3_ct'                                            ,'';
   'CloudFrc'           ,'Combined layer cloud fraction'                       ,'DU';
   'CloudFrc_sdev'      ,'CloudFrc_sdev'                                       ,'DU';
   'CloudFrc_ct'        ,'CloudFrc_sdev'                                       ,'';
   'CloudTopPres'       ,'Combined cloud top pressure'                         ,'mbar';   
   'CloudTopPres_sdev'  ,'CloudTopPres_sdev'                                   ,'mbar';
   'CloudTopPres_ct'    ,'CloudTopPres_ct'                                     ,'';
   'CloudTopTemp'       ,'Combined cloud top temperature'                      ,'K';
   'CloudTopTemp_sdev'  ,'CloudTopTemp_sdev'                                   ,'K';
   'CloudTopTemp_ct'    ,'CloudTopTemp_ct'                                     ,'';
   'OLR'                ,'Outgoing long-wave radiation flux'                   ,'Watts m-2';
   'OLR_sdev'           ,'OLR_sdev'                                            ,'Watts m-2';
   'OLR_ct'             ,'OLR_ct'                                              ,'';
   'ClrOLR'             ,'Clear-sky outgoing long-wave radiation flux'         ,'Watts m-2';
   'ClrOLR_sdev'        ,'ClrOLR_sdev'                                         ,'Watts m-2';
   'ClrOLR_ct'          ,'ClrOLR_ct'                                           ,'';
   'Temperature'        ,'Atmospheric temperature profile in 24 standard pressure levels from 1000 to 1.0 mbar' ,'K';
   'Temperature_sdev'   ,'Temperature_sdev'                                                                     ,'K';
   'Temperature_ct'     ,'Temperature_ct'                                                                       ,'';
   'Temperature_MW'     ,'Microwave-only atmospheric temperature profile in 24 standard pressure levels from 1000 to 1.0 mbar' ,'K';
   'Temperature_MW_sdev','Temperature_MW_sdev'                                                                                 ,'K';
   'Temperature_MW_ct'  ,'Temperature_MW_ct'                                                                                   ,'';
   'GPHeight'           ,'Geopotential height at 24 standard pressure levels from 1000 to 1.0 mbar' ,'m';
   'GPHeight_sdev'      ,'GPHeight_sdev'                                                            ,'m';
   'GPHeight_ct'        ,'GPHeight_ct'                                                              ,'';
   'GPHeight_MW'        ,'MW-Only geopotential height in meters at 24 standard pressure levels from 1000 to 1.0 mbar'   ,'m';
   'GPHeight_MW_sdev'   ,'GPHeight_MW_sdev'                                                                             ,'m';
   'GPHeight_MW_ct'     ,'GPHeight_MW_ct'                                                                               ,'';
   'RelHumid'           ,'Relative humidity profile in 12 standard pressure levels from 1000 to 100 mbar'   ,'percent';
   'RelHumid_sdev'      ,'RelHumid_sdev'                                                                    ,'percent';
   'RelHumid_ct'        ,'RelHumid_ct'                                                                      ,'';
   'H2OVapMMR'          ,'Water vapor mass mixing ratio at 12 standard pressure levels from 1000 to 100 mbar'   ,'gm/kg dry air';
   'H2OVapMMR_sdev'     ,'H2OVapMMR_sdev'                                                                       ,'gm/kg dry air';
   'H2OVapMMR_ct'       ,'H2OVapMMR_ct'                                                                         ,'';
   'EmisIR'             ,'IR surface emissivity on a frequency grid 832, 961, 1203, 2616 cm-1'       ,'';
   'EmisIR_sdev'        ,'EmisIR_sdev'                                                               ,'';
   'EmisIR_ct'          ,'EmisIR_ct'                                                                 ,'';
   'EmisMW_MW'          ,'Microwave spectral emissivity on a frequency grid 23.8, 50.3 and 89.0 GHz' ,'';
   'EmisMW_MW_sdev'     ,'EmisMW_MW_sdev'                                                            ,'';
   'EmisMW_MW_ct'       ,'EmisMW_MW_ct'                                                              ,'';
   'LandMask'           ,'Land Sea Mask'                                                             ,'';
};

DSI.Fields         =FieldsTable(:,1);
DSI.Fields_NameDS  =FieldsTable(:,1);
DSI.Fields_NameMenu=FieldsTable(:,1);
DSI.Fields_NameLong=FieldsTable(:,2);
%DSI.Fields_Dim     =FieldsTable(:,5);
%DSI.Fields_Coord   =FieldsTable(:,6);
%DSI.Fields_DimDS   =FieldsTable(:,7);
%DSI.Fields_CoordDS =FieldsTable(:,8);
%DSI.Fields_CBNum   =FieldsTable(:,9);
DSI.Fields_Units   =FieldsTable(:,3);

for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
eval(['DSI.Attributes.' FIELD '.long_name=FieldsTable{k,2};']);
eval(['DSI.Attributes.' FIELD '.units    =FieldsTable{k,3};']);
end

CoordinatesTable={
'latitude','latitude','degrees_north';
'longitude','longitude','degrees_east';
'time','time','days since 0000-01-01 00:00:00';
'TempPresLvls','Pressure Levels for Temperature','mbar';
'H2OPresLvls','Pressure Levels for Water Vapor','mbar';
'IREmisFreqs','InfraRed Emissivity Frequencies','cm-1';
'MWEmisFreqs','MicroWave Emissivity Frequencies','GHz';
};

DSI.Coordinates=CoordinatesTable(:,1);
for k=1:length(DSI.Coordinates)
    COORD=DSI.Coordinates{k};
eval(['DSI.Attributes.' COORD '.long_name=CoordinatesTable{k,2};']);
eval(['DSI.Attributes.' COORD '.units    =CoordinatesTable{k,3};']);
end

DSI.Latitude=(89.5:-1:-89.5)';
DSI.Longitude=(-179.5:1:179.5)';
%TempPresLvls, in mbar 
DSI.TempPresLvls=[1000.0;925.0;850.0;700.0;600.0;500.0;400.0;300.0;250.0;200.0;150.0;100.0;70.0;50.0;30.0;20.0;15.0;10.0;7.0;5.0;3.0;2.0;1.5;1.0];
%H2OPresLvls, in mbar
DSI.H2OPresLvls=DSI.TempPresLvls(1:12);
%IREmisFreqs Frequencies for emissivities 
% reported in the AIRS Level-3 product in cm-1
DSI.IREmisFreqs=[832; 961; 1203; 2616];
% MWEmisFreqs Frequencies for microwave emissivity products 
% reported in AIRS Level-3 in GHz.
DSI.MWEmisFreqs=[23.8; 50.3; 89.0];

DSI.Range.Latitude1=min(DSI.Latitude);
DSI.Range.Latitude2=max(DSI.Latitude);
DSI.Range.Longitude1=min(DSI.Longitude);
DSI.Range.Longitude2=max(DSI.Longitude);
R.TAVG=DSI.TAVG;R.DATEINCR=Inf;U=cat_airs(R);
DSI.Range.Time1=U.TimeMin;
DSI.Range.Time2=U.TimeMax;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DSI=DSI_AIRS_weekly
% create DataSetInventory of AIRS
DSI=DSI_AIRS_daily;
% weekly
%http://acdisc.sci.gsfc.nasa.gov/opendap/Aqua_AIRS_Level3/AIRX3ST8.005/2007/AIRS.2007.01.02.L3.RetStd008.v5.0.14.0.G07219003002.hdf?TotCldLiqH2O_D[88:1:89][135:1:136],TotCldLiqH2O_D_sdev[88:1:89][135:1:136],TotCldLiqH2O_D_ct[88:1:89][135:1:136],TotH2OVap_D[88:1:89][135:1:136],TotH2OVap_D_sdev[88:1:89][135:1:136],TotH2OVap_D_ct[88:1:89][135:1:136],TotH2OVap_MW_D[88:1:89][135:1:136],TotH2OVap_MW_D_sdev[88:1:89][135:1:136],TotH2OVap_MW_D_ct[88:1:89][135:1:136],SurfAirTemp_D[88:1:89][135:1:136],SurfAirTemp_D_sdev[88:1:89][135:1:136],SurfAirTemp_D_ct[88:1:89][135:1:136],SurfSkinTemp_D[88:1:89][135:1:136],SurfSkinTemp_D_sdev[88:1:89][135:1:136],SurfSkinTemp_D_ct[88:1:89][135:1:136],SurfPres_D[88:1:89][135:1:136],SurfPres_D_sdev[88:1:89][135:1:136],SurfPres_D_ct[88:1:89][135:1:136],TotO3_D[88:1:89][135:1:136],TotO3_D_sdev[88:1:89][135:1:136],TotO3_D_ct[88:1:89][135:1:136],CloudFrc_D[88:1:89][135:1:136],CloudFrc_D_sdev[88:1:89][135:1:136],CloudFrc_D_ct[88:1:89][135:1:136],CloudTopPres_D[88:1:89][135:1:136],CloudTopPres_D_sdev[88:1:89][135:1:136],CloudTopPres_D_ct[88:1:89][135:1:136],CloudTopTemp_D[88:1:89][135:1:136],CloudTopTemp_D_sdev[88:1:89][135:1:136],CloudTopTemp_D_ct[88:1:89][135:1:136],OLR_D[88:1:89][135:1:136],OLR_D_sdev[88:1:89][135:1:136],OLR_D_ct[88:1:89][135:1:136],ClrOLR_D[88:1:89][135:1:136],ClrOLR_D_sdev[88:1:89][135:1:136],ClrOLR_D_ct[88:1:89][135:1:136],Temperature_D[0:1:2][88:1:89][135:1:136],Temperature_D_sdev[0:1:2][88:1:89][135:1:136],Temperature_D_ct[0:1:2][88:1:89][135:1:136],Temperature_MW_D[0:1:2][88:1:89][135:1:136],Temperature_MW_D_sdev[0:1:2][88:1:89][135:1:136],Temperature_MW_D_ct[0:1:2][88:1:89][135:1:136],GPHeight_D[0:1:2][88:1:89][135:1:136],GPHeight_D_sdev[0:1:2][88:1:89][135:1:136],GPHeight_D_ct[0:1:2][88:1:89][135:1:136],GPHeight_MW_D[0:1:2][88:1:89][135:1:136],GPHeight_MW_D_sdev[0:1:2][88:1:89][135:1:136],GPHeight_MW_D_ct[0:1:2][88:1:89][135:1:136],RelHumid_D[0:1:2][88:1:89][135:1:136],RelHumid_D_sdev[0:1:2][88:1:89][135:1:136],RelHumid_D_ct[0:1:2][88:1:89][135:1:136],H2OVapMMR_D[0:1:2][88:1:89][135:1:136],H2OVapMMR_D_sdev[0:1:2][88:1:89][135:1:136],H2OVapMMR_D_ct[0:1:2][88:1:89][135:1:136],EmisIR_D[0:1:2][88:1:89][135:1:136],EmisIR_D_sdev[0:1:2][88:1:89][135:1:136],EmisIR_D_ct[0:1:2][88:1:89][135:1:136],EmisMW_MW_D[0:1:2][88:1:89][135:1:136],EmisMW_MW_D_sdev[0:1:2][88:1:89][135:1:136],EmisMW_MW_D_ct[0:1:2][88:1:89][135:1:136],LandSeaMask[88:1:89][135:1:136]
%http://acdisc.sci.gsfc.nasa.gov/opendap/Aqua_AIRS_Level3/AIRX3ST8.005/2007/AIRS.2007.01.02.L3.RetStd008.v5.0.14.0.G07219003002.hdf?TotCldLiqH2O_A[88:1:89][135:1:136],TotCldLiqH2O_A_sdev[88:1:89][135:1:136],TotCldLiqH2O_A_ct[88:1:89][135:1:136],TotH2OVap_A[88:1:89][135:1:136],TotH2OVap_A_sdev[88:1:89][135:1:136],TotH2OVap_A_ct[88:1:89][135:1:136],TotH2OVap_MW_A[88:1:89][135:1:136],TotH2OVap_MW_A_sdev[88:1:89][135:1:136],TotH2OVap_MW_A_ct[88:1:89][135:1:136],SurfAirTemp_A[88:1:89][135:1:136],SurfAirTemp_A_sdev[88:1:89][135:1:136],SurfAirTemp_A_ct[88:1:89][135:1:136],SurfSkinTemp_A[88:1:89][135:1:136],SurfSkinTemp_A_sdev[88:1:89][135:1:136],SurfSkinTemp_A_ct[88:1:89][135:1:136],SurfPres_A[88:1:89][135:1:136],SurfPres_A_sdev[88:1:89][135:1:136],SurfPres_A_ct[88:1:89][135:1:136],TotO3_A[88:1:89][135:1:136],TotO3_A_sdev[88:1:89][135:1:136],TotO3_A_ct[88:1:89][135:1:136],CloudFrc_A[88:1:89][135:1:136],CloudFrc_A_sdev[88:1:89][135:1:136],CloudFrc_A_ct[88:1:89][135:1:136],CloudTopPres_A[88:1:89][135:1:136],CloudTopPres_A_sdev[88:1:89][135:1:136],CloudTopPres_A_ct[88:1:89][135:1:136],CloudTopTemp_A[88:1:89][135:1:136],CloudTopTemp_A_sdev[88:1:89][135:1:136],CloudTopTemp_A_ct[88:1:89][135:1:136],OLR_A[88:1:89][135:1:136],OLR_A_sdev[88:1:89][135:1:136],OLR_A_ct[88:1:89][135:1:136],ClrOLR_A[88:1:89][135:1:136],ClrOLR_A_sdev[88:1:89][135:1:136],ClrOLR_A_ct[88:1:89][135:1:136],Temperature_A[0:1:2][88:1:89][135:1:136],Temperature_A_sdev[0:1:2][88:1:89][135:1:136],Temperature_A_ct[0:1:2][88:1:89][135:1:136],Temperature_MW_A[0:1:2][88:1:89][135:1:136],Temperature_MW_A_sdev[0:1:2][88:1:89][135:1:136],Temperature_MW_A_ct[0:1:2][88:1:89][135:1:136],GPHeight_A[0:1:2][88:1:89][135:1:136],GPHeight_A_sdev[0:1:2][88:1:89][135:1:136],GPHeight_A_ct[0:1:2][88:1:89][135:1:136],GPHeight_MW_A[0:1:2][88:1:89][135:1:136],GPHeight_MW_A_sdev[0:1:2][88:1:89][135:1:136],GPHeight_MW_A_ct[0:1:2][88:1:89][135:1:136],RelHumid_A[0:1:2][88:1:89][135:1:136],RelHumid_A_sdev[0:1:2][88:1:89][135:1:136],RelHumid_A_ct[0:1:2][88:1:89][135:1:136],H2OVapMMR_A[0:1:2][88:1:89][135:1:136],H2OVapMMR_A_sdev[0:1:2][88:1:89][135:1:136],H2OVapMMR_A_ct[0:1:2][88:1:89][135:1:136],EmisIR_A[0:1:2][88:1:89][135:1:136],EmisIR_A_sdev[0:1:2][88:1:89][135:1:136],EmisIR_A_ct[0:1:2][88:1:89][135:1:136],EmisMW_MW_A[0:1:2][88:1:89][135:1:136],EmisMW_MW_A_sdev[0:1:2][88:1:89][135:1:136],EmisMW_MW_A_ct[0:1:2][88:1:89][135:1:136],LandSeaMask[88:1:89][135:1:136]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DSI.DataSetBranch='AIRS_weekly';
DSI.URLFILE='/AIRX3ST8.005/';
DSI.TAVG='weekly';
R.TAVG=DSI.TAVG;R.DATEINCR=Inf;U=cat_airs(R);
DSI.Range.Time1=U.TimeMin;
DSI.Range.Time2=U.TimeMax;


function DSI=DSI_AIRS_monthly
% create DataSetInventory of AIRS
DSI=DSI_AIRS_daily;
% monthly
%http://acdisc.sci.gsfc.nasa.gov/opendap/Aqua_AIRS_Level3/AIRX3STM.005/2007/AIRS.2007.01.01.L3.RetStd031.v5.0.14.0.G07348161630.hdf?TotCldLiqH2O_D[88:1:89][135:1:136],TotCldLiqH2O_D_sdev[88:1:89][135:1:136],TotCldLiqH2O_D_ct[88:1:89][135:1:136],TotH2OVap_D[88:1:89][135:1:136],TotH2OVap_D_sdev[88:1:89][135:1:136],TotH2OVap_D_ct[88:1:89][135:1:136],TotH2OVap_MW_D[88:1:89][135:1:136],TotH2OVap_MW_D_sdev[88:1:89][135:1:136],TotH2OVap_MW_D_ct[88:1:89][135:1:136],SurfAirTemp_D[88:1:89][135:1:136],SurfAirTemp_D_sdev[88:1:89][135:1:136],SurfAirTemp_D_ct[88:1:89][135:1:136],SurfSkinTemp_D[88:1:89][135:1:136],SurfSkinTemp_D_sdev[88:1:89][135:1:136],SurfSkinTemp_D_ct[88:1:89][135:1:136],SurfPres_D[88:1:89][135:1:136],SurfPres_D_sdev[88:1:89][135:1:136],SurfPres_D_ct[88:1:89][135:1:136],TotO3_D[88:1:89][135:1:136],TotO3_D_sdev[88:1:89][135:1:136],TotO3_D_ct[88:1:89][135:1:136],CloudFrc_D[88:1:89][135:1:136],CloudFrc_D_sdev[88:1:89][135:1:136],CloudFrc_D_ct[88:1:89][135:1:136],CloudTopPres_D[88:1:89][135:1:136],CloudTopPres_D_sdev[88:1:89][135:1:136],CloudTopPres_D_ct[88:1:89][135:1:136],CloudTopTemp_D[88:1:89][135:1:136],CloudTopTemp_D_sdev[88:1:89][135:1:136],CloudTopTemp_D_ct[88:1:89][135:1:136],OLR_D[88:1:89][135:1:136],OLR_D_sdev[88:1:89][135:1:136],OLR_D_ct[88:1:89][135:1:136],ClrOLR_D[88:1:89][135:1:136],ClrOLR_D_sdev[88:1:89][135:1:136],ClrOLR_D_ct[88:1:89][135:1:136],Temperature_D[0:1:2][88:1:89][135:1:136],Temperature_D_sdev[0:1:2][88:1:89][135:1:136],Temperature_D_ct[0:1:2][88:1:89][135:1:136],Temperature_MW_D[0:1:2][88:1:89][135:1:136],Temperature_MW_D_sdev[0:1:2][88:1:89][135:1:136],Temperature_MW_D_ct[0:1:2][88:1:89][135:1:136],GPHeight_D[0:1:2][88:1:89][135:1:136],GPHeight_D_sdev[0:1:2][88:1:89][135:1:136],GPHeight_D_ct[0:1:2][88:1:89][135:1:136],GPHeight_MW_D[0:1:2][88:1:89][135:1:136],GPHeight_MW_D_sdev[0:1:2][88:1:89][135:1:136],GPHeight_MW_D_ct[0:1:2][88:1:89][135:1:136],RelHumid_D[0:1:2][88:1:89][135:1:136],RelHumid_D_sdev[0:1:2][88:1:89][135:1:136],RelHumid_D_ct[0:1:2][88:1:89][135:1:136],H2OVapMMR_D[0:1:2][88:1:89][135:1:136],H2OVapMMR_D_sdev[0:1:2][88:1:89][135:1:136],H2OVapMMR_D_ct[0:1:2][88:1:89][135:1:136],EmisIR_D[0:1:2][88:1:89][135:1:136],EmisIR_D_sdev[0:1:2][88:1:89][135:1:136],EmisIR_D_ct[0:1:2][88:1:89][135:1:136],EmisMW_MW_D[0:1:2][88:1:89][135:1:136],EmisMW_MW_D_sdev[0:1:2][88:1:89][135:1:136],EmisMW_MW_D_ct[0:1:2][88:1:89][135:1:136],LandSeaMask[88:1:89][135:1:136]
%http://acdisc.sci.gsfc.nasa.gov/opendap/Aqua_AIRS_Level3/AIRX3STM.005/2007/AIRS.2007.01.01.L3.RetStd031.v5.0.14.0.G07348161630.hdf?TotCldLiqH2O_A[88:1:89][135:1:136],TotCldLiqH2O_A_sdev[88:1:89][135:1:136],TotCldLiqH2O_A_ct[88:1:89][135:1:136],TotH2OVap_A[88:1:89][135:1:136],TotH2OVap_A_sdev[88:1:89][135:1:136],TotH2OVap_A_ct[88:1:89][135:1:136],TotH2OVap_MW_A[88:1:89][135:1:136],TotH2OVap_MW_A_sdev[88:1:89][135:1:136],TotH2OVap_MW_A_ct[88:1:89][135:1:136],SurfAirTemp_A[88:1:89][135:1:136],SurfAirTemp_A_sdev[88:1:89][135:1:136],SurfAirTemp_A_ct[88:1:89][135:1:136],SurfSkinTemp_A[88:1:89][135:1:136],SurfSkinTemp_A_sdev[88:1:89][135:1:136],SurfSkinTemp_A_ct[88:1:89][135:1:136],SurfPres_A[88:1:89][135:1:136],SurfPres_A_sdev[88:1:89][135:1:136],SurfPres_A_ct[88:1:89][135:1:136],TotO3_A[88:1:89][135:1:136],TotO3_A_sdev[88:1:89][135:1:136],TotO3_A_ct[88:1:89][135:1:136],CloudFrc_A[88:1:89][135:1:136],CloudFrc_A_sdev[88:1:89][135:1:136],CloudFrc_A_ct[88:1:89][135:1:136],CloudTopPres_A[88:1:89][135:1:136],CloudTopPres_A_sdev[88:1:89][135:1:136],CloudTopPres_A_ct[88:1:89][135:1:136],CloudTopTemp_A[88:1:89][135:1:136],CloudTopTemp_A_sdev[88:1:89][135:1:136],CloudTopTemp_A_ct[88:1:89][135:1:136],OLR_A[88:1:89][135:1:136],OLR_A_sdev[88:1:89][135:1:136],OLR_A_ct[88:1:89][135:1:136],ClrOLR_A[88:1:89][135:1:136],ClrOLR_A_sdev[88:1:89][135:1:136],ClrOLR_A_ct[88:1:89][135:1:136],Temperature_A[0:1:2][88:1:89][135:1:136],Temperature_A_sdev[0:1:2][88:1:89][135:1:136],Temperature_A_ct[0:1:2][88:1:89][135:1:136],Temperature_MW_A[0:1:2][88:1:89][135:1:136],Temperature_MW_A_sdev[0:1:2][88:1:89][135:1:136],Temperature_MW_A_ct[0:1:2][88:1:89][135:1:136],GPHeight_A[0:1:2][88:1:89][135:1:136],GPHeight_A_sdev[0:1:2][88:1:89][135:1:136],GPHeight_A_ct[0:1:2][88:1:89][135:1:136],GPHeight_MW_A[0:1:2][88:1:89][135:1:136],GPHeight_MW_A_sdev[0:1:2][88:1:89][135:1:136],GPHeight_MW_A_ct[0:1:2][88:1:89][135:1:136],RelHumid_A[0:1:2][88:1:89][135:1:136],RelHumid_A_sdev[0:1:2][88:1:89][135:1:136],RelHumid_A_ct[0:1:2][88:1:89][135:1:136],H2OVapMMR_A[0:1:2][88:1:89][135:1:136],H2OVapMMR_A_sdev[0:1:2][88:1:89][135:1:136],H2OVapMMR_A_ct[0:1:2][88:1:89][135:1:136],EmisIR_A[0:1:2][88:1:89][135:1:136],EmisIR_A_sdev[0:1:2][88:1:89][135:1:136],EmisIR_A_ct[0:1:2][88:1:89][135:1:136],EmisMW_MW_A[0:1:2][88:1:89][135:1:136],EmisMW_MW_A_sdev[0:1:2][88:1:89][135:1:136],EmisMW_MW_A_ct[0:1:2][88:1:89][135:1:136],LandSeaMask[88:1:89][135:1:136]
%http://acdisc.sci.gsfc.nasa.gov/opendap/Aqua_AIRS_Level3/AIRX3STM.005/2007/AIRS.2007.02.01.L3.RetStd028.v5.0.14.0.G07351000036.hdf?TotCldLiqH2O_D[88:1:89][135:1:136],TotCldLiqH2O_D_sdev[88:1:89][135:1:136],TotCldLiqH2O_D_ct[88:1:89][135:1:136],TotH2OVap_D[88:1:89][135:1:136],TotH2OVap_D_sdev[88:1:89][135:1:136],TotH2OVap_D_ct[88:1:89][135:1:136],TotH2OVap_MW_D[88:1:89][135:1:136],TotH2OVap_MW_D_sdev[88:1:89][135:1:136],TotH2OVap_MW_D_ct[88:1:89][135:1:136],SurfAirTemp_D[88:1:89][135:1:136],SurfAirTemp_D_sdev[88:1:89][135:1:136],SurfAirTemp_D_ct[88:1:89][135:1:136],SurfSkinTemp_D[88:1:89][135:1:136],SurfSkinTemp_D_sdev[88:1:89][135:1:136],SurfSkinTemp_D_ct[88:1:89][135:1:136],SurfPres_D[88:1:89][135:1:136],SurfPres_D_sdev[88:1:89][135:1:136],SurfPres_D_ct[88:1:89][135:1:136],TotO3_D[88:1:89][135:1:136],TotO3_D_sdev[88:1:89][135:1:136],TotO3_D_ct[88:1:89][135:1:136],CloudFrc_D[88:1:89][135:1:136],CloudFrc_D_sdev[88:1:89][135:1:136],CloudFrc_D_ct[88:1:89][135:1:136],CloudTopPres_D[88:1:89][135:1:136],CloudTopPres_D_sdev[88:1:89][135:1:136],CloudTopPres_D_ct[88:1:89][135:1:136],CloudTopTemp_D[88:1:89][135:1:136],CloudTopTemp_D_sdev[88:1:89][135:1:136],CloudTopTemp_D_ct[88:1:89][135:1:136],OLR_D[88:1:89][135:1:136],OLR_D_sdev[88:1:89][135:1:136],OLR_D_ct[88:1:89][135:1:136],ClrOLR_D[88:1:89][135:1:136],ClrOLR_D_sdev[88:1:89][135:1:136],ClrOLR_D_ct[88:1:89][135:1:136],Temperature_D[0:1:2][88:1:89][135:1:136],Temperature_D_sdev[0:1:2][88:1:89][135:1:136],Temperature_D_ct[0:1:2][88:1:89][135:1:136],Temperature_MW_D[0:1:2][88:1:89][135:1:136],Temperature_MW_D_sdev[0:1:2][88:1:89][135:1:136],Temperature_MW_D_ct[0:1:2][88:1:89][135:1:136],GPHeight_D[0:1:2][88:1:89][135:1:136],GPHeight_D_sdev[0:1:2][88:1:89][135:1:136],GPHeight_D_ct[0:1:2][88:1:89][135:1:136],GPHeight_MW_D[0:1:2][88:1:89][135:1:136],GPHeight_MW_D_sdev[0:1:2][88:1:89][135:1:136],GPHeight_MW_D_ct[0:1:2][88:1:89][135:1:136],RelHumid_D[0:1:2][88:1:89][135:1:136],RelHumid_D_sdev[0:1:2][88:1:89][135:1:136],RelHumid_D_ct[0:1:2][88:1:89][135:1:136],H2OVapMMR_D[0:1:2][88:1:89][135:1:136],H2OVapMMR_D_sdev[0:1:2][88:1:89][135:1:136],H2OVapMMR_D_ct[0:1:2][88:1:89][135:1:136],EmisIR_D[0:1:2][88:1:89][135:1:136],EmisIR_D_sdev[0:1:2][88:1:89][135:1:136],EmisIR_D_ct[0:1:2][88:1:89][135:1:136],EmisMW_MW_D[0:1:2][88:1:89][135:1:136],EmisMW_MW_D_sdev[0:1:2][88:1:89][135:1:136],EmisMW_MW_D_ct[0:1:2][88:1:89][135:1:136],LandSeaMask[88:1:89][135:1:136]
%http://acdisc.sci.gsfc.nasa.gov/opendap/Aqua_AIRS_Level3/AIRX3STM.005/2007/AIRS.2007.02.01.L3.RetStd028.v5.0.14.0.G07351000036.hdf?TotCldLiqH2O_A[88:1:89][135:1:136],TotCldLiqH2O_A_sdev[88:1:89][135:1:136],TotCldLiqH2O_A_ct[88:1:89][135:1:136],TotH2OVap_A[88:1:89][135:1:136],TotH2OVap_A_sdev[88:1:89][135:1:136],TotH2OVap_A_ct[88:1:89][135:1:136],TotH2OVap_MW_A[88:1:89][135:1:136],TotH2OVap_MW_A_sdev[88:1:89][135:1:136],TotH2OVap_MW_A_ct[88:1:89][135:1:136],SurfAirTemp_A[88:1:89][135:1:136],SurfAirTemp_A_sdev[88:1:89][135:1:136],SurfAirTemp_A_ct[88:1:89][135:1:136],SurfSkinTemp_A[88:1:89][135:1:136],SurfSkinTemp_A_sdev[88:1:89][135:1:136],SurfSkinTemp_A_ct[88:1:89][135:1:136],SurfPres_A[88:1:89][135:1:136],SurfPres_A_sdev[88:1:89][135:1:136],SurfPres_A_ct[88:1:89][135:1:136],TotO3_A[88:1:89][135:1:136],TotO3_A_sdev[88:1:89][135:1:136],TotO3_A_ct[88:1:89][135:1:136],CloudFrc_A[88:1:89][135:1:136],CloudFrc_A_sdev[88:1:89][135:1:136],CloudFrc_A_ct[88:1:89][135:1:136],CloudTopPres_A[88:1:89][135:1:136],CloudTopPres_A_sdev[88:1:89][135:1:136],CloudTopPres_A_ct[88:1:89][135:1:136],CloudTopTemp_A[88:1:89][135:1:136],CloudTopTemp_A_sdev[88:1:89][135:1:136],CloudTopTemp_A_ct[88:1:89][135:1:136],OLR_A[88:1:89][135:1:136],OLR_A_sdev[88:1:89][135:1:136],OLR_A_ct[88:1:89][135:1:136],ClrOLR_A[88:1:89][135:1:136],ClrOLR_A_sdev[88:1:89][135:1:136],ClrOLR_A_ct[88:1:89][135:1:136],Temperature_A[0:1:2][88:1:89][135:1:136],Temperature_A_sdev[0:1:2][88:1:89][135:1:136],Temperature_A_ct[0:1:2][88:1:89][135:1:136],Temperature_MW_A[0:1:2][88:1:89][135:1:136],Temperature_MW_A_sdev[0:1:2][88:1:89][135:1:136],Temperature_MW_A_ct[0:1:2][88:1:89][135:1:136],GPHeight_A[0:1:2][88:1:89][135:1:136],GPHeight_A_sdev[0:1:2][88:1:89][135:1:136],GPHeight_A_ct[0:1:2][88:1:89][135:1:136],GPHeight_MW_A[0:1:2][88:1:89][135:1:136],GPHeight_MW_A_sdev[0:1:2][88:1:89][135:1:136],GPHeight_MW_A_ct[0:1:2][88:1:89][135:1:136],RelHumid_A[0:1:2][88:1:89][135:1:136],RelHumid_A_sdev[0:1:2][88:1:89][135:1:136],RelHumid_A_ct[0:1:2][88:1:89][135:1:136],H2OVapMMR_A[0:1:2][88:1:89][135:1:136],H2OVapMMR_A_sdev[0:1:2][88:1:89][135:1:136],H2OVapMMR_A_ct[0:1:2][88:1:89][135:1:136],EmisIR_A[0:1:2][88:1:89][135:1:136],EmisIR_A_sdev[0:1:2][88:1:89][135:1:136],EmisIR_A_ct[0:1:2][88:1:89][135:1:136],EmisMW_MW_A[0:1:2][88:1:89][135:1:136],EmisMW_MW_A_sdev[0:1:2][88:1:89][135:1:136],EmisMW_MW_A_ct[0:1:2][88:1:89][135:1:136],LandSeaMask[88:1:89][135:1:136]
DSI.DataSetBranch='AIRS_monthly';
DSI.URLFILE='/AIRX3STM.005/';
DSI.TAVG='monthly';
R.TAVG=DSI.TAVG;R.DATEINCR=Inf;U=cat_airs(R);
DSI.Range.Time1=U.TimeMin;
DSI.Range.Time2=U.TimeMax;

