function dsi_oaflux(varargin)
% updates the DataSetInventory structures
% without arguments for all the following DataSetBranches

if nargin < 1
%DSI_OAFlux
DataSetBranchList={
'OAFlux_daily';
'OAFlux_monthly';
};

dsi_filename ='dsi_oaflux.mat';
%if (~isdeployed)
pname = strrep(mfilename('fullpath'),mfilename,'');
dsi_filename = [pname,dsi_filename];
%end

if ~exist(dsi_filename,'file')
save(dsi_filename,'DataSetBranchList')
end
save(dsi_filename,'DataSetBranchList','-append')

for k=1:length(DataSetBranchList)    
eval(['DSI_' DataSetBranchList{k} ';'])
end

elseif nargin > 1
disp('Incorrect number of arguments.'); return
else
% with one argument for the specified DataSetBranch only
eval(varargin{1});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DSI=DSI_OAFlux_daily
% create DataSetInventory of OAFlux
%A=loaddap('-A','http://apdrc.soest.hawaii.edu:80/dods/public_data/WHOI_OAFlux/version3/daily/evapr_oaflux')
%
%A = 
%
%                evapr: [1x1 struct]
%                  err: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]
%    Global_Attributes: [1x1 struct]

DSI.DataSetName='OAFlux';
DSI.DataSetBranch='OAFlux_daily';

% a general request employed by get_DataSet program has the following pattern
% loaddap
DSI.LOADDAPURL='LOADDAPURL=[URLSITE URLPATH URLFILE ''?'' CVAR CTIME CLAT CLON]';
DSI.URLSITE='http://apdrc.soest.hawaii.edu:80/dods/public_data/WHOI_OAFlux/version3';
DSI.URLPATH='/daily/';
DSI.URLFILE='URLFILE=DSI.Fields_URLFile{j};';
DSI.CVAR='CVAR=DSI.Fields_NameDS{j};';

% define list of fields or variables available in the dataset
LaLo={'latitude','longitude'};
TLaLo={'time','lat','lon'};
FieldsTable={
%http://apdrc.soest.hawaii.edu:80/dods/public_data/WHOI_OAFlux/version3/daily/evapr_oaflux
%http://apdrc.soest.hawaii.edu:80/dods/public_data/WHOI_OAFlux/version3/daily/lh_oaflux
%http://apdrc.soest.hawaii.edu:80/dods/public_data/WHOI_OAFlux/version3/daily/lw_isccp

% opendap name  ,dataset name  , menu name                                         , long name,                                                                         URLFILE
    'nlwrs'     ,'nlwrs'       ,'net surface longwave radiation flux (nlwrs)'              ,'net surface longwave radiation flux (positive upward) [w/m/m]'   ,2,LaLo,3,TLaLo,1,'w/m/m','lw_isccp';
    'nswrs'     ,'nswrs'       ,'net surface shortwave radiation flux (nswrs)'             ,'net surface shortwave radiation flux (positive downward) [w/m/m]',2,LaLo,3,TLaLo,2,'w/m/m','sw_isccp'; 
    'qnet'      ,'qnet'        ,'net surface heat flux (qnet)'                             ,'net surface heat flux (positive down) [w/m/m]'                   ,2,LaLo,3,TLaLo,3,'w/m/m','qnet';
    'lhtfl'     ,'lhtfl'       ,'surface latent heat flux (lhtfl)'                         ,'daily mean surface latent heat flux, positive upward [w/m/m]'    ,2,LaLo,3,TLaLo,4,'w/m/m','lh_oaflux';
%  	'lhtfl_err' ,'err'         ,'error of surface latent heat flux (lhtfl_err)'            ,'estimated error of surface latent heat flux [w/m/m]'             ,2,LaLo,3,TLaLo,1,'w/m/m','lh_oaflux';
	'shtfl'     ,'shtfl'       ,'surface sensible heat flux (shtfl)'                       ,'daily mean surface sensible heat flux (positive upward) [w/m/m]' ,2,LaLo,3,TLaLo,5,'w/m/m','sh_oaflux';
%   'shtfl_err' ,'err'         ,'error of surface sensible heat flux (shtfl_err)'          ,'estimated error of surface sensible heat flux [w/m/m]'           ,2,LaLo,3,TLaLo,1,'w/m/m','sh_oaflux';
	'evapr'     ,'evapr'       ,'evaporation rate (evapr)'                                 ,'daily mean evaporation rate [cm/yr]'                             ,2,LaLo,3,TLaLo,6,'cm/yr','evapr_oaflux';
%   'evapr_err' ,'err'         ,'error of evaporation rate (evapr_err)'                    ,'daily mean estimated error of evaporation rate [cm/yr]'          ,2,LaLo,3,TLaLo,1,'cm/yr','evapr_oaflux';
    'hum2m'     ,'hum2m'       ,'specific humidity (hum2m)'                                ,'daily mean specific humidity at 2m [g/kg]'                       ,2,LaLo,3,TLaLo,7,'g/kg' ,'qa_oaflux';
% 	'hum2m_err' ,'err'         ,'error of specific humidity (hum2m_err)'                   ,'daily mean estimated error of specific humidity at 2m [g/kg]'    ,2,LaLo,3,TLaLo,1,'g/kg' ,'qa_oaflux';
    'tmp2m'     ,'tmp2m'       ,'air temperature at 2m (tmp2m)'                            ,'daily mean air temperature at 2m [deg C]'                         ,2,LaLo,3,TLaLo,8,'w/m/m','ta_oaflux';
% 	'tmp2m_err' ,'err'         ,'estimated error of air temperature at 2m (tmp2m_err)'     ,'daily mean estimated error of air temperature at 2m [deg C]'      ,2,LaLo,3,TLaLo,1,'deg C' ,'ta_oaflux';
    'tmpsf'     ,'tmpsf'       ,'sea surface temperature (tmpsf)'                          ,'daily mean sea surface temperature [deg C]'                       ,2,LaLo,3,TLaLo,9,'deg C' ,'ts_oaflux';
%   'tmpsf_err' ,'err'         ,'estimated error of sea surface temperature (tmpsf_err)'   ,'daily mean estimated error of sea surface temperature[deg C]'     ,2,LaLo,3,TLaLo,1,'deg C' ,'ts_oaflux';
 	'wnd10'     ,'wnd10'       ,'neutral wind speed at 10m (wnd10)'                        ,'daily mean neutral wind speed at 10m [m/s]'                      ,2,LaLo,3,TLaLo,10,'m/s' ,'ws_oaflux';
%   'wnd10_err' ,'err'         ,'estimated error of neutral wind speed at 10m (wnd10_err)' ,'daily mean estimated error of neutral wind speed at 10m [m/s]'   ,2,LaLo,3,TLaLo,1,'m/s'  ,'ws_oaflux';
};

DSI.Fields         =FieldsTable(:,1);
DSI.Fields_NameDS  =FieldsTable(:,2);
DSI.Fields_NameMenu=FieldsTable(:,3);
DSI.Fields_NameLong=FieldsTable(:,4);
DSI.Fields_Dim     =FieldsTable(:,5);
DSI.Fields_Coord   =FieldsTable(:,6);
DSI.Fields_DimDS   =FieldsTable(:,7);
DSI.Fields_CoordDS =FieldsTable(:,8);
DSI.Fields_CBNum   =FieldsTable(:,9);
DSI.Fields_Units   =FieldsTable(:,10);
DSI.Fields_URLFile =FieldsTable(:,11);
DSI.Fields_ToolTip =FieldsTable(:,4);


% Rename returned variables and change dimensions order
%f=squeeze(nlwrs.nlwrs);
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
    FREN=['f=squeeze(' DSI.Fields_NameDS{k} '.' DSI.Fields_NameDS{k} ');'];
eval(['DSI.FormulaRen.' FIELD '=''' FREN ''';']);      
end

% missing_value: 1.0E20 % same for all fields
%DSI.FormulaNaN='f(find(f==1.0e20))=NaN;'; % does not work
DSI.FormulaNaN='f(find(f>0.9e20))=NaN;';

% use evapr to get time other vars on same grid
%Time:00Z01JAN1985 to 00Z31DEC2006 (8035 points, avg. res. 1.0 days)
URL=[DSI.URLSITE DSI.URLPATH 'evapr_oaflux?time']
loaddap(['+v','-e'],URL)
if dods_err
%    disp(dods_err_msg)
    msgbox(['Could not access ' DSI.DataSetName ' site. Sorry, try again later.']);
    return
end

% probably mistake in calculating datenum of 0001-01-01 use reference to
% the first element instead
t0=datenum('1985-01-01 00:00:00','yyyy-mm-dd HH:MM:SS');
t=time-time(1)+t0; % in datenum
DSI.Time=t;

DSI.Latitude=(-89.5:1:89.5)';
DSI.Longitude=(0.5:1:359.5)';

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

DSI.Range.Time1=min(DSI.Time);
DSI.Range.Time2=max(DSI.Time);
DSI.Range.Latitude1=min(DSI.Latitude);
DSI.Range.Latitude2=max(DSI.Latitude);
DSI.Range.Longitude1=min(DSI.Longitude);
DSI.Range.Longitude2=max(DSI.Longitude);
DSI.Range.TRes='daily';
DSI.Range.HRes='1 deg';

DSI.Range.Time1_nlwrs=datenum(1985,01,01);
DSI.Range.Time1_nswrs=datenum(1985,01,01);
DSI.Range.Time1_qnet =datenum(1985,01,01);

DSI.Range.Time2_nlwrs=datenum(2004,12,31);
DSI.Range.Time2_nswrs=datenum(2004,12,31);
DSI.Range.Time2_qnet =datenum(2004,12,31);

S1=['time in seconds since reference_time;'];
S2=['missing values replaced with NaN;'];
S3=['variables converted to physical units.'];
%DSI.Readme=[S1 '/' S2 '/' S3];
DSI.Readme=[S1 '/' S2];

eval(['DSI_' DSI.DataSetBranch '=DSI;'])
dsi_filename ='dsi_oaflux.mat';
%if (~isdeployed)
pname = strrep(mfilename('fullpath'),mfilename,'');
dsi_filename = [pname,dsi_filename];
%end
if ~exist(dsi_filename,'file')
save(dsi_filename,['DSI_' DSI.DataSetBranch])
end
save(dsi_filename,['DSI_' DSI.DataSetBranch],'-append')

function DSI=DSI_OAFlux_monthly
% create DataSetInventory of OAFlux
%A=loaddap('-A','http://apdrc.soest.hawaii.edu:80/dods/public_data/WHOI_OAFlux/version3/monthly/evapr_oaflux')
%
%A = 
%
%                evapr: [1x1 struct]
%                  err: [1x1 struct]
%                 time: [1x1 struct]
%                  lat: [1x1 struct]
%                  lon: [1x1 struct]
%    Global_Attributes: [1x1 struct]

DSI.DataSetName='OAFlux';
DSI.DataSetBranch='OAFlux_monthly';

% a general request employed by get_DataSet program has the following pattern
% loaddap
DSI.LOADDAPURL='LOADDAPURL=[URLSITE URLPATH URLFILE ''?'' CVAR CTIME CLAT CLON]';
DSI.URLSITE='http://apdrc.soest.hawaii.edu:80/dods/public_data/WHOI_OAFlux/version3';
DSI.URLPATH='/monthly/';
DSI.URLFILE='URLFILE=DSI.Fields_URLFile{j};';
DSI.CVAR='CVAR=DSI.Fields_NameDS{j};';

% define list of fields or variables available in the dataset
LaLo={'latitude','longitude'};
TLaLo={'time','lat','lon'};
FieldsTable={
% opendap name  ,dataset name  , menu name                                         , long name,                                                                         URLFILE
    'nlwrs'     ,'nlwrs'       ,'net surface longwave radiation flux (nlwrs)'              ,'net surface longwave radiation flux (positive upward) [w/m/m]'     ,2,LaLo,3,TLaLo,1,'w/m/m','lw_isccp';
    'nswrs'     ,'nswrs'       ,'net surface shortwave radiation flux (nswrs)'             ,'net surface shortwave radiation flux (positive downward) [w/m/m]'  ,2,LaLo,3,TLaLo,2,'w/m/m','sw_isccp'; 
    'qnet'      ,'qnet'        ,'net surface heat flux (qnet)'                             ,'net surface heat flux (positive down) [w/m/m]'                     ,2,LaLo,3,TLaLo,3,'w/m/m','qnet';
    'lhtfl'     ,'lhtfl'       ,'surface latent heat flux (lhtfl)'                         ,'monthly mean surface latent heat flux, positive upward [w/m/m]'    ,2,LaLo,3,TLaLo,4,'w/m/m','lh_oaflux';
%  	'lhtfl_err' ,'err'         ,'error of surface latent heat flux (lhtfl_err)'            ,'estimated error of surface latent heat flux [w/m/m]'               ,2,LaLo,3,TLaLo,1,'w/m/m','lh_oaflux';
	'shtfl'     ,'shtfl'       ,'surface sensible heat flux (shtfl)'                       ,'monthly mean surface sensible heat flux (positive upward) [w/m/m]' ,2,LaLo,3,TLaLo,5,'w/m/m','sh_oaflux';
%   'shtfl_err' ,'err'         ,'error of surface sensible heat flux (shtfl_err)'          ,'estimated error of surface sensible heat flux [w/m/m]'             ,2,LaLo,3,TLaLo,1,'w/m/m','sh_oaflux';
	'evapr'     ,'evapr'       ,'evaporation rate (evapr)'                                 ,'monthly mean evaporation rate [cm/yr]'                             ,2,LaLo,3,TLaLo,6,'cm/yr','evapr_oaflux';
%   'evapr_err' ,'err'         ,'error of evaporation rate (evapr_err)'                    ,'monthly mean estimated error of evaporation rate [cm/yr]'          ,2,LaLo,3,TLaLo,1,'cm/yr','evapr_oaflux';
    'hum2m'     ,'hum2m'       ,'specific humidity (hum2m)'                                ,'monthly mean specific humidity at 2m [g/kg]'                       ,2,LaLo,3,TLaLo,7,'g/kg' ,'qa_oaflux';
% 	'hum2m_err' ,'err'         ,'error of specific humidity (hum2m_err)'                   ,'monthly mean estimated error of specific humidity at 2m [g/kg]'    ,2,LaLo,3,TLaLo,1,'g/kg' ,'qa_oaflux';
    'tmp2m'     ,'tmp2m'       ,'air temperature at 2m (tmp2m)'                            ,'monthly mean air temperature at 2m [deg C]'                         ,2,LaLo,3,TLaLo,8,'w/m/m','ta_oaflux';
% 	'tmp2m_err' ,'err'         ,'estimated error of air temperature at 2m (tmp2m_err)'     ,'monthly mean estimated error of air temperature at 2m [deg C]'      ,2,LaLo,3,TLaLo,1,'deg C' ,'ta_oaflux';
    'tmpsf'     ,'tmpsf'       ,'sea surface temperature (tmpsf)'                          ,'monthly mean sea surface temperature [deg C]'                       ,2,LaLo,3,TLaLo,9,'deg C' ,'ts_oaflux';
%   'tmpsf_err' ,'err'         ,'estimated error of sea surface temperature (tmpsf_err)'   ,'monthly mean estimated error of sea surface temperature[deg C]'     ,2,LaLo,3,TLaLo,1,'deg C' ,'ts_oaflux';
 	'wnd10'     ,'wnd10'       ,'neutral wind speed at 10m (wnd10)'                        ,'monthly mean neutral wind speed at 10m [m/s]'                      ,2,LaLo,3,TLaLo,10,'m/s' ,'ws_oaflux';
%   'wnd10_err' ,'err'         ,'estimated error of neutral wind speed at 10m (wnd10_err)' ,'monthly mean estimated error of neutral wind speed at 10m [m/s]'   ,2,LaLo,3,TLaLo,1,'m/s'  ,'ws_oaflux';
};

DSI.Fields         =FieldsTable(:,1);
DSI.Fields_NameDS  =FieldsTable(:,2);
DSI.Fields_NameMenu=FieldsTable(:,3);
DSI.Fields_NameLong=FieldsTable(:,4);
DSI.Fields_Dim     =FieldsTable(:,5);
DSI.Fields_Coord   =FieldsTable(:,6);
DSI.Fields_DimDS   =FieldsTable(:,7);
DSI.Fields_CoordDS =FieldsTable(:,8);
DSI.Fields_CBNum   =FieldsTable(:,9);
DSI.Fields_Units   =FieldsTable(:,10);
DSI.Fields_URLFile =FieldsTable(:,11);
DSI.Fields_ToolTip =FieldsTable(:,4);


% Rename returned variables and change dimensions order
%f=squeeze(nlwrs.nlwrs);
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
    FREN=['f=squeeze(' DSI.Fields_NameDS{k} '.' DSI.Fields_NameDS{k} ');'];
eval(['DSI.FormulaRen.' FIELD '=''' FREN ''';']);      
end

% missing_value: 1.0E20 % same for all fields
%DSI.FormulaNaN='f(find(f==1.0e20))=NaN;'; % does not work
DSI.FormulaNaN='f(find(f>0.9e20))=NaN;';

% use evapr to get time other vars on same grid
%Time:00Z15JAN1958 to 00Z15DEC2006 (588 points, avg. res. 30.436 days) 
URL=[DSI.URLSITE DSI.URLPATH 'evapr_oaflux?time']
loaddap(['+v','-e'],URL)
if dods_err
%    disp(dods_err_msg)
    msgbox(['Could not access ' DSI.DataSetName ' site. Sorry, try again later.']);
    return
end

% probably mistake in calculating datenum of 0001-01-01 use reference to
% the first element instead
t0=datenum('1958-01-15 00:00:00','yyyy-mm-dd HH:MM:SS');
t=time-time(1)+t0; % in datenum
DSI.Time=t;

DSI.Latitude=(-89.5:1:89.5)';
DSI.Longitude=(0.5:1:359.5)';

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

DSI.Range.Time1=min(DSI.Time);
DSI.Range.Time2=max(DSI.Time);
DSI.Range.Latitude1=min(DSI.Latitude);
DSI.Range.Latitude2=max(DSI.Latitude);
DSI.Range.Longitude1=min(DSI.Longitude);
DSI.Range.Longitude2=max(DSI.Longitude);
DSI.Range.TRes='monthly';
DSI.Range.HRes='1 deg';

S1=['time in seconds since reference_time;'];
S2=['missing values replaced with NaN;'];
S3=['variables converted to physical units.'];
%DSI.Readme=[S1 '/' S2 '/' S3];
DSI.Readme=[S1 '/' S2];

eval(['DSI_' DSI.DataSetBranch '=DSI;'])
dsi_filename ='dsi_oaflux.mat';
%if (~isdeployed)
pname = strrep(mfilename('fullpath'),mfilename,'');
dsi_filename = [pname,dsi_filename];
%end
if ~exist(dsi_filename,'file')
save(dsi_filename,['DSI_' DSI.DataSetBranch])
end
save(dsi_filename,['DSI_' DSI.DataSetBranch],'-append')


