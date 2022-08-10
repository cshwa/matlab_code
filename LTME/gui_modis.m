function varargout = gui_modis(varargin)
% Boots-up MODIS GUI, interacts with user via interface.
% gui_modis M-file for gui_modis.fig
%
% INPUT(S):
%
%   [none]
%
% OUTPUT(S):
%
%   h -- (optional) GUI handle
%
% OPeNDAP Science Team
% Copyright 2007,2008,2009
% $Version 3.0.0$

% Meri Sheremet
% May 2009

warning('off','MATLAB:dispatcher:InexactMatch');
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_modis_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_modis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
warning('on','MATLAB:dispatcher:InexactMatch');

% --- Executes just before gui_modis is made visible.
function gui_modis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_modis (see VARARGIN)

DSI=load('dsi_modis.mat');
handles.DSI=DSI;

% Define handles "openfigure" field.
if ~isfield(handles,'openfigure')
    handles.openfigure = 0; %set equal to zero on the first time.
end

% Set interface name.
filename = 'get_modis';
txt = [''];
tmp1 = ['MODIS GUI'];
tmp2 = '';
txt = [txt,tmp1];
if (~isdeployed)
    tmp2 = strrep(op_version('mfile',filename),'Version ','');
    txt = [txt,' (Ver ',tmp2,')'];
end
set(hObject,'Name',txt);


% Set "SaveFiles" off. (THIS IS TEMPORARY ... 08 MAY 2007)
%set(handles.radiobuttonToFiles,'Enable','Off');

% Set Version Number.
filename = strrep(mfilename,'gui','get');
info = op_version(filename);
set(handles.textVer,'String',info);

% default GUI settings

%Place logos on the GUI    
% determine screen resolution and size of character in pixels
Units0=get(0,'Units'); % save initial settings
set(0,'Units','characters');
SSC=get(0,'ScreenSize'); % screen size in characters
set(0,'Units','pixels');
SSP=get(0,'ScreenSize'); % screen size in pixels
set(0,'Units',Units0); % restore to original units
ppcx=SSP(3)/SSC(3); % pixels per character along x
ppcy=SSP(4)/SSC(4); % pixels per character along y

[img,MAP]=imread('logo_opendap.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonOpendap,'CData',img)
pos = get(handles.pushbuttonOpendap,'position'); %units will be in characters
set(handles.pushbuttonOpendap,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_nasa_jpl.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonNASAjpl,'CData',img)
pos = get(handles.pushbuttonNASAjpl,'position'); %units will be in characters
set(handles.pushbuttonNASAjpl,'position',[pos(1),pos(2),CX,CY]);

set(handles.radiobuttonPlatformAqua,'Value',1)
set(handles.radiobuttonPassDay,'Value',1) 
set(handles.radiobuttonPassNight,'Value',1) 

set(handles.editTimeStep,'String',num2str(1));
set(handles.editLatStep,'String',num2str(1));
set(handles.editLonStep,'String',num2str(1));
%set(handles.editFileNamePrefix,'String','opendap_');

if get(handles.radiobuttonInWksp,'Value')==0 & get(handles.radiobuttonToFiles,'Value')==0
set(handles.radiobuttonInWksp,'Value',1)
end

if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end

set(handles.pushbuttonGetData,'String','Get Data')

% set time range
CheckRanges(hObject, eventdata, handles)

% Make 'Load Last Request' button inactive when file gui_xxx_request.mat is absent.
if ~exist('gui_modis_request.mat','file') 
    set(handles.pushbuttonLoadLastRequest,'Enable','Off')
end

% Choose default command line output for gui_modis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_modis wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_modis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function checkbox1_Callback(hObject, eventdata, handles)
function checkbox2_Callback(hObject, eventdata, handles)
function checkbox3_Callback(hObject, eventdata, handles)
function checkbox4_Callback(hObject, eventdata, handles)
function checkbox5_Callback(hObject, eventdata, handles)
function checkbox6_Callback(hObject, eventdata, handles)
function checkbox7_Callback(hObject, eventdata, handles)
function checkbox8_Callback(hObject, eventdata, handles)

function editTime1yyyy_CreateFcn(hObject, eventdata, handles)
function editTime1yyyy_Callback(hObject, eventdata, handles)
function editTime1mm_CreateFcn(hObject, eventdata, handles)
function editTime1mm_Callback(hObject, eventdata, handles)
function editTime1dd_CreateFcn(hObject, eventdata, handles)
function editTime1dd_Callback(hObject, eventdata, handles)
function editTime2yyyy_CreateFcn(hObject, eventdata, handles)
function editTime2yyyy_Callback(hObject, eventdata, handles)
function editTime2mm_CreateFcn(hObject, eventdata, handles)
function editTime2mm_Callback(hObject, eventdata, handles)
function editTime2dd_CreateFcn(hObject, eventdata, handles)
function editTime2dd_Callback(hObject, eventdata, handles)
function editLat1_CreateFcn(hObject, eventdata, handles)
function editLat1_Callback(hObject, eventdata, handles)
function editLat2_CreateFcn(hObject, eventdata, handles)
function editLat2_Callback(hObject, eventdata, handles)
function editLon1_CreateFcn(hObject, eventdata, handles)
function editLon1_Callback(hObject, eventdata, handles)
function editLon2_CreateFcn(hObject, eventdata, handles)
function editLon2_Callback(hObject, eventdata, handles)

function popupmenuTemporalAverages_CreateFcn(hObject, eventdata, handles)
function popupmenuTemporalAverages_Callback(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)
function popupmenuSpatialResolution_CreateFcn(hObject, eventdata, handles)
function popupmenuSpatialResolution_Callback(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)

function radiobuttonPlatformTerra_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonPlatformTerra,'Value') == 1 
set(handles.radiobuttonPlatformAqua,'Value',0);
else
set(handles.radiobuttonPlatformAqua,'Value',1);
end
CheckRanges(hObject, eventdata, handles)
function radiobuttonPlatformAqua_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonPlatformAqua,'Value') == 1 
set(handles.radiobuttonPlatformTerra,'Value',0);
else
set(handles.radiobuttonPlatformTerra,'Value',1);
end
CheckRanges(hObject, eventdata, handles)

function radiobuttonPassDay_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonPassDay,'Value') == 0 & get(handles.radiobuttonPassNight,'Value') == 0 
set(handles.radiobuttonPassDay,'Value',1);set(handles.radiobuttonPassNight,'Value',1);
end
CheckRanges(hObject, eventdata, handles)
function radiobuttonPassNight_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonPassDay,'Value') == 0 & get(handles.radiobuttonPassNight,'Value') == 0 
set(handles.radiobuttonPassDay,'Value',1);set(handles.radiobuttonPassNight,'Value',1);
end
CheckRanges(hObject, eventdata, handles)

function editLatStep_CreateFcn(hObject, eventdata, handles)
function editLatStep_Callback(hObject, eventdata, handles)
function editLonStep_CreateFcn(hObject, eventdata, handles)
function editLonStep_Callback(hObject, eventdata, handles)
function editDirectoryName_CreateFcn(hObject, eventdata, handles)
function editDirectoryName_Callback(hObject, eventdata, handles)
function editFileNamePrefix_CreateFcn(hObject, eventdata, handles)
function editFileNamePrefix_Callback(hObject, eventdata, handles)
function editTimeStep_Callback(hObject, eventdata, handles)
function editTimeStep_CreateFcn(hObject, eventdata, handles)

function radiobuttonInWksp_Callback(hObject, eventdata, handles)
function radiobuttonToFiles_Callback(hObject, eventdata, handles)
SaveFilesOnOff(hObject, eventdata, handles)

function popupmenuMode_Callback(hObject, eventdata, handles)
S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');
R.SaveMode=S{j};
% disable saving to workspace for netcdf mode
set(handles.radiobuttonInWksp,'Enable','On');
if strcmp(R.SaveMode,'netcdf')
    set(handles.radiobuttonInWksp,'Value',0,'Enable','Off');
    set(handles.radiobuttonToFiles,'Value',1);
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end
function popupmenuMode_CreateFcn(hObject, eventdata, handles)


function pushbuttonBrowse_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonToFiles,'Value') == 1
DIRNAME=uigetdir(pwd,'Select a Directory for Downloading Files');
set(handles.editDirectoryName,'String',DIRNAME)
end

% The main button that starts data acquisition
function pushbuttonGetData_Callback(hObject, eventdata, handles)
%% Turn off ability to re-issue "Get Data."
set(hObject,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off','String','Please wait...');
pause(.1)

%% GET INPUT FROM GUI 
request=InputGet(hObject, eventdata, handles);

% Save request to .mat file.
filename = 'gui_modis_request.mat';
%if (~isdeployed)
    pname = mfilename('fullpath');
    pname = pname(1:length(pname)-length(mfilename)); %mfilename does not have .m in it
    filename = [pname,filename];
%end
save(filename,'request','-mat')
set(handles.pushbuttonLoadLastRequest,'Enable','On')

% If it passes inspection, send down to user's workspace
% as variable "request". Then the user can get the data
assignin('base','request',request);

get_modis(request);

set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25]);

function pushbuttonOpendap_Callback(hObject, eventdata, handles)
web('http://www.opendap.org/');
%web 'http://www.opendap.org/' -browser
function pushbuttonNASAjpl_Callback(hObject, eventdata, handles)
web('http://podaac.jpl.nasa.gov:2031/DATASET_DOCS/modis_sst.html');
%web 'http://podaac.jpl.nasa.gov:2031/DATASET_DOCS/modis_sst.html' -browser

function request=InputGet(hObject, eventdata, handles)
S1=deblank(get(handles.editTime1yyyy,'String'));
S2=deblank(get(handles.editTime1mm,'String'));
S3=deblank(get(handles.editTime1dd,'String'));
S4=deblank(get(handles.editTime2yyyy,'String'));
S5=deblank(get(handles.editTime2mm,'String'));
S6=deblank(get(handles.editTime2dd,'String'));

DATE1='00000000';
DATE1(4-length(S1)+1:4)=S1;
DATE1(6-length(S2)+1:6)=S2;
DATE1(8-length(S3)+1:8)=S3;
DATE2='00000000';
DATE2(4-length(S4)+1:4)=S4;
DATE2(6-length(S5)+1:6)=S5;
DATE2(8-length(S6)+1:8)=S6;

LAT1=str2num(deblank(get(handles.editLat1,'String')));
LAT2=str2num(deblank(get(handles.editLat2,'String')));
LON1=str2num(deblank(get(handles.editLon1,'String')));
LON2=str2num(deblank(get(handles.editLon2,'String')));
LATINCR=str2num(deblank(get(handles.editLatStep,'String')));
LONINCR=str2num(deblank(get(handles.editLonStep,'String')));
DATEINCR=str2num(deblank(get(handles.editTimeStep,'String')));

DIRNAME =deblank(get(handles.editDirectoryName,'String'));
if isempty(DIRNAME) DIRNAME=''; end
FNPREFIX=deblank(get(handles.editFileNamePrefix,'String'));

SaveWorkspace='n'; 
if get(handles.radiobuttonInWksp,'Value') == 1 SaveWorkspace='y'; end
SaveFiles='n'; 
if get(handles.radiobuttonToFiles,'Value') == 1 SaveFiles='y'; end

Fields={};
if get(handles.checkbox1,'Value')==1 Fields{end+1}='sst'; end
if get(handles.checkbox2,'Value')==1 Fields{end+1}='sst4'; end
if get(handles.checkbox3,'Value')==1 Fields{end+1}='sst_qual'; end
if get(handles.checkbox4,'Value')==1 Fields{end+1}='sst4_qual'; end

ApplyLandMask='n';
%if get(handles.checkbox44 ,'Value')==1 ApplyLandMask='y'; end

ApplyQualMask='n';
%if get(handles.checkbox33 ,'Value')==1 
%S=get(handles.popupmenu11,'String');j=get(handles.popupmenu11,'Value');
%ApplyQualMask=S{j}(3); %'n','4','7'
%end

CalcAnomaly='n';
%if get(handles.checkbox7 ,'Value')==1 CalcAnomaly='y'; end

S=get(handles.popupmenuTemporalAverages,'String');j=get(handles.popupmenuTemporalAverages,'Value');
TAvg=deblank(S{j});
S=get(handles.popupmenuSpatialResolution,'String');j=get(handles.popupmenuSpatialResolution,'Value');
HRes=deblank(S{j});

if get(handles.radiobuttonPlatformTerra,'Value') == 1 Platform='terra'; end % terra
if get(handles.radiobuttonPlatformAqua,'Value') == 1 Platform='aqua'; end % aqua

if get(handles.radiobuttonPassDay,'Value') == 1 & get(handles.radiobuttonPassNight,'Value') == 0 
    Pass='day';  % day
elseif get(handles.radiobuttonPassDay,'Value') == 0 & get(handles.radiobuttonPassNight,'Value') == 1 
    Pass='night'; % night
else
    Pass='all';
end
    


S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');
SaveMode=deblank(S{j});

%% Get default command line output from handles structure
R.DataSetBranch=['MODIS_' Platform '_' HRes '_' TAvg '_' Pass];
R.DATE1=DATE1;
R.DATE2=DATE2;
R.DATEINCR=DATEINCR;
R.LAT1=LAT1;
R.LAT2=LAT2;
R.LATINCR=LATINCR;
R.LON1=LON1;
R.LON2=LON2;
R.LONINCR=LONINCR;
R.Fields=Fields;
R.Coordinates={'latitude','longitude','time'};
R.Platform=Platform;
R.Pass=Pass;
R.TAvg=TAvg;
R.HRes=HRes;
R.ApplyLandMask=ApplyLandMask;
R.ApplyQualMask=ApplyQualMask;
R.CalcAnomaly=CalcAnomaly;
R.SaveWorkspace=SaveWorkspace;
R.SaveMode=SaveMode;
R.SaveFiles=SaveFiles;
R.DIRNAME=DIRNAME;
R.FNPREFIX=FNPREFIX;
%R.requestDate=datestr(now);
request=R;
handles.request=request;
   
% Load Last request
function pushbuttonLoadLastRequest_Callback(hObject, eventdata, handles)
pname = strrep(mfilename('fullpath'),mfilename,'');
fname = [mfilename,'_request','.mat'];
if exist([pname,fname],'file')
load([pname,fname])
else
msgbox(['File ',pname,fname,' does not exist.']);
end   

R=request;
DATE1=R.DATE1;
DATE2=R.DATE2;
DATEINCR=R.DATEINCR;
LAT1=R.LAT1;
LAT2=R.LAT2;
LATINCR=R.LATINCR;
LON1=R.LON1;
LON2=R.LON2;
LONINCR=R.LONINCR;
Fields=R.Fields;
Platform=R.Platform;
Pass=R.Pass;
TAvg=R.TAvg;
HRes=R.HRes;
ApplyLandMask=R.ApplyLandMask;
ApplyQualMask=R.ApplyQualMask;
CalcAnomaly=R.CalcAnomaly;
SaveWorkspace=R.SaveWorkspace;
SaveMode=R.SaveMode;
SaveFiles=R.SaveFiles;
DIRNAME=R.DIRNAME;
FNPREFIX=R.FNPREFIX;

% select only checked variables
for k=1:length(Fields)
    FIELD=Fields{k};
    if strcmp(FIELD,'sst') set(handles.checkbox1,'Value',1); end
    if strcmp(FIELD,'sst4') set(handles.checkbox2,'Value',1); end
    if strcmp(FIELD,'sst_qual') set(handles.checkbox3,'Value',1); end
    if strcmp(FIELD,'sst4_qual') set(handles.checkbox4,'Value',1); end
end

% get strings from GUI menu and convert them to numbers
set(handles.editTime1yyyy,'String',DATE1(1:4));
set(handles.editTime1mm,'String',DATE1(5:6)); 
set(handles.editTime1dd,'String',DATE1(7:8)); 
set(handles.editTime2yyyy,'String',DATE2(1:4));
set(handles.editTime2mm,'String',DATE2(5:6)); 
set(handles.editTime2dd,'String',DATE2(7:8)); 
set(handles.editLat1,'String',num2str(LAT1));
set(handles.editLat2,'String',num2str(LAT2));
set(handles.editLon1,'String',num2str(LON1));
set(handles.editLon2,'String',num2str(LON2));
set(handles.editLatStep,'String',num2str(LATINCR));
set(handles.editLonStep,'String',num2str(LONINCR));
set(handles.editDirectoryName,'String',DIRNAME);
set(handles.editFileNamePrefix,'String',FNPREFIX);
set(handles.editTimeStep,'String',num2str(DATEINCR));

S=get(handles.popupmenuTemporalAverages,'String');
j=find(strcmp(S,TAvg));
set(handles.popupmenuTemporalAverages,'Value',j);

S=get(handles.popupmenuSpatialResolution,'String');
j=find(strcmp(S,HRes));
set(handles.popupmenuSpatialResolution,'Value',j);

set(handles.radiobuttonPlatformTerra,'Value',0);
set(handles.radiobuttonPlatformAqua,'Value',0);
set(handles.radiobuttonPassDay,'Value',0);
set(handles.radiobuttonPassNight,'Value',0);
if strcmp(Platform,'terra') set(handles.radiobuttonPlatformTerra,'Value',1), set(handles.radiobuttonPlatformAqua,'Value',0); end
if strcmp(Platform,'aqua')  set(handles.radiobuttonPlatformAqua,'Value',1), set(handles.radiobuttonPlatformTerra,'Value',0); end
if strcmp(Pass,'all') set(handles.radiobuttonPassDay,'Value',1);set(handles.radiobuttonPassNight,'Value',1); end
if strcmp(Pass,'day') set(handles.radiobuttonPassDay,'Value',1);set(handles.radiobuttonPassNight,'Value',0); end
if strcmp(Pass,'night') set(handles.radiobuttonPassDay,'Value',0);set(handles.radiobuttonPassNight,'Value',1); end

set(handles.radiobuttonInWksp,'Value',0);
if SaveWorkspace=='y'
    set(handles.radiobuttonInWksp,'Value',1); 
end
set(handles.radiobuttonToFiles,'Value',0);
if SaveFiles=='y'
    set(handles.radiobuttonToFiles,'Value',1); 
end

S=get(handles.popupmenuMode,'String');
j=find(strcmp(S,SaveMode));
if isempty(j) j=1; end
set(handles.popupmenuMode,'Value',j);

% disable saving to workspace for netcdf mode
set(handles.radiobuttonInWksp,'Enable','On');
if strcmp(R.SaveMode,'netcdf')
    set(handles.radiobuttonInWksp,'Value',0,'Enable','Off');
    set(handles.radiobuttonToFiles,'Value',1);
end

SaveFilesOnOff(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)

function SaveFilesOnOff(hObject, eventdata, handles)   
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CheckRanges(hObject, eventdata, handles)
% depending on the checked GUI buttons update appearence of GUI menu
R.Platform='aqua';
if get(handles.radiobuttonPlatformTerra,'Value') == 1 R.Platform='terra'; end % terra
if get(handles.radiobuttonPlatformAqua,'Value') == 1 R.Platform='aqua'; end % aqua
S=get(handles.popupmenuTemporalAverages,'String');j=get(handles.popupmenuTemporalAverages,'Value');R.TAvg=deblank(S{j});
S=get(handles.popupmenuSpatialResolution,'String');j=get(handles.popupmenuSpatialResolution,'Value');R.HRes=deblank(S{j});
if get(handles.radiobuttonPassDay,'Value') == 1 & get(handles.radiobuttonPassNight,'Value') == 0 
    R.Pass='day';  % day
elseif get(handles.radiobuttonPassDay,'Value') == 0 & get(handles.radiobuttonPassNight,'Value') == 1 
    R.Pass='night'; % night
else
    R.Pass='all';
end
R.DataSetBranch=['MODIS_' R.Platform '_' R.HRes '_' R.TAvg '_' R.Pass]; 

%example: load DSI_OceanColor_MODIS_4km_8day_hdf.mat
%S=['load dsi_modis.mat DSI_' R.DataSetBranch];
%eval(['load dsi_modis.mat DSI_' R.DataSetBranch])
eval(['DSI=handles.DSI.DSI_' R.DataSetBranch ';'])
%eval(['clear DSI_' R.DataSetBranch ';'])

set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2,'yyyy-mm-dd') ']'])
set(handles.textLatitudeRange,'String',...
['Available latitude range: [' num2str(DSI.Range.Latitude1) ' to ' num2str(DSI.Range.Latitude2) ']'])
set(handles.textLongitudeRange,'String',...
['Available longitude range: [' num2str(DSI.Range.Longitude1) ' to ' num2str(DSI.Range.Longitude2) ']'])

