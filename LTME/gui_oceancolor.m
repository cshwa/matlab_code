function varargout = gui_oceancolor(varargin)
% Boots-up OceanColor GUI, interacts with user via interface.
% gui_oceancolor M-file for gui_oceancolor.fig
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
% Copyright 2007, 2008, 2009
% $Version 3.0.0$

%==========================================================================
% Meri Sheremet
% May 2009
%
% REVISION HISTORY:
% 2007/11/20 created
% 2008/06/15 modified for compiled version, ceb
%==========================================================================

% Last Modified by GUIDE v2.5 04-Feb-2009 07:25:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_oceancolor_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_oceancolor_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_oceancolor is made visible.
function gui_oceancolor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_oceancolor (see VARARGIN)

DSI=load('dsi_oceancolor.mat');
handles.DSI=DSI;

% Set interface name.
filename = 'get_oceancolor';
txt = ['NASA/REASoN Ocean Data Portal - '];
tmp1 = ['Ocean Color GUI'];
tmp2 = '';
txt = [txt,tmp1];
if (~isdeployed)
    tmp2 = strrep(op_version('mfile',filename),'Version ','');
    txt = [txt,' (Ver ',tmp2,')'];
end
set(hObject,'Name',txt);

% Place logos on the GUI with links
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
set(handles.pushbutton7,'CData',img)
pos = get(handles.pushbutton7,'position'); %units will be in characters
set(handles.pushbutton7,'position',[pos(1),pos(2),CX,CY]);  

[img,MAP]=imread('logo_nasa1.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonNASA,'CData',img)
pos = get(handles.pushbuttonNASA,'position'); %units will be in characters
set(handles.pushbuttonNASA,'position',[pos(1),pos(2),CX,CY]);  

[img,MAP]=imread('logo_nasa_gsfc.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonOceanColor,'CData',img)
pos = get(handles.pushbuttonOceanColor,'position'); %units will be in characters
set(handles.pushbuttonOceanColor,'position',[pos(1),pos(2),CX,CY]);  

[img,MAP]=imread('logo_csiro.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonCSIRO,'CData',img)
pos = get(handles.pushbuttonCSIRO,'position'); %units will be in characters
set(handles.pushbuttonCSIRO,'position',[pos(1),pos(2),CX,CY]);  

set(handles.editTimeStep,'String',num2str(1));
set(handles.editLatStep,'String',num2str(1));
set(handles.editLonStep,'String',num2str(1));
%set(handles.editFileNamePrefix,'String','opendap_');

txt=['Help'];
op_tooltipwrap(handles.pushbuttonInfo,txt,50);

[img,MAP]=imread('logo_questionmark.jpg');
% [NY,NX,N3]=size(img);
% CX=(NX+8)/ppcx;
% CY=(NY+12)/ppcy;
% set(handles.pushbuttonInfo,'CData',img)
% pos = get(handles.pushbuttonInfo,'position'); %units will be in characters
% set(handles.pushbuttonInfo,'position',[pos(1),pos(2),CX,CY]);
[N1,N2,N3]=size(img);
set(handles.pushbuttonInfo,'string','')
set(handles.pushbuttonInfo,'CData',img)
set(handles.pushbuttonInfo,'Units','pixels')
pos = get(handles.pushbuttonInfo,'position');
set(handles.pushbuttonInfo,'position',[pos(1) pos(2) N1 N2])

if get(handles.radiobuttonCZCS,'Value')==0 & get(handles.radiobuttonMODIS,'Value')==0 ...
 & get(handles.radiobuttonSeaWiFS,'Value')==0 & get(handles.radiobuttonMergedMODISSeaWiFS,'Value')==0 
    set(handles.radiobuttonCZCS,'Value',1)
end
    
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

% Make 'Load Last Request' button inactive when file gui_xxx_request.mat is
% absent.
if ~exist('gui_oceancolor_request.mat','file') 
    set(handles.pushbuttonLoadLastRequest,'Enable','Off')
end

set(handles.checkbox210,'Visible','off','Enable','off') % disable temporary
CheckRanges(hObject, eventdata, handles)
% Choose default command line output for gui_oceancolor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_oceancolor wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_oceancolor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function radiobuttonCZCS_Callback(hObject, eventdata, handles)
%CZCS
CheckExclusive(hObject, eventdata, handles,{'radiobuttonCZCS'},{'radiobuttonMODIS';'radiobuttonSeaWiFS';'radiobuttonMergedMODISSeaWiFS'})
set(handles.popupmenuTemporalAverages,'String',{'8day';'monthly';},'Value',1);
set(handles.popupmenuSpatialResolution,'String',{'9km';},'Value',1);
CheckRanges(hObject, eventdata, handles)

function radiobuttonMODIS_Callback(hObject, eventdata, handles)
% MODIS
CheckExclusive(hObject, eventdata, handles,{'radiobuttonMODIS'},{'radiobuttonCZCS';'radiobuttonSeaWiFS';'radiobuttonMergedMODISSeaWiFS'})
set(handles.popupmenuTemporalAverages,'String',{'3day';'8day';},'Value',1);
set(handles.popupmenuSpatialResolution,'String',{'4km';'9km';},'Value',1);
CheckRanges(hObject, eventdata, handles)

function radiobuttonSeaWiFS_Callback(hObject, eventdata, handles)
% SeaWiFS
CheckExclusive(hObject, eventdata, handles,{'radiobuttonSeaWiFS'},{'radiobuttonMODIS';'radiobuttonCZCS';'radiobuttonMergedMODISSeaWiFS'})
set(handles.popupmenuTemporalAverages,'String',{'8day';},'Value',1);
set(handles.popupmenuSpatialResolution,'String',{'9km';},'Value',1);
CheckRanges(hObject, eventdata, handles)

function radiobuttonMergedMODISSeaWiFS_Callback(hObject, eventdata, handles)
% Merged MODIS & SeaWiFS
CheckExclusive(hObject, eventdata, handles,{'radiobuttonMergedMODISSeaWiFS'},{'radiobuttonMODIS';'radiobuttonSeaWiFS';'radiobuttonCZCS'})
set(handles.popupmenuTemporalAverages,'String',{'1day';'8day';},'Value',1);
set(handles.popupmenuSpatialResolution,'String',{'9km';},'Value',1);
CheckRanges(hObject, eventdata, handles)

function radiobuttonInWksp_Callback(hObject, eventdata, handles)
function radiobuttonToFiles_Callback(hObject, eventdata, handles)
SaveFilesOnOff(hObject, eventdata, handles)

function checkbox101_Callback(hObject, eventdata, handles)
function checkbox102_Callback(hObject, eventdata, handles)
function checkbox201_Callback(hObject, eventdata, handles)
function checkbox202_Callback(hObject, eventdata, handles)
function checkbox203_Callback(hObject, eventdata, handles)
function checkbox204_Callback(hObject, eventdata, handles)
function checkbox205_Callback(hObject, eventdata, handles)
function checkbox206_Callback(hObject, eventdata, handles)
function checkbox207_Callback(hObject, eventdata, handles)
function checkbox208_Callback(hObject, eventdata, handles)
function checkbox209_Callback(hObject, eventdata, handles)
function checkbox210_Callback(hObject, eventdata, handles)
function checkbox301_Callback(hObject, eventdata, handles)
function checkbox302_Callback(hObject, eventdata, handles)
function checkbox303_Callback(hObject, eventdata, handles)
function checkbox304_Callback(hObject, eventdata, handles)
function checkbox305_Callback(hObject, eventdata, handles)
function checkbox306_Callback(hObject, eventdata, handles)
function checkbox307_Callback(hObject, eventdata, handles)
function checkbox401_Callback(hObject, eventdata, handles)

function popupmenuTemporalAverages_Callback(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)
S=get(handles.popupmenuTemporalAverages,'String');j=get(handles.popupmenuTemporalAverages,'Value');S=deblank(S{j});
TXT=['Temporal resolution is ' S];
set(handles.textTemporalResolution,'String',TXT,'ForegroundColor',[0 0 1])
function popupmenuTemporalAverages_CreateFcn(hObject, eventdata, handles)
function popupmenuSpatialResolution_Callback(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)
S=get(handles.popupmenuSpatialResolution,'String');j=get(handles.popupmenuSpatialResolution,'Value');S=deblank(S{j});
TXT=['Spatial resolution is ' S];
set(handles.textSpatialResolution,'String',TXT,'ForegroundColor',[0 0 1])
function popupmenuSpatialResolution_CreateFcn(hObject, eventdata, handles)
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

function editTime1yyyy_Callback(hObject, eventdata, handles)
function editTime1yyyy_CreateFcn(hObject, eventdata, handles)
function editTime1mm_Callback(hObject, eventdata, handles)
function editTime1mm_CreateFcn(hObject, eventdata, handles)
function editTime1dd_Callback(hObject, eventdata, handles)
function editTime1dd_CreateFcn(hObject, eventdata, handles)
function editTime2yyyy_Callback(hObject, eventdata, handles)
function editTime2yyyy_CreateFcn(hObject, eventdata, handles)
function editTime2mm_Callback(hObject, eventdata, handles)
function editTime2mm_CreateFcn(hObject, eventdata, handles)
function editTime2dd_Callback(hObject, eventdata, handles)
function editTime2dd_CreateFcn(hObject, eventdata, handles)
function editLat1_Callback(hObject, eventdata, handles)
function editLat1_CreateFcn(hObject, eventdata, handles)
function editLat2_Callback(hObject, eventdata, handles)
function editLat2_CreateFcn(hObject, eventdata, handles)
function editLon1_Callback(hObject, eventdata, handles)
function editLon1_CreateFcn(hObject, eventdata, handles)
function editLon2_Callback(hObject, eventdata, handles)
function editLon2_CreateFcn(hObject, eventdata, handles)

function editTimeStep_Callback(hObject, eventdata, handles)
function editTimeStep_CreateFcn(hObject, eventdata, handles)
function editLatStep_Callback(hObject, eventdata, handles)
function editLatStep_CreateFcn(hObject, eventdata, handles)
function editLonStep_Callback(hObject, eventdata, handles)
function editLonStep_CreateFcn(hObject, eventdata, handles)

function editDirectoryName_Callback(hObject, eventdata, handles)
function editDirectoryName_CreateFcn(hObject, eventdata, handles)
function editFileNamePrefix_Callback(hObject, eventdata, handles)
function editFileNamePrefix_CreateFcn(hObject, eventdata, handles)

function pushbuttonGetData_Callback(hObject, eventdata, handles)
% Main button that starts data acquisition.

% Turn off ability to re-issue "Get Data."
set(hObject,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off','String','Please wait...');
pause(.1)

% Get input from user.
request=InputGet(hObject, eventdata, handles);

% Save request to .mat file.
filename = 'gui_oceancolor_request.mat';
%if (~isdeployed)
    pname = mfilename('fullpath');
    pname = pname(1:length(pname)-length(mfilename)); %mfilename does not have .m in it
    filename = [pname,filename];
%end
save(filename,'request','-mat')
set(handles.pushbuttonLoadLastRequest,'Enable','On')

R=request;

%    if isempty(R.LAT1) 
%        set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25])
%        return; 
%    end

% If it passes inspection, send down to user's workspace
% as variable "request". Then the user can get the data
%assignin('base','request',request);

get_oceancolor(request);

set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25])

function pushbuttonBrowse_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonToFiles,'Value') == 1
DIRNAME=uigetdir(pwd,'Select a Directory for Downloading Files');
set(handles.editDirectoryName,'String',DIRNAME)
end
function pushbuttonOceanColor_Callback(hObject, eventdata, handles)
web('http://oceancolor.gsfc.nasa.gov/');
%web 'http://oceancolor.gsfc.nasa.gov/' -browser
function pushbutton7_Callback(hObject, eventdata, handles)
web('http://www.opendap.org/');
%web 'http://www.opendap.org/' -browser
function pushbuttonUpdateInventory_Callback(hObject, eventdata, handles)
dsi_oceancolor
function pushbuttonCSIRO_Callback(hObject, eventdata, handles)
web('http://www.marine.csiro.au/remotesensing/restricted/un001.html');
%web 'http://www.marine.csiro.au/remotesensing/restricted/un001.html' -browser 
function pushbuttonNASA_Callback(hObject, eventdata, handles)
web('http://www.nasa.gov');
%web 'http://www.nasa.gov' -browser
%web 'http://www.nasa.gov/topics/earth/index.html' -browser


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get Input from GUI menu
function request=InputGet(hObject, eventdata, handles)

S1=deblank(get(handles.editTime1yyyy,'String'));
S2=deblank(get(handles.editTime1mm,'String'));
S3=deblank(get(handles.editTime1dd,'String'));
S4=deblank(get(handles.editTime2yyyy,'String'));
S5=deblank(get(handles.editTime2mm,'String'));
S6=deblank(get(handles.editTime2dd,'String'));

R.DATE1='00000000';
R.DATE1(4-length(S1)+1:4)=S1;
R.DATE1(6-length(S2)+1:6)=S2;
R.DATE1(8-length(S3)+1:8)=S3;
R.DATE2='00000000';
R.DATE2(4-length(S4)+1:4)=S4;
R.DATE2(6-length(S5)+1:6)=S5;
R.DATE2(8-length(S6)+1:8)=S6;

R.LAT1=str2num(deblank(get(handles.editLat1,'String')));
R.LAT2=str2num(deblank(get(handles.editLat2,'String')));
R.LON1=str2num(deblank(get(handles.editLon1,'String')));
R.LON2=str2num(deblank(get(handles.editLon2,'String')));

R.DATEINCR=str2num(deblank(get(handles.editTimeStep,'String')));
R.LATINCR=str2num(deblank(get(handles.editLatStep,'String')));
R.LONINCR=str2num(deblank(get(handles.editLonStep,'String')));

if get(handles.radiobuttonCZCS,'Value') == 1 DSN='OceanColor_CZCS';               DSNum=1; end
if get(handles.radiobuttonMODIS,'Value') == 1 DSN='OceanColor_MODIS';              DSNum=2; end
if get(handles.radiobuttonSeaWiFS,'Value') == 1 DSN='OceanColor_SeaWiFS';            DSNum=3; end
if get(handles.radiobuttonMergedMODISSeaWiFS,'Value') == 1 DSN='OceanColor_MergedMODISSeaWiFS'; DSNum=4; end

S=get(handles.popupmenuTemporalAverages,'String');
j=get(handles.popupmenuTemporalAverages,'Value');
R.TAvg=deblank(S{j}); % only one value
S=get(handles.popupmenuSpatialResolution,'String');
j=get(handles.popupmenuSpatialResolution,'Value');
R.HRes=deblank(S{j}); % only one value

R.DataSetBranch=[DSN '_' R.HRes '_' R.TAvg '_' 'hdf'];
%example: load DSI_OceanColor_MODIS_4km_8day_hdf.mat
%eval(['load dsi_oceancolor.mat DSI_' R.DataSetBranch])
eval(['DSI=handles.DSI.DSI_' R.DataSetBranch ';'])
%eval(['clear DSI_' R.DataSetBranch ';'])

R.Fields={};
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
    CBNUM=num2str(100*DSNum+k);
    if eval(['get(handles.checkbox' CBNUM ',''Value'')==1']) R.Fields{end+1}=FIELD; end
end

R.Coordinates={'latitude','longitude','time'};

SaveModeL={'opendap';'native';'netcdf'};
j=get(handles.popupmenuMode,'Value');
R.SaveMode=SaveModeL{j};

R.SaveWorkspace='n'; 
if get(handles.radiobuttonInWksp,'Value') == 1 R.SaveWorkspace='y'; end
R.SaveFiles='n'; 
if get(handles.radiobuttonToFiles,'Value') == 1 R.SaveFiles='y'; end
R.DIRNAME =deblank(get(handles.editDirectoryName,'String'));
if isempty(R.DIRNAME) R.DIRNAME=''; end
R.FNPREFIX=deblank(get(handles.editFileNamePrefix,'String'));

%% Get default command line output from handles structure.   
%R.RequestDate=datestr(now,'yyyy-mm-dd HH:MM:SS');
request=R;

%handles.output = hObject;
%guidata(hObject, handles);

% Load Last request
function pushbuttonLoadLastRequest_Callback(hObject, eventdata, handles)
pname = strrep(mfilename('fullpath'),mfilename,'');
fname = [mfilename,'_request','.mat'];
if exist([pname,fname],'file')
load([pname,fname]);
else
msgbox(['File ',pname,fname,' does not exist.']);
return
end
if ~exist('request','var')
    request = Request;
end
R=request;

set(handles.editTime1yyyy,'String',R.DATE1(1:4));
set(handles.editTime1mm,'String',R.DATE1(5:6)); 
set(handles.editTime1dd,'String',R.DATE1(7:8)); 
set(handles.editTime2yyyy,'String',R.DATE2(1:4));
set(handles.editTime2mm,'String',R.DATE2(5:6)); 
set(handles.editTime2dd,'String',R.DATE2(7:8)); 
set(handles.editLat1,'String',num2str(R.LAT1));
set(handles.editLat2,'String',num2str(R.LAT2));
set(handles.editLon1,'String',num2str(R.LON1));
set(handles.editLon2,'String',num2str(R.LON2));
set(handles.editTimeStep,'String',num2str(R.DATEINCR));
set(handles.editLatStep,'String',num2str(R.LATINCR));
set(handles.editLonStep,'String',num2str(R.LONINCR));

if strcmp(R.DataSetBranch(12:15),'CZCS')
CheckExclusive(hObject, eventdata, handles,{'radiobuttonCZCS'},{'radiobuttonMODIS';'radiobuttonSeaWiFS';'radiobuttonMergedMODISSeaWiFS'})
DSNum=1; TAvgL={'8day';'monthly';};HResL={'9km';}; 
end 
if strcmp(R.DataSetBranch(12:16),'MODIS')    
CheckExclusive(hObject, eventdata, handles,{'radiobuttonMODIS'},{'radiobuttonCZCS';'radiobuttonSeaWiFS';'radiobuttonMergedMODISSeaWiFS'})
DSNum=2; TAvgL={'3day';'8day';};   HResL={'4km';'9km';}; 
end 
if strcmp(R.DataSetBranch(12:18),'SeaWiFS')  
CheckExclusive(hObject, eventdata, handles,{'radiobuttonSeaWiFS'},{'radiobuttonMODIS';'radiobuttonCZCS';'radiobuttonMergedMODISSeaWiFS'})
DSNum=3; TAvgL={'8day';};          HResL={'9km';}; 
end 
if strcmp(R.DataSetBranch(12:17),'Merged')   
CheckExclusive(hObject, eventdata, handles,{'radiobuttonMergedMODISSeaWiFS'},{'radiobuttonMODIS';'radiobuttonSeaWiFS';'radiobuttonCZCS'})
DSNum=4; TAvgL={'1day';'8day';};   HResL={'9km';}; 
end 

set(handles.popupmenuTemporalAverages,'String',TAvgL);
j=find(strcmp(TAvgL,R.TAvg));
if isempty(j) j=1; end
set(handles.popupmenuTemporalAverages,'Value',j);
set(handles.popupmenuSpatialResolution,'String',HResL);
j=find(strcmp(HResL,R.HRes));
if isempty(j) j=1; end
set(handles.popupmenuSpatialResolution,'Value',j);

%example: load DSI_OceanColor_MODIS_4km_8day_hdf.mat
%S=['load dsi_oceancolor.mat DSI_' R.DataSetBranch];
%eval(['load dsi_oceancolor.mat DSI_' R.DataSetBranch])
eval(['DSI=handles.DSI.DSI_' R.DataSetBranch ';'])
%eval(['clear DSI_' R.DataSetBranch ';'])
set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2,'yyyy-mm-dd') ']'])
set(handles.textLatitudeRange,'String',...
['Available latitude range: [' num2str(DSI.Range.Latitude1) ' to ' num2str(DSI.Range.Latitude2) ']'])
set(handles.textLongitudeRange,'String',...
['Available longitude range: [' num2str(DSI.Range.Longitude1) ' to ' num2str(DSI.Range.Longitude2) ']'])

for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
    CBNUM=num2str(100*DSNum+k);
    eval(['set(handles.checkbox' CBNUM ',''Value'',0)'])
    if strmatch(FIELD,R.Fields) eval(['set(handles.checkbox' CBNUM ',''Value'',1)']) ; end
end

CheckRanges(hObject, eventdata, handles)

SaveModeL={'opendap';'native';'netcdf'};
set(handles.popupmenuMode,'String',SaveModeL);
j=find(strcmp(SaveModeL,R.SaveMode));
if isempty(j) j=1; end
set(handles.popupmenuMode,'Value',j);

if R.SaveWorkspace=='y' set(handles.radiobuttonInWksp,'Value',1); end
if R.SaveFiles=='y' set(handles.radiobuttonToFiles,'Value',1); end
set(handles.editDirectoryName,'String',R.DIRNAME);
set(handles.editFileNamePrefix,'String',R.FNPREFIX);
S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');
R.SaveMode=S{j};
% disable saving to workspace for netcdf mode
set(handles.radiobuttonInWksp,'Enable','On');
if strcmp(R.SaveMode,'netcdf')
    set(handles.radiobuttonInWksp,'Value',0,'Enable','Off');
    set(handles.radiobuttonToFiles,'Value',1);
end
SaveFilesOnOff(hObject, eventdata, handles)



function CheckExclusive(hObject, eventdata, handles, Us, Them)
% sets value of Us to 1 and Them to 0
% example: Us={'radiobuttonCZCS'}; Them={'radiobuttonMODIS';'radiobuttonSeaWiFS'}
eval(['set(handles.' Us{1} ',''Value'',1)'])
for k=1:length(Them)
eval(['set(handles.' Them{k} ',''Value'',0)'])
end

function SaveFilesOnOff(hObject, eventdata, handles)   
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end

function CheckRanges(hObject, eventdata, handles)

DSNum=0;
if get(handles.radiobuttonCZCS,'Value') == 1 
    DSN='OceanColor_CZCS';
    DSNum=1; TAvgL={'8day';'monthly';};HResL={'9km';};  
end
if get(handles.radiobuttonMODIS,'Value') == 1 
    DSN='OceanColor_MODIS';
    DSNum=2; TAvgL={'3day';'8day';};   HResL={'4km';'9km';};                
end
if get(handles.radiobuttonSeaWiFS,'Value') == 1 
    DSN='OceanColor_SeaWiFS';
    DSNum=3; TAvgL={'8day';};          HResL={'9km';};  
end
if get(handles.radiobuttonMergedMODISSeaWiFS,'Value') == 1 
    DSN='OceanColor_MergedMODISSeaWiFS'; 
    DSNum=4; TAvgL={'1day';'8day';};   HResL={'9km';};  
end
if DSNum==0 return; end

set(handles.popupmenuTemporalAverages,'String',TAvgL);
set(handles.popupmenuSpatialResolution,'String',HResL);
   
S=get(handles.popupmenuTemporalAverages,'String');
j=get(handles.popupmenuTemporalAverages,'Value');
if isempty(j) return; end
R.TAvg=deblank(S{j}); % only one value
S=get(handles.popupmenuSpatialResolution,'String');
j=get(handles.popupmenuSpatialResolution,'Value');
if isempty(j) return; end
R.HRes=deblank(S{j}); % only one value

R.DataSetBranch=[DSN '_' R.HRes '_' R.TAvg '_' 'hdf'];
%CheckTLLRanges(hObject, eventdata, handles,R.DataSetBranch)
%example: load DSI_OceanColor_MODIS_4km_8day_hdf.mat
%S=['load dsi_oceancolor.mat DSI_' R.DataSetBranch];
%eval(['load dsi_oceancolor.mat DSI_' R.DataSetBranch])
eval(['DSI=handles.DSI.DSI_' R.DataSetBranch ';'])
%eval(['clear DSI_' R.DataSetBranch ';'])
set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2,'yyyy-mm-dd') ']'])
set(handles.textLatitudeRange,'String',...
['Available latitude range: [' num2str(DSI.Range.Latitude1) ' to ' num2str(DSI.Range.Latitude2) ']'])
set(handles.textLongitudeRange,'String',...
['Available longitude range: [' num2str(DSI.Range.Longitude1) ' to ' num2str(DSI.Range.Longitude2) ']'])

% make visible and active only checkboxes corresponding to existing
% variables
L=[(101:102) (201:210) (301:307) (401)]; % all field checkboxes
for k=1:length(L)
CBNUM=num2str(L(k));
%eval(['set(handles.checkbox' CBNUM ',''Visible'',''off'',''Enable'',''off'')'])
eval(['set(handles.checkbox' CBNUM ',''Enable'',''off'')'])
end
for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
    CBNUM=num2str(100*DSNum+k);
%    eval(['set(handles.checkbox' CBNUM ',''Visible'',''on'',''Enable'',''on'',''String'',''' FIELD ''')'])
    eval(['set(handles.checkbox' CBNUM ',''Enable'',''on'',''String'',''' FIELD ''')'])
%    if R.Fields{k,end}=='y' eval(['set(handles.checkbox' CBNUM ',''Value'',1)']) ; end
end

S=get(handles.popupmenuTemporalAverages,'String');j=get(handles.popupmenuTemporalAverages,'Value');S=deblank(S{j});
TXT=['Temporal resolution is ' S];
set(handles.textTemporalResolution,'String',TXT,'ForegroundColor',[0 0 1])
S=get(handles.popupmenuSpatialResolution,'String');j=get(handles.popupmenuSpatialResolution,'Value');S=deblank(S{j});
TXT=['Spatial resolution is ' S];
set(handles.textSpatialResolution,'String',TXT,'ForegroundColor',[0 0 1])

function pushbuttonInfo_Callback(hObject, eventdata, handles)
h = op_timemessage('polar',0);