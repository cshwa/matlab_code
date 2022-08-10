function varargout = gui_hycom(varargin)
% Boots-up HYCOM GUI, interacts with user via interface.
% gui_hycom M-file for gui_hycom.fig
%
% INPUT(S):
%
%   [none] -- DataSetDATES is input to GUI to avoid using globals.
%
% OUTPUT(S):
%
%   h -- (optional) GUI handle
%
% OPeNDAP Science Team
% Copyright 2007, 2008
% $Version 3.0.0$

%warning('off','MATLAB:dispatcher:InexactMatch');
%warning('on','MATLAB:dispatcher:InexactMatch');

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_hycom_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_hycom_OutputFcn, ...
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



% --- Executes just before gui_hycom is made visible.
function gui_hycom_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_hycom (see VARARGIN)

DSI=load('dsi_hycom.mat');
handles.DSI=DSI;

% Set interface name.
filename = 'get_hycom';
txt = [''];
tmp1 = ['HYCOM GUI'];
tmp2 = '';
txt = [txt,tmp1];
if (~isdeployed)
    tmp2 = strrep(op_version('mfile',filename),'Version ','');
    txt = [txt,' (Ver ',tmp2,')'];
end
set(hObject,'Name',txt);

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

[img,MAP]=imread('logo_nopp.jpg');
[NY,NX,N3]=size(img);
CX=(NX+4)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+4)/ppcy;
set(handles.pushbuttonNOPP,'CData',img)
pos = get(handles.pushbuttonNOPP,'position'); %units will be in characters
set(handles.pushbuttonNOPP,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_hycom.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;
CY=(NY+12)/ppcy;
set(handles.pushbuttonHycom,'CData',img)
pos = get(handles.pushbuttonHycom,'position'); %units will be in characters
set(handles.pushbuttonHycom,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_coaps.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;
CY=(NY+12)/ppcy;
set(handles.pushbuttonHycomCoaps,'CData',img)
pos = get(handles.pushbuttonHycomCoaps,'position'); %units will be in characters
set(handles.pushbuttonHycomCoaps,'position',[pos(1),pos(2),CX,CY]);

set(handles.editTimeStep,'String',num2str(1));
set(handles.editLatStep,'String',num2str(1));
set(handles.editLonStep,'String',num2str(1));
set(handles.editDepthStep,'String',num2str(1));
%set(handles.editFileNamePrefix,'String','opendap_');

set(handles.popupmenuDatasetBranch,'Value',2);
CheckRanges(hObject, eventdata, handles)

LevelsOnOff(hObject, eventdata, handles)

%SaveFilesOnOff(hObject, eventdata, handles)
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

set(handles.pushbuttonGetData,'String','Get Data','BackgroundColor',[0 1. 0.25])

% Make 'Load Last Request' button inactive when file gui_hycom_request.mat is
% absent
pname = strrep(mfilename('fullpath'),mfilename,'');
fname = [mfilename,'_request','.mat'];
%'gui_hycom_request.mat';
if ~exist([pname,fname],'file') 
    set (handles.pushbuttonLoadLastRequest,'Enable','Off')
end

% Update handles "openfigure" field.
handles.openfigure = 1;

% Choose default command line output for gui_hycom
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_hycom wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_hycom_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

function checkbox1_Callback(hObject, eventdata, handles)
LevelsOnOff(hObject, eventdata, handles)
function checkbox2_Callback(hObject, eventdata, handles)
LevelsOnOff(hObject, eventdata, handles)
function checkbox3_Callback(hObject, eventdata, handles)
LevelsOnOff(hObject, eventdata, handles)
function checkbox4_Callback(hObject, eventdata, handles)
LevelsOnOff(hObject, eventdata, handles)
function checkbox5_Callback(hObject, eventdata, handles)
LevelsOnOff(hObject, eventdata, handles)

function checkbox6_Callback(hObject, eventdata, handles)
function checkbox7_Callback(hObject, eventdata, handles)
function checkbox8_Callback(hObject, eventdata, handles)
function checkbox9_Callback(hObject, eventdata, handles)
function checkbox10_Callback(hObject, eventdata, handles)
function checkbox11_Callback(hObject, eventdata, handles)
function checkbox12_Callback(hObject, eventdata, handles)
function checkbox13_Callback(hObject, eventdata, handles)
function checkbox14_Callback(hObject, eventdata, handles)
function checkbox15_Callback(hObject, eventdata, handles)
function checkbox16_Callback(hObject, eventdata, handles)
function checkbox17_Callback(hObject, eventdata, handles)

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
function editDepth1_Callback(hObject, eventdata, handles)
function editDepth1_CreateFcn(hObject, eventdata, handles)
function editDepth2_Callback(hObject, eventdata, handles)
function editDepth2_CreateFcn(hObject, eventdata, handles)
function editTimeStep_Callback(hObject, eventdata, handles)
function editTimeStep_CreateFcn(hObject, eventdata, handles)
function editLatStep_Callback(hObject, eventdata, handles)
function editLatStep_CreateFcn(hObject, eventdata, handles)
function editLonStep_Callback(hObject, eventdata, handles)
function editLonStep_CreateFcn(hObject, eventdata, handles)
function editDepthStep_Callback(hObject, eventdata, handles)
function editDepthStep_CreateFcn(hObject, eventdata, handles)
function editDirectoryName_Callback(hObject, eventdata, handles)
function editDirectoryName_CreateFcn(hObject, eventdata, handles)
function editFileNamePrefix_Callback(hObject, eventdata, handles)
function editFileNamePrefix_CreateFcn(hObject, eventdata, handles)

function radiobuttonInWksp_Callback(hObject, eventdata, handles)
function radiobuttonToFiles_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end

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
%SaveFilesOnOff(hObject, eventdata, handles)

function pushbuttonSelectAll4D_Callback(hObject, eventdata, handles)
S=get(handles.pushbuttonSelectAll4D,'String');
if S(1) == 'S' % Select All
set(handles.checkbox1,'Value',1)
set(handles.checkbox2,'Value',1)
set(handles.checkbox3,'Value',1)
set(handles.checkbox4,'Value',1)
set(handles.checkbox5,'Value',1)
set(handles.pushbuttonSelectAll4D,'String','Clear All')
elseif S(1) == 'C' % Clear All
set(handles.checkbox1,'Value',0)
set(handles.checkbox2,'Value',0)
set(handles.checkbox3,'Value',0)
set(handles.checkbox4,'Value',0)
set(handles.checkbox5,'Value',0)
set(handles.pushbuttonSelectAll4D,'String','Select All')
end  
LevelsOnOff(hObject, eventdata, handles)

function pushbuttonSelectAll3D1_Callback(hObject, eventdata, handles)
S=get(handles.pushbuttonSelectAll3D1,'String');
if S(1) == 'S' % Select All
set(handles.checkbox6,'Value',1)
set(handles.checkbox7,'Value',1) 
set(handles.checkbox8,'Value',1)
set(handles.checkbox9,'Value',1)
set(handles.checkbox10,'Value',1)
set(handles.checkbox11,'Value',1)
set(handles.pushbuttonSelectAll3D1,'String','Clear All')
elseif S(1) == 'C' % Clear All
set(handles.checkbox6,'Value',0)
set(handles.checkbox7,'Value',0) 
set(handles.checkbox8,'Value',0)
set(handles.checkbox9,'Value',0)
set(handles.checkbox10,'Value',0)
set(handles.checkbox11,'Value',0)
set(handles.pushbuttonSelectAll3D1,'String','Select All')
end 

function pushbuttonSelectAll3D2_Callback(hObject, eventdata, handles)
S=get(handles.pushbuttonSelectAll3D2,'String');
if S(1) == 'S' % Select All
set(handles.checkbox12,'Value',1)
set(handles.checkbox13,'Value',1)
set(handles.checkbox14,'Value',1)
set(handles.checkbox15,'Value',1)
set(handles.checkbox16,'Value',1)
set(handles.checkbox17,'Value',1)
set(handles.pushbuttonSelectAll3D2,'String','Clear All')
elseif S(1) == 'C' % Clear All
set(handles.checkbox12,'Value',0)
set(handles.checkbox13,'Value',0)
set(handles.checkbox14,'Value',0)
set(handles.checkbox15,'Value',0)
set(handles.checkbox16,'Value',0)
set(handles.checkbox17,'Value',0)
set(handles.pushbuttonSelectAll3D2,'String','Select All')
end

function pushbuttonBrowse_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonToFiles,'Value') == 1
DIRNAME=uigetdir(pwd,'Select a Directory for Downloading Files');
set(handles.editDirectoryName,'String',DIRNAME)
end

% The main button that starts data acquisition
function pushbuttonGetData_Callback(hObject, eventdata, handles)
% set(handles.pushbuttonGetData,'String','For a new request reissue get_hycom')
% uiresume(handles.figure1)
%% Turn off ability to re-issue "Get Data"
set(hObject,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off','String','Please wait ...');
pause(.1)

%% GET INPUT FROM GUI 
request=InputGet(hObject, eventdata, handles);

%save gui_hycom_request.mat request -MAT %-V6 %-APPEND
% Remember request for future sessions.
pname = strrep(mfilename('fullpath'),mfilename,'');
fname = [mfilename,'_request','.mat'];
save([pname,fname],'request','-mat')

%'gui_hycom_request.mat';
if exist([pname,fname],'file') 
    set (handles.pushbuttonLoadLastRequest,'Enable','On')
end

% If it passes inspection, send down to user's workspace
% as variable "request". Then the user can get the data
%assignin('base','request',request);

get_hycom(request);

set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25]);

function pushbuttonHycom_Callback(hObject, eventdata, handles)
web('http://hycom.coaps.fsu.edu/data/atlantic_info.html');
%web 'http://hycom.coaps.fsu.edu/data/atlantic_info.html' -browser
%web 'http://oceanmodeling.rsmas.miami.edu/hycom/' -browser
function pushbuttonOpendap_Callback(hObject, eventdata, handles)
web('http://www.opendap.org/');
%web 'http://www.opendap.org/' -browser
function pushbuttonHycomCoaps_Callback(hObject, eventdata, handles)
web('http://hycom.coaps.fsu.edu/thredds/dodsC/atl_ops.info');
%web 'http://hycom.coaps.fsu.edu/thredds/dodsC/atl_ops.info' -browser
function pushbuttonNOPP_Callback(hObject, eventdata, handles)
web('http://www.nopp.org/');
%web 'http://www.nopp.org/' -browser

function request=InputGet(hObject, eventdata, handles)
DSI=LoadDSI(hObject, eventdata, handles);
% select only checked variables
R.Fields={};
for k=1:length(DSI.Fields_CBNum)
FIELD=DSI.Fields{k};    
CBNUM=num2str(DSI.Fields_CBNum{k});
    if eval(['get(handles.checkbox' CBNUM ',''Value'')']) == 1
    R.Fields{end+1}=FIELD;    
    end
end
L4D=get(handles.checkbox1,'Value')|get(handles.checkbox2,'Value')|...
    get(handles.checkbox3,'Value')|get(handles.checkbox4,'Value')|...
    get(handles.checkbox5,'Value');   
if L4D==0
    R.Coordinates={'latitude','longitude','time'};
else
    R.Coordinates={'latitude','longitude','depth','time'};
end


% get strings from GUI menu and convert them to numbers
yyyy=deblank(get(handles.editTime1yyyy,'String'));
mm=deblank(get(handles.editTime1mm,'String')); if length(mm) == 1 mm=['0' mm]; end
dd=deblank(get(handles.editTime1dd,'String')); if length(dd) == 1 dd=['0' dd]; end
R.DATE1=[yyyy mm dd];

yyyy=deblank(get(handles.editTime2yyyy,'String'));
mm=deblank(get(handles.editTime2mm,'String')); if length(mm) == 1 mm=['0' mm]; end
dd=deblank(get(handles.editTime2dd,'String')); if length(dd) == 1 dd=['0' dd]; end
R.DATE2=[yyyy mm dd];

R.LAT1=str2num(deblank(get(handles.editLat1,'String')));
R.LAT2=str2num(deblank(get(handles.editLat2,'String')));
R.LON1=str2num(deblank(get(handles.editLon1,'String')));
R.LON2=str2num(deblank(get(handles.editLon2,'String')));
R.DEPTH1=str2num(deblank(get(handles.editDepth1,'String')));
R.DEPTH2=str2num(deblank(get(handles.editDepth2,'String')));
R.DATEINCR=str2num(deblank(get(handles.editTimeStep,'String')));
R.LATINCR=str2num(deblank(get(handles.editLatStep,'String')));
R.LONINCR=str2num(deblank(get(handles.editLonStep,'String')));
R.LEVELINCR=str2num(deblank(get(handles.editDepthStep,'String')));
R.DIRNAME=deblank(get(handles.editDirectoryName,'String'));
R.FNPREFIX=deblank(get(handles.editFileNamePrefix,'String'));

R.SaveWorkspace='n';
if get(handles.radiobuttonInWksp,'Value') == 1
R.SaveWorkspace='y'; % save extracted data in workspace variables
end
R.SaveFiles='n';
if get(handles.radiobuttonToFiles,'Value') == 1
R.SaveFiles='y'; % save extracted data in disk files
end

S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');
R.SaveMode=deblank(S{j});

S=get(handles.popupmenuDatasetBranch,'String');j=get(handles.popupmenuDatasetBranch,'Value');
R.DataSetBranch=S{j};

% R.RequestDate=datestr(now,'yyyy-mm-dd HH:MM:SS');
request=R;
handles.request=request;

% Load Last request
function pushbuttonLoadLastRequest_Callback(hObject, eventdata, handles)
pname = strrep(mfilename('fullpath'),mfilename,'');
fname = [mfilename,'_request','.mat'];
if exist([pname,fname],'file')
load([pname,fname]);
else
msgbox(['File ',pname,fname,' .mat does not exist.']);
return
end

if ~exist('request','var')
    request = Request;
end

R=request;

%load dsi_hycom.mat % List
List=handles.DSI.List;
j = strmatch(R.DataSetBranch,List.DataSetBranch);
set(handles.popupmenuDatasetBranch,'Value',j);
CheckRanges(hObject, eventdata, handles)
DSI=LoadDSI(hObject, eventdata, handles);

% select only checked variables
for k=1:length(DSI.Fields)
FIELD=DSI.Fields{k};    
CBNUM=num2str(DSI.Fields_CBNum{k});    
eval(['set(handles.checkbox' CBNUM ',''Value'',0)'])
    for kk=1:length(R.Fields) 
        if strcmp(R.Fields{kk},FIELD) 
        eval(['set(handles.checkbox' CBNUM ',''Value'',1)']) 
        end
    end
end

% get strings from GUI menu and convert them to numbers
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
set(handles.editDepth1,'String',num2str(R.DEPTH1));
set(handles.editDepth2,'String',num2str(R.DEPTH2));
set(handles.editTimeStep,'String',num2str(R.DATEINCR));
set(handles.editLatStep,'String',num2str(R.LATINCR));
set(handles.editLonStep,'String',num2str(R.LONINCR));
set(handles.editDepthStep,'String',num2str(R.LEVELINCR));
set(handles.editDirectoryName,'String',R.DIRNAME);
set(handles.editFileNamePrefix,'String',R.FNPREFIX);

LevelsOnOff(hObject, eventdata, handles)

set(handles.radiobuttonInWksp,'Value',0);
if R.SaveWorkspace=='y'
    set(handles.radiobuttonInWksp,'Value',1); 
end
set(handles.radiobuttonToFiles,'Value',0);
if R.SaveFiles=='y'
    set(handles.radiobuttonToFiles,'Value',1); 
end

S=get(handles.popupmenuMode,'String');
j=find(strcmp(S,R.SaveMode));
if isempty(j) j=1; end
set(handles.popupmenuMode,'Value',j);

% disable saving to workspace for netcdf mode
set(handles.radiobuttonInWksp,'Enable','On');
if strcmp(R.SaveMode,'netcdf')
    set(handles.radiobuttonInWksp,'Value',0,'Enable','Off');
    set(handles.radiobuttonToFiles,'Value',1);
end

SaveFilesOnOff(hObject, eventdata, handles)

function SaveFilesOnOff(hObject, eventdata, handles)   
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end

function LevelsOnOff(hObject, eventdata, handles)   
L4D=get(handles.checkbox1,'Value')|get(handles.checkbox2,'Value')|...
    get(handles.checkbox3,'Value')|get(handles.checkbox4,'Value')|...
    get(handles.checkbox5,'Value');   
if L4D==0
    set(handles.editDepth1,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editDepth2,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editDepthStep,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif L4D==1
    set(handles.editDepth1,'BackgroundColor',[1.  1.  1. ],'Enable','On')
    set(handles.editDepth2,'BackgroundColor',[1.  1.  1. ],'Enable','On')
    set(handles.editDepthStep,'BackgroundColor',[1.  1.  1. ],'Enable','On')
end  

% --- Executes on selection change in popupmenuDatasetBranch.
function popupmenuDatasetBranch_Callback(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)
function popupmenuDatasetBranch_CreateFcn(hObject, eventdata, handles)
function CheckRanges(hObject, eventdata, handles)
% depending on the checked GUI buttons update appearence of GUI menu
load dsi_hycom.mat % List
List=handles.DSI.List;

set(handles.popupmenuDatasetBranch,'String',List.DataSetBranch(:,1)); % long name for menu

DSI=LoadDSI(hObject, eventdata, handles);

set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2,'yyyy-mm-dd') ']'])
set(handles.textLatitudeRange,'String',...
['Available latitude range: [' num2str(DSI.Range.Latitude1) ' to ' num2str(DSI.Range.Latitude2) ']'])
if ~strcmp(DSI.DataSetBranch,'HYCOM_Atlantic_Ocean_Prediction_System')
set(handles.textLatitudeRange,'String',...
['Available latitude range: [' num2str(DSI.Range.Latitude1) ' to ' num2str(DSI.Range.Latitude2) '] plus polar patch'])
end    
set(handles.textLongitudeRange,'String',...
['Available longitude range: [' num2str(DSI.Range.Longitude1) ' to ' num2str(DSI.Range.Longitude2) ']'])

%% make visible and active only checkboxes corresponding to existing variables
L=[1:1:17]; % all field checkboxes
for k=1:length(L)
CBNUM=num2str(L(k));
%eval(['set(handles.checkbox' CBNUM ',''Visible'',''off'',''Enable'',''off'')'])
eval(['set(handles.checkbox' CBNUM ',''Enable'',''off'')'])
end

for k=1:length(DSI.Fields(:,1))
    FIELDMENU=DSI.Fields_NameMenu{k};
    CBNUM=num2str(DSI.Fields_CBNum{k});
 %  eval(['set(handles.checkbox' CBNUM ',''Visible'',''on'',''Enable'',''on'',''String'',''' FIELDMENU ''')'])
    eval(['set(handles.checkbox' CBNUM ',''Enable'',''on'',''String'',''' FIELDMENU ''')'])

end
%
%S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');S=deblank(S{j});
%TXT=['Temporal resolution is ' S];
%set(handles.text39,'String',TXT,'ForegroundColor',[0 0 1])
%S=get(handles.popupmenuDatasetBranch,'String');j=get(handles.popupmenuDatasetBranch,'Value');S=deblank(S{j});
%TXT=['Spatial resolution is ' S];
%set(handles.text40,'String',TXT,'ForegroundColor',[0 0 1])

function DSI=LoadDSI(hObject, eventdata, handles)
S=get(handles.popupmenuDatasetBranch,'String');j=get(handles.popupmenuDatasetBranch,'Value');
R.DataSetBranch=S{j};

%example: load DSI_OceanColor_MODIS_4km_8day_hdf.mat
%S=['load dsi_hycom.mat DSI_' R.DataSetBranch];
%eval(['load dsi_hycom.mat DSI_' R.DataSetBranch])
eval(['DSI=handles.DSI.DSI_' R.DataSetBranch ';'])
%eval(['clear DSI_' R.DataSetBranch ';'])
