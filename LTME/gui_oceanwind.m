function varargout = gui_oceanwind(varargin)
% GUI_OCEANWIND M-file for gui_oceanwind.fig
%      GUI_OCEANWIND, by itself, creates a new GUI_OCEANWIND or raises the existing
%      singleton*.
%
%      H = GUI_OCEANWIND returns the handle to a new GUI_OCEANWIND or the handle to
%      the existing singleton*.
%
%      GUI_OCEANWIND('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_OCEANWIND.M with the given input arguments.
%
%      GUI_OCEANWIND('Property','Value',...) creates a new GUI_OCEANWIND or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_oceanwind_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_oceanwind_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_oceanwind

% Last Modified by GUIDE v2.5 06-Nov-2009 14:20:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_oceanwind_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_oceanwind_OutputFcn, ...
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


% --- Executes just before gui_oceanwind is made visible.
function gui_oceanwind_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_oceanwind (see VARARGIN)

DSI=load('dsi_oceanwind.mat');
handles.DSI=DSI;

% Define handles "openfigure" field.
if ~isfield(handles,'openfigure')
    handles.openfigure = 0; %set equal to zero on the first time.
end

% Set interface name.
filename = 'get_oceanwind';
txt = [''];
tmp1 = ['OceanWind GUI'];
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

[img,MAP]=imread('logo_nasa1.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonNASA,'CData',img)
pos = get(handles.pushbuttonNASA,'position'); %units will be in characters
set(handles.pushbuttonNASA,'position',[pos(1),pos(2),CX,CY]);  

[img,MAP]=imread('logo_nasa_jpl.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;
CY=(NY+12)/ppcy;
set(handles.pushbuttonNASAJPL,'CData',img)
pos = get(handles.pushbuttonNASAJPL,'position'); %units will be in characters
set(handles.pushbuttonNASAJPL,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_podaac.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;
CY=(NY+12)/ppcy;
set(handles.pushbuttonPODAAC,'CData',img)
pos = get(handles.pushbuttonPODAAC,'position'); %units will be in characters
set(handles.pushbuttonPODAAC,'position',[pos(1),pos(2),CX,CY]);

txt=['Help'];
op_tooltipwrap(handles.pushbuttonInfo,txt,50);

[img,MAP]=imread('logo_questionmark.jpg');
% [NY,NX,N3]=size(img);
% CX=(NX+8)/ppcx;
% CY=(NY+12)/ppcy;
% set(handles.pushbutton10,'CData',img)
% pos = get(handles.pushbutton10,'position'); %units will be in characters
% set(handles.pushbutton10,'position',[pos(1),pos(2),CX,CY]);
[N1,N2,N3]=size(img);
set(handles.pushbuttonInfo,'string','')
set(handles.pushbuttonInfo,'CData',img)
set(handles.pushbuttonInfo,'Units','pixels')
pos = get(handles.pushbuttonInfo,'position');
set(handles.pushbuttonInfo,'position',[pos(1) pos(2) N1 N2])

% Set ToolTipString.
%txt=['The user can choose between QuikSCAT or SeaWinds data. ',...
%     'Currently, only QuikSCAT data are available. In addition, ',...
%     'wind stress products are unavailable.'];
%op_tooltipwrap(handles.radiobuttonQuikSCAT,txt,50);

set(handles.radiobuttonQuikSCAT,'Value',0)
set(handles.radiobuttonSeaWinds,'Value',0)
set(handles.radiobuttonCCMP,'Value',0)

% disable wind stress fields
L=200+[1:1:8]; % all field checkboxes
for k=1:length(L)
CBNUM=num2str(L(k));
%eval(['set(handles.checkbox' CBNUM ',''Visible'',''off'',''Enable'',''off'')'])
eval(['set(handles.checkbox' CBNUM ',''Enable'',''off'')'])
end
set(handles.checkboxDailyTimes1,'Enable','off','Visible','off')
set(handles.checkboxDailyTimes2,'Enable','off','Visible','off')
set(handles.checkboxDailyTimes3,'Enable','off','Visible','off')
set(handles.checkboxDailyTimes4,'Enable','off','Visible','off')
set(handles.textDailyTimes,'String','')


set(handles.editTimeStep,'String',num2str(1));
set(handles.editLatStep,'String',num2str(1));
set(handles.editLonStep,'String',num2str(1));

SaveFilesOnOff(hObject, eventdata, handles)
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

% Make 'Load Last Request' button inactive when file gui_oceanwind_request.mat is
% absent
pname = strrep(mfilename('fullpath'),mfilename,'');
fname = [mfilename,'_request','.mat'];
%'gui_oceanwind_request.mat';
if ~exist([pname,fname],'file') 
    set (handles.pushbuttonLoadLastRequest,'Enable','Off')
end

set(handles.radiobuttonCCMP,'Value',1)
SaveFilesOnOff(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)

% Update handles "openfigure" field.
handles.openfigure = 1;

% Choose default command line output for gui_oceanwind
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_oceanwind wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_oceanwind_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonOpendap.
function pushbuttonOpendap_Callback(hObject, eventdata, handles)
web 'http://www.opendap.org/' -browser
function pushbuttonNASA_Callback(hObject, eventdata, handles)
web 'http://www.nasa.gov/' -browser

% DataSetBranch selection via popupmenu is disabled
function popupmenu1_Callback(hObject, eventdata, handles)
function popupmenu1_CreateFcn(hObject, eventdata, handles)

function radiobuttonQuikSCAT_Callback(hObject, eventdata, handles)
CheckExclusive(hObject, eventdata, handles, {'radiobuttonQuikSCAT'}, {'radiobuttonSeaWinds','radiobuttonCCMP'})
CheckRanges(hObject, eventdata, handles)
function radiobuttonSeaWinds_Callback(hObject, eventdata, handles)
CheckExclusive(hObject, eventdata, handles, {'radiobuttonSeaWinds'}, {'radiobuttonQuikSCAT','radiobuttonCCMP'})
CheckRanges(hObject, eventdata, handles)
function radiobuttonCCMP_Callback(hObject, eventdata, handles)
CheckExclusive(hObject, eventdata, handles, {'radiobuttonCCMP'}, {'radiobuttonQuikSCAT','radiobuttonSeaWinds'})
CheckRanges(hObject, eventdata, handles)

function checkbox101_Callback(hObject, eventdata, handles)
function checkbox102_Callback(hObject, eventdata, handles)
function checkbox103_Callback(hObject, eventdata, handles)
function checkbox104_Callback(hObject, eventdata, handles)
function checkbox105_Callback(hObject, eventdata, handles)
function checkbox106_Callback(hObject, eventdata, handles)
function checkbox107_Callback(hObject, eventdata, handles)
function checkbox108_Callback(hObject, eventdata, handles)
function checkbox109_Callback(hObject, eventdata, handles)
function checkbox110_Callback(hObject, eventdata, handles)
function checkbox111_Callback(hObject, eventdata, handles)
function checkbox112_Callback(hObject, eventdata, handles)
function checkbox113_Callback(hObject, eventdata, handles)
function checkbox201_Callback(hObject, eventdata, handles)
function checkbox202_Callback(hObject, eventdata, handles)
function checkbox203_Callback(hObject, eventdata, handles)
function checkbox204_Callback(hObject, eventdata, handles)
function checkbox205_Callback(hObject, eventdata, handles)
function checkbox206_Callback(hObject, eventdata, handles)
function checkbox207_Callback(hObject, eventdata, handles)
function checkbox208_Callback(hObject, eventdata, handles)
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

function checkboxDailyTimes1_Callback(hObject, eventdata, handles)
function checkboxDailyTimes2_Callback(hObject, eventdata, handles)
function checkboxDailyTimes3_Callback(hObject, eventdata, handles)
function checkboxDailyTimes4_Callback(hObject, eventdata, handles)

function pushbuttonInfo_Callback(hObject, eventdata, handles)

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

function pushbuttonBrowse_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonToFiles,'Value') == 1
DIRNAME=uigetdir(pwd,'Select a Directory for Downloading Files');
set(handles.editDirectoryName,'String',DIRNAME)
end

% The main button that starts data acquisition
function pushbuttonGetData_Callback(hObject, eventdata, handles)
% set(handles.pushbuttonGetData,'String','For a new request reissue get_oceanwind')
% uiresume(handles.figure1)
%% Turn off ability to re-issue "Get Data"
set(hObject,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off','String','Please wait ...');
pause(.1)

%% GET INPUT FROM GUI 
request=InputGet(hObject, eventdata, handles);

%save gui_oceanwind_request.mat request -MAT %-V6 %-APPEND
% Remember request for future sessions.
pname = strrep(mfilename('fullpath'),mfilename,'');
fname = [mfilename,'_request','.mat'];
save([pname,fname],'request','-mat')

%'gui_oceanwind_request.mat';
if exist([pname,fname],'file') 
    set (handles.pushbuttonLoadLastRequest,'Enable','On')
end

% If it passes inspection, send down to user's workspace
% as variable "request". Then the user can get the data
%assignin('base','request',request);

get_oceanwind(request);

set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25]);

function pushbuttonNASAJPL_Callback(hObject, eventdata, handles)
function pushbuttonPODAAC_Callback(hObject, eventdata, handles)

function request=InputGet(hObject, eventdata, handles)
DSI=LoadDSI(hObject, eventdata, handles);
R.Fields={};
for k=1:length(DSI.Fields_CBNum)
FIELD=DSI.Fields{k};    
CBNUM=num2str(DSI.Fields_CBNum{k});
    if eval(['get(handles.checkbox' CBNUM ',''Value'')']) == 1
    R.Fields{end+1}=FIELD;    
    end
end
R.Coordinates={'latitude','longitude','time'};

R.DailyTimes=[];
if strcmp(get(handles.checkboxDailyTimes1,'Enable'),'on') & get(handles.checkboxDailyTimes1,'Value')==1 R.DailyTimes=[R.DailyTimes 1]; end
if strcmp(get(handles.checkboxDailyTimes2,'Enable'),'on') & get(handles.checkboxDailyTimes2,'Value')==1 R.DailyTimes=[R.DailyTimes 2]; end
if strcmp(get(handles.checkboxDailyTimes3,'Enable'),'on') & get(handles.checkboxDailyTimes3,'Value')==1 R.DailyTimes=[R.DailyTimes 3]; end
if strcmp(get(handles.checkboxDailyTimes4,'Enable'),'on') & get(handles.checkboxDailyTimes4,'Value')==1 R.DailyTimes=[R.DailyTimes 4]; end

R.DataSetBranch=DSI.DataSetBranch;

R.DATE1='00000000';
S=deblank(get(handles.editTime1yyyy,'String'));R.DATE1(4-length(S)+1:4)=S;
S=deblank(get(handles.editTime1mm,'String'));R.DATE1(6-length(S)+1:6)=S;
S=deblank(get(handles.editTime1dd,'String'));R.DATE1(8-length(S)+1:8)=S;
R.DATE2='00000000';
S=deblank(get(handles.editTime2yyyy,'String'));R.DATE2(4-length(S)+1:4)=S;
S=deblank(get(handles.editTime2mm,'String'));R.DATE2(6-length(S)+1:6)=S;
S=deblank(get(handles.editTime2dd,'String'));R.DATE2(8-length(S)+1:8)=S;
R.DATEINCR=str2num(deblank(get(handles.editTimeStep,'String')));

R.LAT1=str2num(deblank(get(handles.editLat1,'String')));
R.LAT2=str2num(deblank(get(handles.editLat2,'String')));
R.LON1=str2num(deblank(get(handles.editLon1,'String')));
R.LON2=str2num(deblank(get(handles.editLon2,'String')));
R.LATINCR=str2num(deblank(get(handles.editLatStep,'String')));
R.LONINCR=str2num(deblank(get(handles.editLonStep,'String')));

R.SaveWorkspace='n'; 
if get(handles.radiobuttonInWksp,'Value') == 1 R.SaveWorkspace='y'; end
S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');
R.SaveMode=deblank(S{j});
R.SaveFiles='n'; 
if get(handles.radiobuttonToFiles,'Value') == 1 R.SaveFiles='y'; end
R.DIRNAME =deblank(get(handles.editDirectoryName,'String'));
if isempty(R.DIRNAME) R.DIRNAME=''; end
R.FNPREFIX=deblank(get(handles.editFileNamePrefix,'String'));

%R.RequestDate=datestr(now,'yyyy-mm-dd HH:MM:SS');
request=R;
handles.request=request;

% Load Last request
function pushbuttonLoadLastRequest_Callback(hObject, eventdata, handles)
filename = 'gui_oceanwind_request.mat';
if exist(filename,'file')
    load(filename)
else
    msgbox(['File ',filename,' does not exist.']);
    return
end
R=request;

set(handles.radiobuttonQuikSCAT,'Value',0)
set(handles.radiobuttonSeaWinds,'Value',0)
set(handles.radiobuttonCCMP,'Value',0)
if strcmp(R.DataSetBranch,'OceanWind_QuikSCAT')
    set(handles.radiobuttonQuikSCAT,'Value',1)
end
if strcmp(R.DataSetBranch,'OceanWind_SeaWinds')
    set(handles.radiobuttonSeaWinds,'Value',1)
end
if strcmp(R.DataSetBranch,'OceanWind_CCMP')
    set(handles.radiobuttonCCMP,'Value',1)
end

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

set(handles.checkboxDailyTimes1,'Value',0)
set(handles.checkboxDailyTimes2,'Value',0)
set(handles.checkboxDailyTimes3,'Value',0)
set(handles.checkboxDailyTimes4,'Value',0)
if ~isempty(find(R.DailyTimes==1)) set(handles.checkboxDailyTimes1,'Value',1); end
if ~isempty(find(R.DailyTimes==2)) set(handles.checkboxDailyTimes2,'Value',1); end
if ~isempty(find(R.DailyTimes==3)) set(handles.checkboxDailyTimes3,'Value',1); end
if ~isempty(find(R.DailyTimes==4)) set(handles.checkboxDailyTimes4,'Value',1); end

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
set(handles.editLatStep,'String',num2str(R.LATINCR));
set(handles.editLonStep,'String',num2str(R.LONINCR));
set(handles.editDirectoryName,'String',R.DIRNAME);
set(handles.editFileNamePrefix,'String',R.FNPREFIX);
set(handles.editTimeStep,'String',num2str(R.DATEINCR));

set(handles.radiobuttonInWksp,'Value',0);
if R.SaveWorkspace=='y'
    set(handles.radiobuttonInWksp,'Value',1); 
end
S=get(handles.popupmenuMode,'String');
j=find(strcmp(S,R.SaveMode)); if isempty(j) j=1; end
set(handles.popupmenuMode,'Value',j);
set(handles.radiobuttonToFiles,'Value',0);
if R.SaveFiles=='y'
    set(handles.radiobuttonToFiles,'Value',1); 
end

S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');
R.SaveMode=S{j};
% disable saving to workspace for netcdf mode
set(handles.radiobuttonInWksp,'Enable','On');
if strcmp(R.SaveMode,'netcdf')
    set(handles.radiobuttonInWksp,'Value',0,'Enable','Off');
    set(handles.radiobuttonToFiles,'Value',1);
end
SaveFilesOnOff(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)


function DSI=LoadDSI(hObject, eventdata, handles)
if strcmp(get(handles.radiobuttonQuikSCAT,'Enable'),'on') & get(handles.radiobuttonQuikSCAT,'Value')==1
    DSI=handles.DSI.DSI_OceanWind_QuikSCAT;
end
if strcmp(get(handles.radiobuttonSeaWinds,'Enable'),'on') & get(handles.radiobuttonSeaWinds,'Value')==1
    DSI=handles.DSI.DSI_OceanWind_SeaWinds;
end
if strcmp(get(handles.radiobuttonCCMP,'Enable'),'on') & get(handles.radiobuttonCCMP,'Value')==1
    DSI=handles.DSI.DSI_OceanWind_CCMP;
end

function SaveFilesOnOff(hObject, eventdata, handles)
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end
function CheckExclusive(hObject, eventdata, handles, Us, Them)
% sets value of Us to 1 and Them to 0
% example: Us={'radiobuttonInWksp'}; Them={'radiobuttonToFiles';'radiobutton3'}
for k=1:length(Us)
eval(['set(handles.' Us{k} ',''Value'',1)'])
end
for k=1:length(Them)
eval(['set(handles.' Them{k} ',''Value'',0)'])
end

function CheckRanges(hObject, eventdata, handles)
DSI=LoadDSI(hObject, eventdata, handles);
% make visible and active only checkboxes corresponding to existing variables
%% make visible and active only checkboxes corresponding to existing variables
L=100+[1:1:13]; % all field checkboxes
for k=1:length(L)
CBNUM=num2str(L(k));
eval(['set(handles.checkbox' CBNUM ',''Visible'',''off'',''Enable'',''off'')'])
%eval(['set(handles.checkbox' CBNUM ',''Enable'',''off'')'])
end
for k=1:length(DSI.Fields(:,1))
    FIELDMENU=DSI.Fields_NameMenu{k};
    CBNUM=num2str(DSI.Fields_CBNum{k});
    eval(['set(handles.checkbox' CBNUM ',''Visible'',''on'',''Enable'',''on'',''String'',''' FIELDMENU ''')'])
 %   eval(['set(handles.checkbox' CBNUM ',''Enable'',''on'',''String'',''' FIELDMENU ''')'])
    S=DSI.Fields_ToolTip{k}; 
    eval(['op_tooltipwrap(handles.checkbox' CBNUM ',S,50)'])
end

set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2,'yyyy-mm-dd') ']'],'ForegroundColor',[0 0 1])
set(handles.textLatitudeRange,'String',...
['Available latitude range: [' num2str(DSI.Range.Latitude1) ' to ' num2str(DSI.Range.Latitude2) ']'],'ForegroundColor',[0 0 1])
set(handles.textLongitudeRange,'String',...
['Available longitude range: [' num2str(DSI.Range.Longitude1) ' to ' num2str(DSI.Range.Longitude2) ']'],'ForegroundColor',[0 0 1])
set(handles.textTemporalResolution,'String', ['Temporal resolution is ' DSI.Range.TRes],'ForegroundColor',[0 0 1])
set(handles.textSpatialResolution ,'String', ['Spatial resolution is ' DSI.Range.HRes],'ForegroundColor',[0 0 1])

set(handles.checkboxDailyTimes1,'Enable','off','Visible','off')
set(handles.checkboxDailyTimes2,'Enable','off','Visible','off')
set(handles.checkboxDailyTimes3,'Enable','off','Visible','off')
set(handles.checkboxDailyTimes4,'Enable','off','Visible','off')
if strcmp(DSI.DataSetBranch,'OceanWind_QuikSCAT')
set(handles.textDailyTimes,'String','Pass:')
set(handles.checkboxDailyTimes1,'Enable','on','Visible','on','String','ascending:  ~6am LST')
set(handles.checkboxDailyTimes2,'Enable','on','Visible','on','String','descending: ~6pm LST')
end    
if strcmp(DSI.DataSetBranch,'OceanWind_SeaWinds')
set(handles.textDailyTimes,'String','Pass:')
set(handles.checkboxDailyTimes1,'Enable','on','Visible','on','String','ascending:  ~3am LST')
set(handles.checkboxDailyTimes2,'Enable','on','Visible','on','String','descending: ~3pm LST')
end    
if strcmp(DSI.DataSetBranch,'OceanWind_CCMP')
set(handles.textDailyTimes,'String','Daily Times:')
set(handles.checkboxDailyTimes1,'Enable','on','Visible','on','String','00h')
set(handles.checkboxDailyTimes2,'Enable','on','Visible','on','String','06h')
set(handles.checkboxDailyTimes3,'Enable','on','Visible','on','String','12h')
set(handles.checkboxDailyTimes4,'Enable','on','Visible','on','String','18h')
end    

