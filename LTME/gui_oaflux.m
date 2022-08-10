function varargout = gui_oaflux(varargin)
% Boots-up OAFlux GUI, interacts with user via interface.
% gui_oaflux M-file for gui_oaflux.fig
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
% 2007/06/01 created, ms
% 2008/06/15 modified for compiled version, ceb
%==========================================================================

warning('off','MATLAB:dispatcher:InexactMatch');
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_oaflux_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_oaflux_OutputFcn, ...
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

% --- Executes just before gui_oaflux is made visible.
function gui_oaflux_OpeningFcn(hObject, eventdata, handles, varargin)

DSI=load('dsi_oaflux.mat');
handles.DSI=DSI;

% Define handles "openfigure" field.
if ~isfield(handles,'openfigure')
    handles.openfigure = 0; %set equal to zero on the first time.
end

% Set interface name.
filename = 'get_oaflux';
txt = ['NASA/REASoN Ocean Data Portal - '];
tmp1 = ['OAFlux GUI'];
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

[img,MAP]=imread('logo_whoi.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonWHOI,'CData',img)
pos = get(handles.pushbuttonWHOI,'position'); %units will be in characters
set(handles.pushbuttonWHOI,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_apdrc.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonapdrc,'CData',img)
pos = get(handles.pushbuttonapdrc,'position'); %units will be in characters
set(handles.pushbuttonapdrc,'position',[pos(1),pos(2),CX,CY]);

set(handles.editTimeStep,'String',num2str(1));
set(handles.editLatStep,'String',num2str(1));
set(handles.editLonStep,'String',num2str(1));
%set(handles.editFileNamePrefix,'String','opendap_');

if get(handles.radiobuttonInWksp,'Value')==0 & get(handles.radiobuttonToFiles,'Value')==0
    set(handles.radiobuttonInWksp,'Value',1)
end

% Make 'Load Last Request' button inactive when file gui_xxx_request.mat is
% absent.
if ~exist('gui_oaflux_request.mat','file') 
    set(handles.pushbuttonLoadLastRequest,'Enable','Off')
end

% set monthly as default
set(handles.popupmenuTemporalAverages,'Value',2);


SaveFilesOnOff(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)

%set(handles.pushbuttonGetData,'String','Get Data','BackgroundColor',[0 1. 0.25])
% Update handles "openfigure" field.
handles.openfigure = 1;
% Choose default command line output for gui_oaflux
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_oaflux wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_oaflux_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;

function checkbox1_Callback(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)
function checkbox2_Callback(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)
function checkbox3_Callback(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)
function checkbox4_Callback(hObject, eventdata, handles)
function checkbox5_Callback(hObject, eventdata, handles)
function checkbox6_Callback(hObject, eventdata, handles)
function checkbox7_Callback(hObject, eventdata, handles)
function checkbox8_Callback(hObject, eventdata, handles)
function checkbox9_Callback(hObject, eventdata, handles)
function checkbox10_Callback(hObject, eventdata, handles)

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
function editLatStep_Callback(hObject, eventdata, handles)
function editLatStep_CreateFcn(hObject, eventdata, handles)
function editLonStep_Callback(hObject, eventdata, handles)
function editLonStep_CreateFcn(hObject, eventdata, handles)
function editDirectoryName_Callback(hObject, eventdata, handles)
function editDirectoryName_CreateFcn(hObject, eventdata, handles)
function editFileNamePrefix_Callback(hObject, eventdata, handles)
function editFileNamePrefix_CreateFcn(hObject, eventdata, handles)
function editTimeStep_Callback(hObject, eventdata, handles)
function editTimeStep_CreateFcn(hObject, eventdata, handles)

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

function popupmenuTemporalAverages_Callback(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)
function popupmenuTemporalAverages_CreateFcn(hObject, eventdata, handles)

function radiobuttonInWksp_Callback(hObject, eventdata, handles)
function radiobuttonToFiles_Callback(hObject, eventdata, handles)
SaveFilesOnOff(hObject, eventdata, handles)

function pushbuttonBrowse_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonToFiles,'Value') == 1
DIRNAME=uigetdir(pwd,'Select a Directory for Downloading Files');
set(handles.editDirectoryName,'String',DIRNAME)
end

% The main button that starts data acquisition
function pushbuttonGetData_Callback(hObject, eventdata, handles)
% uiresume(handles.figure1)
% Turn off ability to re-issue "Get Data."
set(hObject,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off','String','Please wait...');
pause(.1)

% Get input from user.
request=InputGet(hObject, eventdata, handles);

% Save request to .mat file.
filename = 'gui_oaflux_request.mat';
%if (~isdeployed)
    pname = mfilename('fullpath');
    pname = pname(1:length(pname)-length(mfilename)); %mfilename does not have .m in it
    filename = [pname,filename];
%end
save(filename,'request','-mat')
set(handles.pushbuttonLoadLastRequest,'Enable','On')

% If it passes inspection, send down to user's workspace
% as variable "request". Then the user can get the data
% assignin('base','request',request);

get_oaflux(request);

set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25]);

function pushbuttonWHOI_Callback(hObject, eventdata, handles)
web('http://oaflux.whoi.edu/');
%web('http://apdrc.soest.hawaii.edu/datadoc/whoi_oaswlw.htm','-browser');
function pushbuttonapdrc_Callback(hObject, eventdata, handles)
web('http://apdrc.soest.hawaii.edu/w_data/alldata3.htm');
%web 'http://apdrc.soest.hawaii.edu/dods/public_data/WHOI_OAFlux' -browser  
function pushbuttonOpendap_Callback(hObject, eventdata, handles)
web('http://www.opendap.org/');
function pushbuttonNASA_Callback(hObject, eventdata, handles)
%web 'http://www.nasa.gov/topics/earth/index.html' -browser
web('http://www.nasa.gov');

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

S=get(handles.popupmenuTemporalAverages,'String');j=get(handles.popupmenuTemporalAverages,'Value');
R.TAVG=deblank(S{j});
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

% make sure that both dataset requests have same time coverage 
%if (R.Fields{1,2}=='y' | R.Fields{2,2}=='y') & str2num(R.DATE1)<19830701
%    R.DATE1='19830701';
%end

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
filename = 'gui_oaflux_request.mat';
if exist(filename,'file')
    load(filename)
else
    msgbox(['File ',filename,' does not exist.']);
    return
end

R=request;

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
set(handles.editLatStep,'String',num2str(R.LATINCR));
set(handles.editLonStep,'String',num2str(R.LONINCR));
set(handles.editDirectoryName,'String',R.DIRNAME);
set(handles.editFileNamePrefix,'String',R.FNPREFIX);
set(handles.editTimeStep,'String',num2str(R.DATEINCR));

S=get(handles.popupmenuTemporalAverages,'String');
j=find(strcmp(S,R.TAVG)); if isempty(j) j=1; end
set(handles.popupmenuTemporalAverages,'Value',j);

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

function SaveFilesOnOff(hObject, eventdata, handles)   
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end


function DSI=LoadDSI(hObject, eventdata, handles)
%load dsi_oaflux.mat
DataSetBranchList=handles.DSI.DataSetBranchList;
%S=get(handles.popupmenuTemporalAverages,'String');
j=get(handles.popupmenuTemporalAverages,'Value');
eval(['DSI=handles.DSI.DSI_' DataSetBranchList{j} ';'])
%DSI.DataSetBranchList=DataSetBranchList;

function CheckRanges(hObject, eventdata, handles)
DSI=LoadDSI(hObject, eventdata, handles);
% make visible and active only checkboxes corresponding to existing variables
%% make visible and active only checkboxes corresponding to existing variables
L=[1:1:10]; % all field checkboxes
for k=1:length(L)
CBNUM=num2str(L(k));
eval(['set(handles.checkbox' CBNUM ',''Visible'',''on'',''Enable'',''off'')'])
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

% fix for daily fields nlwrs,nswrs,qnet
if strcmp(DSI.DataSetBranch,'OAFlux_daily')
    if get(handles.checkbox1,'Value')==1 
set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2_nlwrs,'yyyy-mm-dd') ']  for nlwrs'],'ForegroundColor',[1 0 0])
    end
    if get(handles.checkbox2,'Value')==1 
set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2_nswrs,'yyyy-mm-dd') ']  for nswrs'],'ForegroundColor',[1 0 0])
    end
    if get(handles.checkbox3,'Value')==1 
set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2_qnet,'yyyy-mm-dd') ']  for qnet'],'ForegroundColor',[1 0 0])
    end       
end
