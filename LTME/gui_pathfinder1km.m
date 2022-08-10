function varargout = gui_pathfinder1km(varargin)
% Boots-up GOES GUI, interacts with user via interface.
% gui_xxx M-file for gui_xxx.fig
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
% Copyright 2007-2009
% $Version 3.0.0$

%==========================================================================
% Meri Sheremet
% 
% REVISION HISTORY:
% 2009 Template for all GUI
%==========================================================================

warning('off','MATLAB:dispatcher:InexactMatch');
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_pathfinder1km_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_pathfinder1km_OutputFcn, ...
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

%**************************************************************************
% standard functions used in all gui
% only small modifications needed: replace xxx with dataset name, enable checkbox#
%
% --- Outputs from this function are returned to the command line.
function varargout = gui_pathfinder1km_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

function checkbox1_Callback(hObject, eventdata, handles)
function checkbox2_Callback(hObject, eventdata, handles)
%function checkbox3_Callback(hObject, eventdata, handles)
% ...

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
function editTimeStep_CreateFcn(hObject, eventdata, handles)
function editTimeStep_Callback(hObject, eventdata, handles)
function editLatStep_CreateFcn(hObject, eventdata, handles)
function editLatStep_Callback(hObject, eventdata, handles)
function editLonStep_CreateFcn(hObject, eventdata, handles)
function editLonStep_Callback(hObject, eventdata, handles)
function editDirectoryName_CreateFcn(hObject, eventdata, handles)
function editDirectoryName_Callback(hObject, eventdata, handles)
function editFileNamePrefix_CreateFcn(hObject, eventdata, handles)
function editFileNamePrefix_Callback(hObject, eventdata, handles)

function textTimeRange_CreateFcn(hObject, eventdata, handles)
function textLatitudeRange_CreateFcn(hObject, eventdata, handles)
function textLongitudeRange_CreateFcn(hObject, eventdata, handles)
function popupmenuRegion_CreateFcn(hObject, eventdata, handles)
function popupmenuRegion_Callback(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)
function radiobuttonPassNight_Callback(hObject, eventdata, handles)
CheckPasses(hObject, eventdata, handles)
function radiobuttonPassDay_Callback(hObject, eventdata, handles)
CheckPasses(hObject, eventdata, handles)
function popupmenuMode_CreateFcn(hObject, eventdata, handles)
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
function radiobuttonInWksp_Callback(hObject, eventdata, handles)
function radiobuttonToFiles_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end
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
filename = 'gui_pathfinder1km_request.mat';
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
% Get data.
get_pathfinder1km(request);
% Re-enable button after obtaining data.
set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25])

function SaveFilesOnOff(hObject, eventdata, handles)   
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end

function CheckPasses(hObject, eventdata, handles)
% make sure that at least one pass is checked
% if none selected => select all
    if get(handles.radiobuttonPassNight,'Value')==0 & get(handles.radiobuttonPassDay,'Value')==0
    set(handles.radiobuttonPassNight,'Value',1);
    set(handles.radiobuttonPassDay,  'Value',1);
    end

function DSI=LoadDSI(hObject, eventdata, handles)
S=get(handles.popupmenuRegion,'String');j=get(handles.popupmenuRegion,'Value');
Region=deblank(S{j});
if strcmp(Region,'Northwest Atlantic') R.Region='NWA'; end
if strcmp(Region,'Northeast Pacific' ) R.Region='NEP'; end
if strcmp(Region,'Northwest Pacific' ) R.Region='NWP'; end

R.DataSetBranch=['Pathfinder1km_' R.Region];
%example: load DSI_OceanColor_MODIS_4km_8day_hdf.mat
%S=['load dsi_pathfinder1km.mat DSI_' R.DataSetBranch];
%eval(['load dsi_pathfinder1km.mat DSI_' R.DataSetBranch])
eval(['DSI=handles.DSI.DSI_' R.DataSetBranch ';'])
%eval(['clear DSI_' R.DataSetBranch ';'])

function CheckRanges(hObject, eventdata, handles)
% depending on the checked GUI buttons update appearence of GUI menu
DSI=LoadDSI(hObject, eventdata, handles);
% set ranges
set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2,'yyyy-mm-dd') ']'])
set(handles.textLatitudeRange,'String',...
['Available latitude range: [' num2str(DSI.Range.Latitude1) ' to ' num2str(DSI.Range.Latitude2) ']'])
set(handles.textLongitudeRange,'String',...
['Available longitude range: [' num2str(DSI.Range.Longitude1) ' to ' num2str(DSI.Range.Longitude2) ']'])
%% make visible and active only checkboxes corresponding to existing variables
L=[1:1:2]; % all field checkboxes
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some gui specific definitions

function pushbuttonGSO_Callback(hObject, eventdata, handles)
web('http://satdat1.gso.uri.edu/opendap/Pathfinder/Pathfinder1km/pathfinder_1km.html');
%web 'http://satdat1.gso.uri.edu/opendap/Pathfinder/Pathfinder1km/pathfinder_1km.html' -browser
%web 'http://satdat1.gso.uri.edu/opendap/reason/pathfinder1km/pathfinder_1km.html' -browser
function pushbuttonOpendap_Callback(hObject, eventdata, handles)
web('http://www.opendap.org/');
%web 'http://www.opendap.org/' -browser
function pushbuttonNASA_Callback(hObject, eventdata, handles)
web('http://www.nasa.gov');
%web 'http://www.nasa.gov' -browser
%web 'http://www.nasa.gov/topics/earth/index.html' -browser
function pushbuttonGSOURI_Callback(hObject, eventdata, handles)
web('http://satdat1.gso.uri.edu/opendap/Pathfinder/');
%web 'http://satdat1.gso.uri.edu/opendap/Pathfinder/' -browser

function pushbuttonInfo_Callback(hObject, eventdata, handles)
h = op_timemessage('polar',0);


%**************************************************************************


% --- Executes just before gui_pathfinder1km is made visible.
function gui_pathfinder1km_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_pathfinder1km (see VARARGIN)

DSI=load('dsi_pathfinder1km.mat');
handles.DSI=DSI;

% Define handles "openfigure" field.
if ~isfield(handles,'openfigure')
    handles.openfigure = 0; %set equal to zero on the first time.
end

% Set interface name.
filename = 'get_pathfinder1km';
txt = ['NASA/REASoN Ocean Data Portal - '];
tmp1 = ['Pathfinder1km GUI'];
tmp2 = '';
txt = [txt,tmp1];
if (~isdeployed)
    tmp2 = strrep(op_version('mfile',filename),'Version ','');
    txt = [txt,' (Ver ',tmp2,')'];
end
set(hObject,'Name',txt);

% Place logos on the GUI
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

[img,MAP]=imread('logo_gso.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;
CY=(NY+12)/ppcy;
set(handles.pushbuttonGSO,'CData',img)
pos = get(handles.pushbuttonGSO,'position'); %units will be in characters
set(handles.pushbuttonGSO,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_gso_uri.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;
CY=(NY+12)/ppcy;
set(handles.pushbuttonGSOURI,'CData',img)
pos = get(handles.pushbuttonGSOURI,'position'); %units will be in characters
set(handles.pushbuttonGSOURI,'position',[pos(1),pos(2),CX,CY]);

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
% set(handles.pushbutton10,'CData',img)
% pos = get(handles.pushbutton10,'position'); %units will be in characters
% set(handles.pushbutton10,'position',[pos(1),pos(2),CX,CY]);
[N1,N2,N3]=size(img);
set(handles.pushbuttonInfo,'string','')
set(handles.pushbuttonInfo,'CData',img)
set(handles.pushbuttonInfo,'Units','pixels')
pos = get(handles.pushbuttonInfo,'position');
set(handles.pushbuttonInfo,'position',[pos(1) pos(2) N1 N2])

CheckPasses(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)
if get(handles.checkbox1,'Value')==0 
set(handles.checkbox1,'Value',1)
end
if get(handles.radiobuttonInWksp,'Value')==0 & get(handles.radiobuttonToFiles,'Value')==0
set(handles.radiobuttonInWksp,'Value',1)
end
SaveFilesOnOff(hObject, eventdata, handles) 
% Make 'Load Last Request' button inactive when file gui_xxx_request.mat is
% absent.
if ~exist('gui_pathfinder1km_request.mat','file') 
    set(handles.pushbuttonLoadLastRequest,'Enable','Off')
end
set(handles.pushbuttonGetData,'String','Get Data','BackgroundColor',[0 1. 0.25])
%set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25]);
% Update handles "openfigure" field.
handles.openfigure = 1;
% Choose default command line output for gui_pathfinder1km
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes gui_pathfinder1km wait for user response (see UIRESUME)
% uiwait(handles.figure1);


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
R.Coordinates={'latitude','longitude','time'};
R.Passes={};
if get(handles.radiobuttonPassNight,'Value')==1 R.Passes{end+1}='night'; end
if get(handles.radiobuttonPassDay,  'Value')==1 R.Passes{end+1}='day';   end
    
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
%R.DEPTH1=str2num(deblank(get(handles.editDepth1,'String')));
%R.DEPTH2=str2num(deblank(get(handles.editDepth2,'String')));
R.DATEINCR=str2num(deblank(get(handles.editTimeStep,'String')));
R.LATINCR=str2num(deblank(get(handles.editLatStep,'String')));
R.LONINCR=str2num(deblank(get(handles.editLonStep,'String')));
%R.LEVELINCR=str2num(deblank(get(handles.editDepthStep,'String')));
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
R.DataSetBranch=DSI.DataSetBranch;
% R.RequestDate=datestr(now,'yyyy-mm-dd HH:MM:SS');
request=R;
handles.request=request;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%load dsi_pathfinder1km.mat List% List
List=handles.DSI.List;
j = strmatch(R.DataSetBranch,List.DataSetBranch);
set(handles.popupmenuRegion,'Value',j);
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

if ~isempty(strmatch('night',R.Passes)) set(handles.radiobuttonPassNight,'Value',1); end
if ~isempty(strmatch('day',  R.Passes)) set(handles.radiobuttonPassDay,'Value',1);   end

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
%set(handles.editDepth1,'String',num2str(R.DEPTH1));
%set(handles.editDepth2,'String',num2str(R.DEPTH2));
set(handles.editTimeStep,'String',num2str(R.DATEINCR));
set(handles.editLatStep,'String',num2str(R.LATINCR));
set(handles.editLonStep,'String',num2str(R.LONINCR));
%set(handles.editDepthStep,'String',num2str(R.LEVELINCR));
set(handles.editDirectoryName,'String',R.DIRNAME);
set(handles.editFileNamePrefix,'String',R.FNPREFIX);

%LevelsOnOff(hObject, eventdata, handles)

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
