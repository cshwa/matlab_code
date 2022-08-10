function varargout = gui_ghrsst(varargin)
% Boots-up GHRSST GUI, interacts with user via interface.
% gui_ghrsst M-file for gui_ghrsst.fig
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
% Copyright 2008
% $Version 3.0.0$

%==========================================================================
% Meri Sheremet
% June 2008
%
% REVISION HISTORY:
% 2008/06/01 created, ms
% 2008/06/15 modified for compiled version, ceb
%==========================================================================

% Last Modified by GUIDE v2.5 04-Feb-2009 05:48:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_ghrsst_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_ghrsst_OutputFcn, ...
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


% --- Executes just before gui_ghrsst is made visible.
function gui_ghrsst_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_ghrsst (see VARARGIN)

DSI=load('dsi_ghrsst.mat');
handles.DSI=DSI;

% Define handles "openfigure" field.
if ~isfield(handles,'openfigure')
    handles.openfigure = 0; %set equal to zero on the first time.
end

% Set interface name.
filename = 'get_ghrsst';
txt = ['NASA/REASoN Ocean Data Portal - '];
tmp1 = ['GHRSST GUI'];
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

[img,MAP]=imread('logo_ghrsst.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonGHRSST,'CData',img)
pos = get(handles.pushbuttonGHRSST,'position'); %units will be in characters
set(handles.pushbuttonGHRSST,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_gdac.jpg');
[NY,NX,N3]=size(img);
CX=(NX+6)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+6)/ppcy;
set(handles.pushbuttonGHRSSTjplnasa,'CData',img)
pos = get(handles.pushbuttonGHRSSTjplnasa,'position'); %units will be in characters
set(handles.pushbuttonGHRSSTjplnasa,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_ltsrf.jpg');
[NY,NX,N3]=size(img);
CX=(NX+6)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+6)/ppcy;
set(handles.pushbuttonGHRSSTnodcnoaa,'CData',img)
pos = get(handles.pushbuttonGHRSSTnodcnoaa,'position'); %units will be in characters
set(handles.pushbuttonGHRSSTnodcnoaa,'position',[pos(1),pos(2),CX,CY]);

set(handles.editTimeStep,'String',num2str(1));
set(handles.editLatStep,'String',num2str(1));
set(handles.editLonStep,'String',num2str(1));
%set(handles.editFileNamePrefix,'String','opendap_');

% Make 'Load Last Request' button inactive when file gui_xxx_request.mat is
% absent.
if ~exist('gui_ghrsst_request.mat','file') 
    set(handles.pushbuttonLoadLastRequest,'Enable','Off')
end

SaveFilesOnOff(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)

% Update handles "openfigure" field.
handles.openfigure = 1;

% Choose default command line output for gui_ghrsst
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_ghrsst wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_ghrsst_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

function popupmenuLevel_Callback(hObject, eventdata, handles)   %Level
CheckRanges(hObject, eventdata, handles)
function popupmenuLevel_CreateFcn(hObject, eventdata, handles)
function popupmenuProduct_Callback(hObject, eventdata, handles)   %Product
CheckRanges(hObject, eventdata, handles)
function popupmenuProduct_CreateFcn(hObject, eventdata, handles)
function popupmenuRegion_Callback(hObject, eventdata, handles)   %Region
CheckRanges(hObject, eventdata, handles)
function popupmenuRegion_CreateFcn(hObject, eventdata, handles)

function popupmenuMode_Callback(hObject, eventdata, handles)   %Mode
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

function popupmenuDatasetBranch_Callback(hObject, eventdata, handles)   %Dataset Branch
function popupmenuDatasetBranch_CreateFcn(hObject, eventdata, handles)
function popupmenu6_Callback(hObject, eventdata, handles)   %Temporal averages
function popupmenu6_CreateFcn(hObject, eventdata, handles)
function popupmenu7_Callback(hObject, eventdata, handles)   %Spatial resolution
function popupmenu7_CreateFcn(hObject, eventdata, handles)

function checkbox1_Callback(hObject, eventdata, handles)
function checkbox2_Callback(hObject, eventdata, handles)
function checkbox3_Callback(hObject, eventdata, handles)
function checkbox4_Callback(hObject, eventdata, handles)
function checkbox5_Callback(hObject, eventdata, handles)
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
function checkbox18_Callback(hObject, eventdata, handles)
function checkbox19_Callback(hObject, eventdata, handles)

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
filename = 'gui_ghrsst_request.mat';
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

get_ghrsst(request);

set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25]);

function request=InputGet(hObject, eventdata, handles)

S=get(handles.popupmenuDatasetBranch,'String');j=get(handles.popupmenuDatasetBranch,'Value');
if ~iscell(S) S={S}; end
DataSetBranch=S{j};
R.DataSetBranch=DataSetBranch;

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

%S=get(handles.popupmenu6,'String');j=get(handles.popupmenu6,'Value');
%if ~iscell(S) S={S}; end
%R.TAvg=S{j};
%S=get(handles.popupmenu7,'String');j=get(handles.popupmenu7,'Value');
%if ~iscell(S) S={S}; end
%R.HRes=S{j};

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

function pushbuttonGHRSST_Callback(hObject, eventdata, handles)
web('http://www.ghrsst-pp.org/');
%web 'http://www.ghrsst-pp.org/' -browser
function pushbuttonOpendap_Callback(hObject, eventdata, handles)
web('http://www.opendap.org/');
%web 'http://www.opendap.org/' -browser
function pushbuttonNASA_Callback(hObject, eventdata, handles)
web('http://www.nasa.gov');
%web 'http://www.nasa.gov' -browser
function pushbuttonGHRSSTjplnasa_Callback(hObject, eventdata, handles)
web('http://ghrsst.jpl.nasa.gov/');
%web 'http://ghrsst.jpl.nasa.gov/' -browser
function pushbuttonGHRSSTnodcnoaa_Callback(hObject, eventdata, handles)
web('http://ghrsst.nodc.noaa.gov/');
%web 'http://ghrsst.nodc.noaa.gov/' browser
   
% Load Last request
function pushbuttonLoadLastRequest_Callback(hObject, eventdata, handles)
filename = 'gui_ghrsst_request.mat';
if exist(filename,'file')
    load(filename)
else
    msgbox(['File ',filename,' does not exist.']);
    return
end

R=request;

%load dsi_ghrsst.mat
List=handles.DSI.List;

DataSetBranch=R.DataSetBranch;
set(handles.popupmenuDatasetBranch,'String',List.DataSetBranch(:,2))
j = strmatch(DataSetBranch,List.DataSetBranch(:,1),'exact');
set(handles.popupmenuDatasetBranch,'Value',j)
% Selected DataSetBranch
eval(['DSI=handles.DSI.DSI_' DataSetBranch ';'])
%DSI.List=List;

CheckLoadLastReq(hObject, eventdata, handles)
%DSI=LoadDSI(hObject, eventdata, handles);

% select only checked variables
for k=1:length(DSI.Fields_CBNum)
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
set(handles.editTimeStep,'String',num2str(R.DATEINCR));
set(handles.editLatStep,'String',num2str(R.LATINCR));
set(handles.editLonStep,'String',num2str(R.LONINCR));
set(handles.editDirectoryName,'String',R.DIRNAME);
set(handles.editFileNamePrefix,'String',R.FNPREFIX);

%S=get(handles.popupmenuProduct,'String');
%j=find(strcmp(S,R.TAVG)); if isempty(j) j=1; end
%set(handles.popupmenuProduct,'Value',j);

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

% disable saving to workspace for netcdf mode
set(handles.radiobuttonInWksp,'Enable','On');
if strcmp(R.SaveMode,'netcdf')
    set(handles.radiobuttonInWksp,'Value',0,'Enable','Off');
    set(handles.radiobuttonToFiles,'Value',1);
end

SaveFilesOnOff(hObject, eventdata, handles)
CheckRanges(hObject, eventdata, handles)


function SaveFilesOnOff(hObject, eventdata, handles)
if get(handles.radiobuttonInWksp,'Value')==0 & get(handles.radiobuttonToFiles,'Value')==0
    set(handles.radiobuttonInWksp,'Value',1)
end

if isempty(get(handles.popupmenuProduct,'Value'))
    set(handles.popupmenuProduct,'Value',1);
end

if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end

function DSI=LoadDSI(hObject, eventdata, handles)
%load dsi_ghrsst.mat
S=get(handles.popupmenuDatasetBranch,'String');j=get(handles.popupmenuDatasetBranch,'Value');
if ~iscell(S) S={S}; end
DataSetBranch=S{j};
eval(['DSI=handles.DSI.DSI_' DataSetBranch ';'])
%DSI.List=List;

function CheckRanges(hObject, eventdata, handles)
% depending on the checked GUI buttons update appearence of GUI menu
%load dsi_ghrsst.mat % List
List=handles.DSI.List;

LevelL=List.Level;
set(handles.popupmenuLevel,'String',LevelL(:,2)); % long name
j=get(handles.popupmenuLevel,'Value');
if isempty(j) | (j>length(LevelL(:,1)))
    j=1; set(handles.popupmenuLevel,'Value');
end

Level=LevelL(j,1); % short name

%Product
[N1,N2]=size(List.DataSetBranch);
jj=[];
for k=1:N1
        if strcmp(List.DataSetBranch{k,3},Level)
        Product=List.DataSetBranch{k,4};
        j = strmatch(Product,List.Product(:,1),'exact');
        if ~isempty(j)
            jj=[jj ;j];
        else
            disp('warning: Product is not on the list')
            Product
            List.Product(:,1)
        end
    end
end
jj=unique(sort(jj)); % sort and remove duplicate entries
ProductL=List.Product(jj,:);
set(handles.popupmenuProduct,'String',ProductL(:,2)); % long name
j=get(handles.popupmenuProduct,'Value');
if isempty(j) | (j>length(ProductL(:,1)))
    j=1; set(handles.popupmenuProduct,'Value',j);
end
Product=ProductL(j,1); % short name

%Region
[N1,N2]=size(List.DataSetBranch);
jj=[];
for k=1:N1
    if strcmp(List.DataSetBranch{k,3},Level) & strcmp(List.DataSetBranch{k,4},Product)
        Region=List.DataSetBranch{k,5};
        j = strmatch(Region,List.Region(:,1),'exact');
        if ~isempty(j)
            jj=[jj ;j];
        else
            disp('warning: Region is not on the list')
            Region
            List.Region(:,1)
        end
    end
end
jj=unique(sort(jj)); % sort and remove duplicate entries
RegionL=List.Region(jj,:);
set(handles.popupmenuRegion,'String',RegionL(:,2)); % long name
j=get(handles.popupmenuRegion,'Value');
if isempty(j) | (j>length(RegionL(:,1)))
    j=1; set(handles.popupmenuRegion,'Value',j);
end
Region=RegionL(j,1); % short name

%DataSetBranch
[N1,N2]=size(List.DataSetBranch);
jj=[];
for k=1:N1
    if strcmp(List.DataSetBranch{k,3},Level) & strcmp(List.DataSetBranch{k,4},Product) & strcmp(List.DataSetBranch{k,5},Region)
        jj=[jj ;k];
    end
end
jj=unique(sort(jj)); % sort and remove duplicate entries
if length(jj)>1
disp('warning: DataSetBranch is not unique')
end
j=jj(1);
DataSetBranch=List.DataSetBranch{j,1};
set(handles.popupmenuDatasetBranch,'String',List.DataSetBranch(:,2)); % long name
set(handles.popupmenuDatasetBranch,'Value',j);

eval(['DSI=handles.DSI.DSI_' DataSetBranch ';'])
DSI.TRes=List.DataSetBranch{j,6};
DSI.HRes=List.DataSetBranch{j,7};

%set(handles.popupmenu6,'String',DSI.TRes,'Value',1);
%set(handles.popupmenu7,'String',DSI.HRes,'Value',1);

% CheckBoxes
L=[1:1:22]; % possible checkbox numbers
for k=1:length(L) % first make all field checkboxes inactive
CBNUM=num2str(L(k));
eval(['set(handles.checkbox' CBNUM ',''Enable'',''off'')'])
end
for k=1:length(DSI.Fields_NameMenu) % make active only checked fields
    FIELDMENU=DSI.Fields_NameMenu{k};
    CBNUM=num2str(DSI.Fields_CBNum{k});
    eval(['set(handles.checkbox' CBNUM ',''Visible'',''on'',''Enable'',''on'',''String'',''' FIELDMENU ''')'])
end

%Ranges
set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2,'yyyy-mm-dd') ']'],'ForegroundColor',[0 0 1])
set(handles.textLatitudeRange,'String',...
['Available latitude range: [' num2str(DSI.Range.Latitude1) ' to ' num2str(DSI.Range.Latitude2) ']'],'ForegroundColor',[0 0 1])
set(handles.textLongitudeRange,'String',...
['Available longitude range: [' num2str(DSI.Range.Longitude1) ' to ' num2str(DSI.Range.Longitude2) ']'],'ForegroundColor',[0 0 1])
set(handles.textTemporalResolution,'String', ['Temporal resolution is ' DSI.TRes],'ForegroundColor',[0 0 1])
set(handles.textSpatialResolution ,'String', ['Spatial resolution is ' DSI.HRes],'ForegroundColor',[0 0 1])

function CheckExclusive(hObject, eventdata, handles, Us, Them)
% sets value of Us to 1 and Them to 0
% example: Us={'radiobuttonInWksp'}; Them={'radiobuttonToFiles';'radiobutton3'}
for k=1:length(Us)
eval(['set(handles.' Us{k} ',''Value'',1)'])
end
for k=1:length(Them)
eval(['set(handles.' Them{k} ',''Value'',0)'])
end

function CheckLoadLastReq(hObject, eventdata, handles)
% depending on the checked GUI buttons update appearence of GUI menu
%load dsi_ghrsst.mat % List
List=handles.DSI.List;

S=get(handles.popupmenuDatasetBranch,'String');j=get(handles.popupmenuDatasetBranch,'Value');
if ~iscell(S) S={S}; end
DataSetBranch=S{j};

jds=j;
Level=List.DataSetBranch{jds,3};

LevelL=List.Level;
set(handles.popupmenuLevel,'String',LevelL(:,2)); % long name
j = strmatch(Level,LevelL(:,1),'exact');
if isempty(j) 
    j=1; 
end
set(handles.popupmenuLevel,'Value',j);

%Product
[N1,N2]=size(List.DataSetBranch);
jj=[];
for k=1:N1
        if strcmp(List.DataSetBranch{k,3},Level)
        Product=List.DataSetBranch{k,4};
        j = strmatch(Product,List.Product(:,1),'exact');
        if ~isempty(j)
            jj=[jj ;j];
        else
            disp('warning: Product is not on the list')
            Product
            List.Product(:,1)
        end
    end
end
jj=unique(sort(jj)); % sort and remove duplicate entries
ProductL=List.Product(jj,:);
set(handles.popupmenuProduct,'String',ProductL(:,2)); % long name
Product=List.DataSetBranch{jds,4};
j = strmatch(Product,ProductL(:,1),'exact');
if isempty(j) 
    j=1; 
end
set(handles.popupmenuProduct,'Value',j);

%Region
[N1,N2]=size(List.DataSetBranch);
jj=[];
for k=1:N1
    if strcmp(List.DataSetBranch{k,3},Level) & strcmp(List.DataSetBranch{k,4},Product)
        Region=List.DataSetBranch{k,5};
        j = strmatch(Region,List.Region(:,1),'exact');
        if ~isempty(j)
            jj=[jj ;j];
        else
            disp('warning: Region is not on the list')
        end
    end
end
jj=unique(sort(jj)); % sort and remove duplicate entries
RegionL=List.Region(jj,:);
set(handles.popupmenuRegion,'String',RegionL(:,2)); % long name
Region=List.DataSetBranch{jds,5};
j = strmatch(Region,RegionL(:,1),'exact');
if isempty(j) 
    j=1; 
end
set(handles.popupmenuRegion,'Value',j);


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21


