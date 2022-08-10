function varargout = gui_pathfinder4km(varargin)
% Boots-up Pathfinder4km GUI, interacts with user via interface.
% gui_pathfinder4km M-file for gui_pathfinder4km.fig
%
% INPUT(S):
%
%   [none]
%
% OUTPUT(S)
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
% 2007/06/01 created
% 2008/06/15 modified for compiled version, ceb
%==========================================================================

warning('off','MATLAB:dispatcher:InexactMatch');
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_pathfinder4km_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_pathfinder4km_OutputFcn, ...
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

% --- Executes just before gui_pathfinder4km is made visible.
function gui_pathfinder4km_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_pathfinder4km (see VARARGIN)

DSI=load('dsi_pathfinder4km.mat');
handles.DSI=DSI;

% Define handles "openfigure" field.
if ~isfield(handles,'openfigure')
    handles.openfigure = 0; %set equal to zero on the first time.
end

% Set interface name.
filename = 'get_pathfinder4km';
txt = ['NASA/REASoN Ocean Data Portal - '];
tmp1 = ['Pathfinder4km GUI'];
tmp2 = '';
txt = [txt,tmp1];
if (~isdeployed)
    tmp2 = strrep(op_version('mfile',filename),'Version ','');
    txt = [txt,' (Ver ',tmp2,')'];
end
set(hObject,'Name',txt);
set(handles.radiobuttonTimeSeries,'Enable','On');

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

[img,MAP]=imread('logo_nasa1.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonNASA,'CData',img)
pos = get(handles.pushbuttonNASA,'position'); %units will be in characters
set(handles.pushbuttonNASA,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_opendap.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonOpendap,'CData',img)
pos = get(handles.pushbuttonOpendap,'position'); %units will be in characters
set(handles.pushbuttonOpendap,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_noaa_nodc.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;
CY=(NY+12)/ppcy;
set(handles.pushbuttonNODCNOAA,'CData',img)
pos = get(handles.pushbuttonNODCNOAA,'position'); %units will be in characters
set(handles.pushbuttonNODCNOAA,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_podaac.jpg');
[NY,NX,N3]=size(img);
CX=(NX+4)/ppcx;
CY=(NY+4)/ppcy;
set(handles.pushbuttonjplnasa,'CData',img)
pos = get(handles.pushbuttonjplnasa,'position'); %units will be in characters
set(handles.pushbuttonjplnasa,'position',[pos(1),pos(2),CX,CY]);

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

if get(handles.radiobuttonClimatology,'Value')==0 & get(handles.radiobuttonTimeSeries,'Value')==0
set(handles.radiobuttonClimatology,'Value',1)
end

if get(handles.radiobuttonInWksp,'Value')==0 & get(handles.radiobuttonToFiles,'Value')==0
set(handles.radiobuttonInWksp,'Value',1)
end

% Make 'Load Last Request' button inactive when file gui_xxx_request.mat is
% absent.
if ~exist('gui_pathfinder4km_request.mat','file') 
    set(handles.pushbuttonLoadLastRequest,'Enable','Off')
end

% Climatology vs Time Series
ClimTimeSeriesCheck(hObject, eventdata, handles)
QualMaskOnOff(hObject, eventdata, handles)
SaveFilesOnOff(hObject, eventdata, handles)

InterimDataCheck(hObject, eventdata, handles)

set(handles.pushbuttonGetData,'String','Get Data','BackgroundColor',[0 1. 0.25])
set(handles.pushbuttonUpdateFileList,'Visible','Off','Enable','Off') % disable this button

% Update handles "openfigure" field.
handles.openfigure = 1;

% Choose default command line output for gui_pathfinder4km
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_pathfinder4km wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_pathfinder4km_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function radiobuttonClimatology_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonClimatology,'Value') == 0
    set(handles.radiobuttonTimeSeries,'Value',1)
else
    set(handles.radiobuttonTimeSeries,'Value',0)
end
ClimTimeSeriesCheck(hObject, eventdata, handles)
InterimDataCheck(hObject, eventdata, handles)

function radiobuttonTimeSeries_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonTimeSeries,'Value') == 0
    set(handles.radiobuttonClimatology,'Value',1)
else
    set(handles.radiobuttonClimatology,'Value',0)
end
ClimTimeSeriesCheck(hObject, eventdata, handles)
InterimDataCheck(hObject, eventdata, handles)

function checkbox1_Callback(hObject, eventdata, handles)
function checkbox2_Callback(hObject, eventdata, handles)
QualMaskOnOff(hObject, eventdata, handles)
function checkbox3_Callback(hObject, eventdata, handles)
function checkbox4_Callback(hObject, eventdata, handles)
function checkbox5_Callback(hObject, eventdata, handles)
function checkbox6_Callback(hObject, eventdata, handles)
QualMaskOnOff(hObject, eventdata, handles)
function checkbox7_Callback(hObject, eventdata, handles)
  if get(handles.checkbox7,'Value') == 1 % if requesting anomaly,
     set(handles.checkbox5,'Value',1)    % request sst as well
  end
function checkbox8_Callback(hObject, eventdata, handles)
QualMaskOnOff(hObject, eventdata, handles)
function checkbox9_Callback(hObject, eventdata, handles)
function checkbox10_Callback(hObject, eventdata, handles)
function checkbox11_Callback(hObject, eventdata, handles)
function checkbox12_Callback(hObject, eventdata, handles)
function checkbox13_Callback(hObject, eventdata, handles)
function checkbox14_Callback(hObject, eventdata, handles)

function editTime1yyyy_CreateFcn(hObject, eventdata, handles)
function editTime1yyyy_Callback(hObject, eventdata, handles)
function editTime1mm_CreateFcn(hObject, eventdata, handles)
function editTime1mm_Callback(hObject, eventdata, handles)
function editTime1dd_CreateFcn(hObject, eventdata, handles)
function editTime1dd_Callback(hObject, eventdata, handles)
function editTime2yyyy_CreateFcn(hObject, eventdata, handles)
function editTime2yyyy_Callback(hObject, eventdata, handles)
InterimDataCheck(hObject, eventdata, handles)
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
function editFileNamePrefix_Callback(hObject, eventdata, handles)
function editFileNamePrefix_CreateFcn(hObject, eventdata, handles)

function popupmenuQualityMask_CreateFcn(hObject, eventdata, handles)
function popupmenuQualityMask_Callback(hObject, eventdata, handles)
function popupmenu2_CreateFcn(hObject, eventdata, handles)
function popupmenu2_Callback(hObject, eventdata, handles)
ClimTimeSeriesCheck(hObject, eventdata, handles)
InterimDataCheck(hObject, eventdata, handles)
function popupmenu3_CreateFcn(hObject, eventdata, handles)
function popupmenu3_Callback(hObject, eventdata, handles)
ClimTimeSeriesCheck(hObject, eventdata, handles)

function radiobutton3_Callback(hObject, eventdata, handles)
Pass{1}='y';
InterimDataCheck(hObject, eventdata, handles)
%function radiobutton4_Callback(hObject, eventdata, handles)
%Pass{2}='y';
%InterimDataCheck(hObject, eventdata, handles)
function radiobutton5_Callback(hObject, eventdata, handles)
Pass{3}='y';
InterimDataCheck(hObject, eventdata, handles)
%function radiobutton6_Callback(hObject, eventdata, handles)
%Pass{4}='y';
%InterimDataCheck(hObject, eventdata, handles)

function pushbuttonSelectAll_Callback(hObject, eventdata, handles)
if strcmp(get(handles.radiobutton3,'Visible'),'on'); set(handles.radiobutton3,'Value',1); Pass{1}='y'; end
%if strcmp(get(handles.radiobutton4,'Visible'),'on'); set(handles.radiobutton4,'Value',1); Pass{2}='y'; end
if strcmp(get(handles.radiobutton5,'Visible'),'on'); set(handles.radiobutton5,'Value',1); Pass{3}='y'; end
%if strcmp(get(handles.radiobutton6,'Visible'),'on'); set(handles.radiobutton6,'Value',1); Pass{4}='y'; end
InterimDataCheck(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%check lines 351-365
S=get(handles.pushbuttonSelectAll,'String');
if S(1) == 'S' % Select All
set(handles.radiobutton3,'Value',1)
%set(handles.radiobutton4,'Value',1) 
set(handles.radiobutton5,'Value',1)
%set(handles.radiobutton6,'Value',1)
set(handles.pushbuttonSelectAll,'String','Clear All')
elseif S(1) == 'C' % Clear All
set(handles.radiobutton3,'Value',0)
%set(handles.radiobutton4,'Value',0) 
set(handles.radiobutton5,'Value',0)
%set(handles.radiobutton6,'Value',0)
set(handles.pushbuttonSelectAll,'String','Select All')
end 
%%%%%%%%%%%%%%%%%%%

function InterimDataCheck(hObject, eventdata, handles)
t0=get(handles.radiobuttonClimatology,'Value'); % Climatology
t1=get(handles.radiobutton3,'Value');
t2=get(handles.radiobutton4,'Value');
t3=get(handles.radiobutton5,'Value');
t4=get(handles.radiobutton6,'Value');
Year1=str2num(deblank(get(handles.editTime1yyyy,'String')));
Year2=str2num(deblank(get(handles.editTime2yyyy,'String')));
S='';
if Year2>=2007 & t0==0
S='From 2007-01-01 data are interim products and will be replaced in the future; from 2009-01-01 data are not quality controlled';
%set(handles.radiobutton3,'ForegroundColor',[1. 0. 0.])
%set(handles.radiobutton5,'ForegroundColor',[1. 0. 0.])
end
TAvg=get(handles.popupmenu2,'String');
j=get(handles.popupmenu2,'Value');
TAvg=deblank(TAvg{j});

if Year1<=2005 & Year2>=2005 & t0==0 & strcmp(TAvg,'yearly')
S='For 2005 yearly averages are not defined in Version5.0 due to switch between passes 2,4 and 1,3; GUI offers interim data';
end

set(handles.textInterimProduct,'String',S)

function ClimTimeSeriesCheck(hObject, eventdata, handles)
if get(handles.radiobuttonClimatology,'Value') == 0 % Time Series
        
    set(handles.checkbox1,'Enable','Off')
    set(handles.checkbox2,'Enable','Off')
    set(handles.checkbox3,'Enable','Off')
    set(handles.checkbox4,'Enable','Off')
    
    set(handles.checkbox5,'Enable','On')
    set(handles.checkbox6,'Enable','On')
    set(handles.checkbox7,'Enable','On')
    set(handles.checkbox8,'Enable','On')
    set(handles.checkbox9,'Enable','On')
    set(handles.checkbox10,'Enable','On')
    set(handles.checkbox11,'Enable','On')
    set(handles.checkbox12,'Enable','On')
    set(handles.checkbox13,'Enable','On')
    set(handles.checkbox14,'Enable','On')
    
    if get(handles.checkbox8,'Value') == 1
    set(handles.popupmenuQualityMask,'BackgroundColor',[1. 1. 1.],'Enable','On')
    end
  
    set(handles.editTime1yyyy,'Visible','On','Enable','On')
    set(handles.editTime2yyyy,'Visible','On','Enable','On')
    set(handles.text12,'Visible','On')
    set(handles.text15,'Visible','On')
    set(handles.textTimeRange,'Visible','On')
    set(handles.popupmenu3,'Visible','Off','Enable','Off')
    set(handles.popupmenu2,'Visible','On','Enable','On')
    
    %set(handles.radiobutton4,'Visible','On','Enable','On')
    %set(handles.radiobutton6,'Visible','On','Enable','On')
%    set(handles.textTimeRange,'String',['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2,'yyyy-mm-dd') ']'])
R.DataSetBranch='Pathfinder4km_daily';
%S=['load dsi_pathfinder4km.mat DSI_' R.DataSetBranch];
%eval(['load dsi_pathfinder4km.mat DSI_' R.DataSetBranch])
eval(['DSI=handles.DSI.DSI_' R.DataSetBranch ';'])
%eval(['clear DSI_' R.DataSetBranch ';'])

set(handles.textTimeRange,'String',...
['Available time range: [' datestr(DSI.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Range.Time2,'yyyy-mm-dd') ']'])

    %set(handles.textTimeRange,'String',['Available time range: [' datestr(datenum(1981,1,236),'yyyy-mm-dd') ' to ' datestr(datenum(2009,1,122),'yyyy-mm-dd') ']' ])


else % Climatology
    set(handles.checkbox1,'Enable','On')
    set(handles.checkbox2,'Enable','On')
    set(handles.checkbox3,'Enable','On')
    set(handles.checkbox4,'Enable','On')
        
    set(handles.checkbox5,'Enable','Off')
    set(handles.checkbox6,'Enable','Off')
    set(handles.checkbox7,'Enable','Off')
    set(handles.checkbox8,'Enable','Off')
    set(handles.checkbox9,'Enable','Off')
    set(handles.checkbox10,'Enable','Off')
    set(handles.checkbox11,'Enable','Off')
    set(handles.checkbox12,'Enable','Off')
    set(handles.checkbox13,'Enable','Off')
    set(handles.checkbox14,'Enable','Off')
    
    set(handles.popupmenuQualityMask,'BackgroundColor',[0.878 0.875 0.891],'Enable','Off')
    
    set(handles.editTime1yyyy,'Visible','Off','Enable','Off')
    set(handles.editTime2yyyy,'Visible','Off','Enable','Off')
    set(handles.text12,'Visible','Off')
    set(handles.text15,'Visible','Off')
    set(handles.textTimeRange,'Visible','Off')
    
    set(handles.popupmenu3,'Visible','On','Enable','On')
    set(handles.popupmenu2,'Visible','Off','Enable','Off')
    
    %set(handles.radiobutton4,'Value',0,'Visible','Off','Enable','Off')
    %set(handles.radiobutton6,'Value',0,'Visible','Off','Enable','Off')
end
if get(handles.radiobuttonClimatology,'Value')==1
S=get(handles.popupmenu3,'String');j=get(handles.popupmenu3,'Value');S=deblank(S{j});
%set(handles.text42,'Visible','Off')
%set(handles.text43,'String','1985-2001')
else 
S=get(handles.popupmenu2,'String');j=get(handles.popupmenu2,'Value');S=deblank(S{j});
%set(handles.text42,'Visible','On')
%set(handles.text43,'String','1985-2002')
end
TXT=['Temporal resolution is ' lower(S)];
set(handles.textTemporalResolution,'String',TXT,'ForegroundColor',[0 0 1])

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
% set(handles.pushbuttonGetData,'String','For a new request reissue Get_Pathfinder')
% uiresume(handles.figure1)

% Turn off ability to re-issue "Get Data."
set(hObject,'BackgroundColor',[0.878 0.875 0.891],'Enable','Off','String','Please wait...');
pause(.1)

% Get input from user.
request=InputGet(hObject, eventdata, handles);

% Save request to .mat file.
filename = 'gui_pathfinder4km_request.mat';
%if (~isdeployed)
    pname = mfilename('fullpath');
    pname = pname(1:length(pname)-length(mfilename)); %mfilename does not have .m in it
    filename = [pname,filename];
%end
save(filename,'request','-mat')
set(handles.pushbuttonLoadLastRequest,'Enable','On')

% If it passes inspection, send down to user's workspace
% as variable "request". Then the user can get the data
%assignin('base','request',request);

get_pathfinder4km(request);

set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25])

function pushbuttonNODCNOAA_Callback(hObject, eventdata, handles)
web('http://www.nodc.noaa.gov/sog/pathfinder4km/userguide.html');
%web 'http://www.nodc.noaa.gov/sog/pathfinder4km/userguide.html' -browser
function pushbuttonOpendap_Callback(hObject, eventdata, handles)
web('http://www.opendap.org/');
%web 'http://www.opendap.org/' -browser
function pushbuttonjplnasa_Callback(hObject, eventdata, handles)
web('http://dods.jpl.nasa.gov/');
%web 'http://dods.jpl.nasa.gov/' -browser
function pushbuttonNASA_Callback(hObject, eventdata, handles)
web('http://www.nasa.gov');
%web 'http://www.nasa.gov' -browser

function pushbuttonInfo_Callback(hObject, eventdata, handles)
h = op_timemessage('polar',0);

function request=InputGet(hObject, eventdata, handles)

S1=deblank(get(handles.editTime1yyyy,'String'));
S2=deblank(get(handles.editTime1mm,'String'));
%if length(S2) == 1 S2=['0' S2]; end
S3=deblank(get(handles.editTime1dd,'String'));
%if length(S3) == 1 S3=['0' S3]; end
S4=deblank(get(handles.editTime2yyyy,'String'));
S5=deblank(get(handles.editTime2mm,'String'));
%if length(S5) == 1 S5=['0' S5]; end
S6=deblank(get(handles.editTime2dd,'String'));
%if length(S6) == 1 S6=['0' S6]; end

R.DATE1='00000000';
R.DATE1(4-length(S1)+1:4)=S1;
R.DATE1(6-length(S2)+1:6)=S2;
R.DATE1(8-length(S3)+1:8)=S3;
R.DATE2='00000000';
R.DATE2(4-length(S4)+1:4)=S4;
R.DATE2(6-length(S5)+1:6)=S5;
R.DATE2(8-length(S6)+1:8)=S6;

%DATE1=[S1 S2 S3];
%DATE2=[S4 S5 S6];
%DATEINCR=1;
R.DATEINCR=str2num(deblank(get(handles.editTimeStep,'String')));

R.LAT1=str2num(deblank(get(handles.editLat1,'String')));
R.LAT2=str2num(deblank(get(handles.editLat2,'String')));
R.LON1=str2num(deblank(get(handles.editLon1,'String')));
R.LON2=str2num(deblank(get(handles.editLon2,'String')));
R.LATINCR=str2num(deblank(get(handles.editLatStep,'String')));
R.LONINCR=str2num(deblank(get(handles.editLonStep,'String')));

R.DIRNAME =deblank(get(handles.editDirectoryName,'String'));
if isempty(R.DIRNAME) R.DIRNAME=''; end
R.FNPREFIX=deblank(get(handles.editFileNamePrefix,'String'));

R.Temporal='Climatology';
if get(handles.radiobuttonClimatology,'Value')==0 R.Temporal='TimeSeries'; end
if get(handles.radiobuttonTimeSeries,'Value')==0 R.Temporal='Climatology'; end


%R.Clim_Fields={'Clim_SST','n';'Clim_StandardDeviation','n';'Clim_Counts','n'};
%R.Fields={'sst','n';'bsst','n';'sdev','n';'mask1','n';'mask2','n';'num','n';'qual','n'};
if strcmp(R.Temporal,'Climatology')
R.Fields={};
if get(handles.checkbox1,'Value')==1 R.Fields{end+1}='Clim_SST'; end
if get(handles.checkbox3,'Value')==1 R.Fields{end+1}='Clim_StandardDeviation'; end
if get(handles.checkbox4,'Value')==1 R.Fields{end+1}='Clim_Counts'; end
else
R.Fields={};
if get(handles.checkbox9 ,'Value')==1 R.Fields{end+1}='qual'; end
if get(handles.checkbox5 ,'Value')==1 R.Fields{end+1}='sst'; end
if get(handles.checkbox10,'Value')==1 R.Fields{end+1}='bsst'; end
if get(handles.checkbox11,'Value')==1 R.Fields{end+1}='sdev'; end
if get(handles.checkbox12,'Value')==1 R.Fields{end+1}='mask1'; end
if get(handles.checkbox13,'Value')==1 R.Fields{end+1}='num'; end
if get(handles.checkbox14,'Value')==1 R.Fields{end+1}='mask2'; end
end
R.Coordinates={'latitude','longitude','time'};

R.ApplyLandMask='n';
if strcmp(R.Temporal,'Climatology') 
if get(handles.checkbox2 ,'Value')==1 R.ApplyLandMask='y'; end
else
if get(handles.checkbox6 ,'Value')==1 R.ApplyLandMask='y'; end
end

R.ApplyQualMask='n';
if get(handles.checkbox8 ,'Value')==1 
S=get(handles.popupmenuQualityMask,'String');j=get(handles.popupmenuQualityMask,'Value');
R.ApplyQualMask=S{j}(3); %'n','4','7'
end

R.CalcAnomaly='n';
if get(handles.checkbox7 ,'Value')==1 R.CalcAnomaly='y'; end

if strcmp(R.Temporal,'Climatology') 
S=get(handles.popupmenu3,'String');j=get(handles.popupmenu3,'Value');
R.TAvg=deblank(S{j});
else
S=get(handles.popupmenu2,'String');j=get(handles.popupmenu2,'Value');
R.TAvg=deblank(S{j});
end

if strcmp(R.Temporal,'Climatology')
    R.DATE1(1:4)='0000';R.DATE2(1:4)='0000';
    % if climatology request spans a new year, as in Dec-Jan period
    if str2num(R.DATE2)<str2num(R.DATE1) R.DATE2(1:4)='0001'; end 
end      

R.Passes={}; 
if get(handles.radiobutton3,'Value') == 1 R.Passes{end+1}='night'; end
%if get(handles.radiobutton4,'Value') == 1 R.Pass{2}='y'; end
if get(handles.radiobutton5,'Value') == 1 R.Passes{end+1}='day'; end
%if get(handles.radiobutton6,'Value') == 1 R.Pass{4}='y'; end
    

R.SaveWorkspace='n'; 
if get(handles.radiobuttonInWksp,'Value') == 1 R.SaveWorkspace='y'; end
R.SaveFiles='n'; 
if get(handles.radiobuttonToFiles,'Value') == 1 R.SaveFiles='y'; end

S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');
R.SaveMode=deblank(S{j});

R.DataSetName='Pathfinder4km';
if strcmp(R.Temporal,'Climatology')
    R.DataSetBranch='Pathfinder4km_Climatology';
else
R.DataSetBranch=['Pathfinder4km_' R.TAvg];
end
%% Get default command line output from handles structure.   
% R.RequestDate=datestr(now,'yyyy-mm-dd HH:MM:SS');
request=R;
handles.request=request;

% Load Last request
function pushbuttonLoadLastRequest_Callback(hObject, eventdata, handles)

filename = 'gui_pathfinder4km_request.mat';
if exist(filename,'file')
    load(filename)
else
    msgbox(['File ',filename,' does not exist.']);
    return
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
set(handles.editDirectoryName,'String',R.DIRNAME);
set(handles.editFileNamePrefix,'String',R.FNPREFIX);

if strcmp(R.Temporal,'TimeSeries')
    set(handles.radiobuttonClimatology,'Value',0)
    set(handles.radiobuttonTimeSeries,'Value',1)
else
    set(handles.radiobuttonClimatology,'Value',1)
    set(handles.radiobuttonTimeSeries,'Value',0)
end

set(handles.checkbox1,'Value',0); 
set(handles.checkbox3,'Value',0); 
set(handles.checkbox4,'Value',0); 
if ~isempty(strmatch('Clim_SST',R.Fields))               set(handles.checkbox1,'Value',1); end
if ~isempty(strmatch('Clim_StandardDeviation',R.Fields)) set(handles.checkbox3,'Value',1); end
if ~isempty(strmatch('Clim_Counts',R.Fields))            set(handles.checkbox4,'Value',1); end

set(handles.checkbox5 ,'Value',0); 
set(handles.checkbox9 ,'Value',0); 
set(handles.checkbox10,'Value',0); 
set(handles.checkbox11 ,'Value',0);
set(handles.checkbox12,'Value',0); 
set(handles.checkbox13,'Value',0); 
set(handles.checkbox14,'Value',0); 
if ~isempty(strmatch('sst',R.Fields)) set(handles.checkbox5,'Value',1); end
if ~isempty(strmatch('qual',R.Fields)) set(handles.checkbox9,'Value',1); end
if ~isempty(strmatch('bsst',R.Fields)) set(handles.checkbox10,'Value',1); end
if ~isempty(strmatch('sdev',R.Fields)) set(handles.checkbox11,'Value',1); end
if ~isempty(strmatch('mask1',R.Fields)) set(handles.checkbox12,'Value',1); end
if ~isempty(strmatch('num',R.Fields)) set(handles.checkbox13,'Value',1); end
if ~isempty(strmatch('mask2',R.Fields)) set(handles.checkbox14,'Value',1); end


set(handles.checkbox2 ,'Value',0);
set(handles.checkbox6 ,'Value',0);
if R.ApplyLandMask=='y' & get(handles.radiobuttonClimatology,'Value')==1 set(handles.checkbox2 ,'Value',1); end
if R.ApplyLandMask=='y' & get(handles.radiobuttonTimeSeries,'Value')==1 set(handles.checkbox6 ,'Value',1); end

%ApplyQualMask=S{j}(3); %'n','4','7'
R.ApplyQualMask;
S=get(handles.popupmenuQualityMask,'String');
S=strvcat(S);
j=strfind(S(:,3)',R.ApplyQualMask);
if ~isempty(j) 
    set(handles.checkbox8 ,'Value',1); 
    set(handles.popupmenuQualityMask,'Value',j);
else
    set(handles.checkbox8 ,'Value',0);
end

set(handles.checkbox7 ,'Value',0);
if R.CalcAnomaly=='y' set(handles.checkbox7 ,'Value',1); end

if strcmp(R.Temporal,'Climatology')
S=get(handles.popupmenu3,'String');
j=find(strcmp(S,R.TAvg));
if isempty(j) j=1; end
set(handles.popupmenu3,'Value',j);
else
S=get(handles.popupmenu2,'String');
j=find(strcmp(S,R.TAvg));
if isempty(j) j=1; end
set(handles.popupmenu2,'Value',j);
end
ClimTimeSeriesCheck(hObject, eventdata, handles)

set(handles.radiobutton3,'Value',0);
%set(handles.radiobutton4,'Value',0);
set(handles.radiobutton5,'Value',0);
%set(handles.radiobutton6,'Value',0);
%if R.Pass{1}=='y' set(handles.radiobutton3,'Value',1); end
%if R.Pass{2}=='y' set(handles.radiobutton4,'Value',1); end
%if R.Pass{3}=='y' set(handles.radiobutton5,'Value',1); end
%if R.Pass{4}=='y' set(handles.radiobutton6,'Value',1); end
if ~isempty(strmatch('night',R.Passes)) set(handles.radiobutton3,'Value',1); end
if ~isempty(strmatch('day',  R.Passes)) set(handles.radiobutton5,'Value',1); end


%if ~isempty(strmatch('night',R.Passes)) set(handles.radiobuttonPassNight,'Value',1); end
%if ~isempty(strmatch('day',  R.Passes)) set(handles.radiobuttonPassDay,'Value',1);   end


set(handles.radiobuttonInWksp,'Value',0);
set(handles.radiobuttonToFiles,'Value',0);
if R.SaveWorkspace=='y' set(handles.radiobuttonInWksp,'Value',1); end
if R.SaveFiles=='y' set(handles.radiobuttonToFiles,'Value',1); end

S=get(handles.popupmenuMode,'String');
j=find(strcmp(S,R.SaveMode));
if isempty(j) j=1; end
set(handles.popupmenuMode,'Value',j);

ClimTimeSeriesCheck(hObject, eventdata, handles)
QualMaskOnOff(hObject, eventdata, handles)
InterimDataCheck(hObject, eventdata, handles)
S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');
R.SaveMode=S{j};
% disable saving to workspace for netcdf mode
set(handles.radiobuttonInWksp,'Enable','On');
if strcmp(R.SaveMode,'netcdf')
    set(handles.radiobuttonInWksp,'Value',0,'Enable','Off');
    set(handles.radiobuttonToFiles,'Value',1);
end
SaveFilesOnOff(hObject, eventdata, handles)

function SaveFilesOnOff(hObject, eventdata, handles)   
% Set color of data documentation button to match figure color.
%BGColor = get(hObject,'Color'); %here hObject is the current figure
BGColor=[0.878 0.875 0.891];
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',BGColor,'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',BGColor,'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end

function QualMaskOnOff(hObject, eventdata, handles)
% Set color of data documentation button to match figure color.
%BGColor = get(hObject,'Color'); %here hObject is the current figure
BGColor=[0.878 0.875 0.891];
if     get(handles.checkbox8,'Value') == 0
   set(handles.popupmenuQualityMask,'BackgroundColor',BGColor,'Enable','Off')
elseif get(handles.checkbox8,'Value') == 1
   set(handles.popupmenuQualityMask,'BackgroundColor',[1.0 1.0 1.0],'Enable','On')
end

S='';
if get(handles.checkbox2,'Value') == 1 
S='A new variable Clim_SST_masked will be created'; end
set(handles.textClimMask,'String',S)

S='';
if get(handles.checkbox8,'Value') == 1 | get(handles.checkbox6,'Value') == 1
S='A new variable sst_masked and/or ssta_masked will be created'; end
set(handles.textTimeSeriesMask,'String',S)

if get(handles.checkbox8,'Value') == 1 % if applyiny quality mask,
     set(handles.checkbox9,'Value',1)    % request qual
  end


% Update File list
function pushbuttonUpdateFileList_Callback(hObject, eventdata, handles)
filename = strcat('gui_pathfinder4km_filelist','.mat');
ButtonName=questdlg(['You requested to update the list of files ',...
    'in the Pathfinder4km archive and overwrite ',...
    filename,'. This may take up to 2 hours.', ...
    'If update fails, restore the backup list from gui_pathfinder4km_filelist0.mat'], ...
    'Pathfinder4km Update File List','Continue','Cancel','Cancel');
if strcmp(ButtonName,'Cancel')
    return
end

