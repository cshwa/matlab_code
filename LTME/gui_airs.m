function varargout = gui_airs(varargin)
% Boots-up AIRS GUI, interacts with user via interface.
% gui_airs M-file for gui_airs.fig 
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
% Copyright 2007
% $Version 2.0.0$

%==========================================================================
% Meri Sheremet
% February 2008
%
% REVISION HISTORY:
% 2008/02/01 created, ms
% 2008/06/15 modified for compiled version, ceb
%==========================================================================

warning('off','MATLAB:dispatcher:InexactMatch');
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_airs_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_airs_OutputFcn, ...
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

% --- Executes just before gui_airs is made visible.
function gui_airs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_airs (see VARARGIN)
DSI=load('dsi_airs.mat');
handles.DSI=DSI;

% Define handles "openfigure" field.
if ~isfield(handles,'openfigure')
    handles.openfigure = 0; %set equal to zero on the first time.  
end %ceb 2008/02/28

% Set interface name.
filename = 'get_airs';
pname = strrep(mfilename('fullpath'),mfilename,'');
filename = [pname,filename];
txt = ['NASA/REASoN Ocean Data Portal - '];
tmp1 = ['AIRS GUI'];
tmp2 = '';
txt = [txt,tmp1];
if (~isdeployed)
%    op_version('mfile',filename)
    tmp2 = strrep(op_version('mfile',filename),'Version ','');
    txt = [txt,' (Ver ',tmp2,')'];
end
set(hObject,'Name',txt);

%set(handles.radiobuttonToFiles,'Enable','Off');

% % Set Version Number.
% filename = strrep(mfilename,'gui','get');
% info = op_version(filename);
% set(handles.text2,'String',info);

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
set(handles.pushbuttonOpendap,'CData',img)
pos = get(handles.pushbuttonOpendap,'position'); %units will be in characters
set(handles.pushbuttonOpendap,'position',[pos(1),pos(2),CX,CY]);  

[img,MAP]=imread('logo_nasa_gsfc.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonNASAGODDARDSFC,'CData',img)
pos = get(handles.pushbuttonNASAGODDARDSFC,'position'); %units will be in characters
set(handles.pushbuttonNASAGODDARDSFC,'position',[pos(1),pos(2),CX,CY]);  

[img,MAP]=imread('logo_nasa1.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx;         % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonNASA,'CData',img)
pos = get(handles.pushbuttonNASA,'position'); %units will be in characters
set(handles.pushbuttonNASA,'position',[pos(1),pos(2),CX,CY]); 

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

S=get(handles.popupmenuTemporalAverages,'String');j=get(handles.popupmenuTemporalAverages,'Value');S=deblank(S{j});
TXT=['Temporal resolution is ' S];
set(handles.textTemporalResolution,'String',TXT,'ForegroundColor',[0 0 1])

% Set ToolTipString.
txt=[' Atmospheric temperature profile in 24 standard pressure levels ',...
     ' from 1000 to 1.0 mbar (K).'];
op_tooltipwrap(handles.checkbox13,txt,38); 
 
txt=[' Microwave-only atmospheric temperature profile in 24 standard ',...
     ' pressure levels from 1000 to 1.0 mbar (K).'];
op_tooltipwrap(handles.checkbox14,txt,38);  

txt=[' Geopotential height at 24 standard pressure levels ',...
     ' from 1000 to 1.0 mbar.(m).'];
op_tooltipwrap(handles.checkbox15,txt,38); 

txt=[' MW-Only geopotential height in meters at 24 standard pressure levels ',...
     ' from 1000 to 1.0 mbar.(m).'];
op_tooltipwrap(handles.checkbox16,txt,38);

txt=[' Relative humidity profile in 12 standard pressure levels ',...
     ' from 1000 to 100 mbar (percent).'];
op_tooltipwrap(handles.checkbox17,txt,38);

txt=[' Water vapor mass mixing ratio at 12 standard pressure levels ',...
     ' from 1000 to 100 mbar (gm/kg dry air).'];
op_tooltipwrap(handles.checkbox18,txt,38);

txt=[' IR surface emissivity on a frequency grid ',...
     ' (832, 961, 1203, 2616 cm-1).'];
op_tooltipwrap(handles.checkbox19,txt,38);


txt=[' Microwave spectral emissivity on a frequency grid ',...
     ' (23.8, 50.3 and 89.0 GHz).'];
op_tooltipwrap(handles.checkbox20,txt,38);

CheckAvailableTimeRange(hObject, eventdata, handles)

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

% Make 'Load Last Request' button inactive when file gui_xxx_request.mat is
% absent.
if ~exist('gui_airs_request.mat','file') 
    set(handles.pushbuttonLoadLastRequest,'Enable','Off')
end

% Update handles "openfigure" field.
handles.openfigure = 1;
   
% Choose default command line output for gui_airs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_airs wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_airs_OutputFcn(hObject, eventdata, handles)
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
function checkbox20_Callback(hObject, eventdata, handles)
function checkbox21_Callback(hObject, eventdata, handles)
function checkbox22_Callback(hObject, eventdata, handles)
function checkbox23_Callback(hObject, eventdata, handles)

function pushbuttonSelectAll2D_Callback(hObject, eventdata, handles)
S=get(handles.pushbuttonSelectAll2D,'String');
if S(1) == 'S' % Select All
set(handles.checkbox1,'Value',1)
set(handles.checkbox2,'Value',1)
set(handles.checkbox3,'Value',1)
set(handles.checkbox4,'Value',1)
set(handles.checkbox5,'Value',1)
set(handles.checkbox6,'Value',1)
set(handles.checkbox7,'Value',1)
set(handles.checkbox8,'Value',1)
set(handles.checkbox9,'Value',1)
set(handles.checkbox10,'Value',1)
set(handles.checkbox11,'Value',1)
set(handles.checkbox12,'Value',1)
set(handles.pushbuttonSelectAll2D,'String','Clear All')
elseif S(1) == 'C' % Clear All
set(handles.checkbox1,'Value',0)
set(handles.checkbox2,'Value',0)
set(handles.checkbox3,'Value',0)
set(handles.checkbox4,'Value',0)
set(handles.checkbox5,'Value',0)
set(handles.checkbox6,'Value',0)
set(handles.checkbox7,'Value',0)
set(handles.checkbox8,'Value',0)
set(handles.checkbox9,'Value',0)
set(handles.checkbox10,'Value',0)
set(handles.checkbox11,'Value',0)
set(handles.checkbox12,'Value',0)
set(handles.pushbuttonSelectAll2D,'String','Select All')
end  

function pushbuttonSelectAll3D_Callback(hObject, eventdata, handles)
S=get(handles.pushbuttonSelectAll3D,'String');
if S(1) == 'S' % Select All
set(handles.checkbox13,'Value',1)
set(handles.checkbox14,'Value',1)
set(handles.checkbox15,'Value',1)
set(handles.checkbox16,'Value',1)
set(handles.checkbox17,'Value',1)
set(handles.checkbox18,'Value',1)
set(handles.checkbox19,'Value',1)
set(handles.checkbox20,'Value',1)
set(handles.pushbuttonSelectAll3D,'String','Clear All')
elseif S(1) == 'C' % Clear All
set(handles.checkbox13,'Value',0)
set(handles.checkbox14,'Value',0)
set(handles.checkbox15,'Value',0)
set(handles.checkbox16,'Value',0)
set(handles.checkbox17,'Value',0)
set(handles.checkbox18,'Value',0)
set(handles.checkbox19,'Value',0)
set(handles.checkbox20,'Value',0)
set(handles.pushbuttonSelectAll3D,'String','Select All')
end  

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

function radiobuttonPassDescending_Callback(hObject, eventdata, handles)
function radiobuttonPassAscending_Callback(hObject, eventdata, handles)

function popupmenuTemporalAverages_Callback(hObject, eventdata, handles)
CheckAvailableTimeRange(hObject, eventdata, handles)
S=get(handles.popupmenuTemporalAverages,'String');j=get(handles.popupmenuTemporalAverages,'Value');S=deblank(S{j});
TXT=['Temporal resolution is ' S];
set(handles.textTemporalResolution,'String',TXT,'ForegroundColor',[0 0 1])

function popupmenuTemporalAverages_CreateFcn(hObject, eventdata, handles)

function editLatStep_Callback(hObject, eventdata, handles)
function editLatStep_CreateFcn(hObject, eventdata, handles)
function editLonStep_Callback(hObject, eventdata, handles)
function editLonStep_CreateFcn(hObject, eventdata, handles)
function editTimeStep_Callback(hObject, eventdata, handles)
function editTimeStep_CreateFcn(hObject, eventdata, handles)

function editLevel1_Callback(hObject, eventdata, handles)
function editLevel1_CreateFcn(hObject, eventdata, handles)
function editLevel2_Callback(hObject, eventdata, handles)
function editLevel2_CreateFcn(hObject, eventdata, handles)
function editLevel3_Callback(hObject, eventdata, handles)
function editLevel3_CreateFcn(hObject, eventdata, handles)
function editLevel4_Callback(hObject, eventdata, handles)
function editLevel4_CreateFcn(hObject, eventdata, handles)

function editDirectoryName_Callback(hObject, eventdata, handles)
function editDirectoryName_CreateFcn(hObject, eventdata, handles)
function editFileNamePrefix_Callback(hObject, eventdata, handles)
function editFileNamePrefix_CreateFcn(hObject, eventdata, handles)

function popupmenuMode_Callback(hObject, eventdata, handles)
S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');
R.SaveMode=S{j};
% disable saving to workspace for netcdf mode
set(handles.radiobuttonInWksp,'Enable','On');
if strcmp(R.SaveMode,'netcdf')
    set(handles.radiobuttonInWksp,'Value',0,'Enable','Off');
    set(handles.radiobuttonToFiles,'Value',1);
end

function popupmenuMode_CreateFcn(hObject, eventdata, handles)

function radiobuttonInWksp_Callback(hObject, eventdata, handles)
function radiobuttonToFiles_Callback(hObject, eventdata, handles)
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end

%---Obtains or creates a new directory for output data files.
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

% Check user input. (no file present)

% Save request to .mat file.
filename = 'gui_airs_request.mat';
%if (~isdeployed)
pname = strrep(mfilename('fullpath'),mfilename,'');
filename = [pname,filename];
%end
save(filename,'request','-mat')
set(handles.pushbuttonLoadLastRequest,'Enable','On')

% If it passes inspection, send down to user's workspace
% as variable "request". Then the user can get the data
% assignin('base','request',request);

% Get data.
get_airs(request);

% Re-enable button after obtaining data.
set(hObject,'Enable','On','String','Get Data','BackgroundColor',[0 1. 0.25]);

function pushbuttonNASAGODDARDSFC_Callback(hObject, eventdata, handles)
web('http://disc.gsfc.nasa.gov/AIRS/airsL3_STD.shtml');
%web 'http://disc.gsfc.nasa.gov/AIRS/airsL3_STD.shtml' -browser
%web 'http://disc.gsfc.nasa.gov/AIRS/documentation/airs_instrument_guide.shtml' -browser
function pushbuttonOpendap_Callback(hObject, eventdata, handles)
web('http://www.opendap.org/');
%web 'http://www.opendap.org/' -browser
function pushbuttonNASA_Callback(hObject, eventdata, handles)
web('http://www.nasa.gov');
%web 'http://www.nasa.gov' -browser

function pushbuttonInfo_Callback(hObject, eventdata, handles)
h = op_timemessage('polar',0);

function request=InputGet(hObject, eventdata, handles)
% begining date
S1=deblank(get(handles.editTime1yyyy,'String'));
S2=deblank(get(handles.editTime1mm,'String'));
S3=deblank(get(handles.editTime1dd,'String'));
R.DATE1='00000000';
R.DATE1(4-length(S1)+1:4)=S1;
R.DATE1(6-length(S2)+1:6)=S2;
R.DATE1(8-length(S3)+1:8)=S3;
% end date
S1=deblank(get(handles.editTime2yyyy,'String'));
S2=deblank(get(handles.editTime2mm,'String'));
S3=deblank(get(handles.editTime2dd,'String'));
R.DATE2='00000000';
R.DATE2(4-length(S1)+1:4)=S1;
R.DATE2(6-length(S2)+1:6)=S2;
R.DATE2(8-length(S3)+1:8)=S3;

R.DATEINCR=1;

R.LAT1=str2num(deblank(get(handles.editLat1,'String')));
R.LAT2=str2num(deblank(get(handles.editLat2,'String')));
R.LON1=str2num(deblank(get(handles.editLon1,'String')));
R.LON2=str2num(deblank(get(handles.editLon2,'String')));
R.LATINCR=str2num(deblank(get(handles.editLatStep,'String')));
R.LONINCR=str2num(deblank(get(handles.editLonStep,'String')));

R.DIRNAME =deblank(get(handles.editDirectoryName,'String'));
if isempty(R.DIRNAME) R.DIRNAME=''; end
R.FNPREFIX=deblank(get(handles.editFileNamePrefix,'String'));

%define variablenames, checked, ndim, third dimension
R.Fields={'TotCldLiqH2O',  'n',   2,    ''              ,{'latitude','longitude'};
          'TotH2OVap',     'n',   2,    ''              ,{'latitude','longitude'};
          'TotH2OVap_MW',  'n',   2,    ''              ,{'latitude','longitude'};
          'SurfAirTemp',   'n',   2,    ''              ,{'latitude','longitude'};
          'SurfSkinTemp',  'n',   2,    ''              ,{'latitude','longitude'};
          'SurfPres',      'n',   2,    ''              ,{'latitude','longitude'};
          'TotO3',         'n',   2,    ''              ,{'latitude','longitude'};
          'CloudFrc',      'n',   2,    ''              ,{'latitude','longitude'};
          'CloudTopPres',  'n',   2,    ''              ,{'latitude','longitude'};
          'CloudTopTemp',  'n',   2,    ''              ,{'latitude','longitude'};
          'OLR',           'n',   2,    ''              ,{'latitude','longitude'};
          'ClrOLR',        'n',   2,    ''              ,{'latitude','longitude'};
          'Temperature',   'n',   3,    'TempPresLvls'   ,{'latitude','longitude','TempPresLvls'};
          'Temperature_MW','n',   3,    'TempPresLvls'   ,{'latitude','longitude','TempPresLvls'};
          'GPHeight',      'n',   3,    'TempPresLvls'   ,{'latitude','longitude','TempPresLvls'};
          'GPHeight_MW'    'n',   3,    'TempPresLvls'   ,{'latitude','longitude','TempPresLvls'};
          'RelHumid',      'n',   3,    'H2OPresLvls'    ,{'latitude','longitude','H2OPresLvls'};
          'H2OVapMMR',     'n',   3,    'H2OPresLvls'    ,{'latitude','longitude','H2OPresLvls'};
          'EmisIR',        'n',   3,    'IREmisFreqs'    ,{'latitude','longitude','IREmisFreqs'};
          'EmisMW_MW',     'n',   3,    'MWEmisFreqs'    ,{'latitude','longitude','MWEmisFreqs'};
          };
                 
% select checked variables by setting 'y' in the second column of 
% Fields cell array
if get(handles.checkbox1,'Value')==1 R.Fields{1,2}='y'; end
if get(handles.checkbox2,'Value')==1 R.Fields{2,2}='y'; end
if get(handles.checkbox3,'Value')==1 R.Fields{3,2}='y'; end
if get(handles.checkbox4,'Value')==1 R.Fields{4,2}='y'; end
if get(handles.checkbox5,'Value')==1 R.Fields{5,2}='y'; end
if get(handles.checkbox6,'Value')==1 R.Fields{6,2}='y'; end
if get(handles.checkbox7,'Value')==1 R.Fields{7,2}='y'; end
if get(handles.checkbox8,'Value')==1 R.Fields{8,2}='y'; end
if get(handles.checkbox9,'Value')==1 R.Fields{9,2}='y'; end
if get(handles.checkbox10,'Value')==1 R.Fields{10,2}='y'; end
if get(handles.checkbox11,'Value')==1 R.Fields{11,2}='y'; end
if get(handles.checkbox12,'Value')==1 R.Fields{12,2}='y'; end
if get(handles.checkbox13,'Value')==1 R.Fields{13,2}='y'; end
if get(handles.checkbox14,'Value')==1 R.Fields{14,2}='y'; end
if get(handles.checkbox15,'Value')==1 R.Fields{15,2}='y'; end
if get(handles.checkbox16,'Value')==1 R.Fields{16,2}='y'; end
if get(handles.checkbox17,'Value')==1 R.Fields{17,2}='y'; end
if get(handles.checkbox18,'Value')==1 R.Fields{18,2}='y'; end
if get(handles.checkbox19,'Value')==1 R.Fields{19,2}='y'; end
if get(handles.checkbox20,'Value')==1 R.Fields{20,2}='y'; end

% problem CloudFrcVis descending pass is missing in dataset - disabled
%R.Fields{10,2}='n';

%FieldsAnc={'sdev','n';'ct','n';'LandMask','n'};
R.FieldsAnc={'sdev',  'n';
           'ct',    'n';
           'LandSeaMask','n';
           'TotalCounts','n'};
if get(handles.checkbox21,'Value')==1 R.FieldsAnc{1,2}='y'; end
if get(handles.checkbox22,'Value')==1 R.FieldsAnc{2,2}='y'; end
if get(handles.checkbox23,'Value')==1 R.FieldsAnc{3,2}='y'; end

R.Pass={'Descending','n';'Ascending','n'};
if get(handles.radiobuttonPassDescending,'Value') == 1 R.Pass{1,2}='y'; end
if get(handles.radiobuttonPassAscending,'Value') == 1 R.Pass{2,2}='y'; end

R.DATEINCR=str2num(deblank(get(handles.editTimeStep,'String')));

%get level constraints and add brakets [ ] if necessary
S=deblank(get(handles.editLevel1,'String'));
if ~isempty(S)&S(1)~='[' S=['[' S]; end; 
if ~isempty(S)&S(end)~=']' S=[S ']']; end;
R.TEMPPRESLVLS=S;
S=deblank(get(handles.editLevel2,'String'));
if ~isempty(S)&S(1)~='[' S=['[' S]; end; 
if ~isempty(S)&S(end)~=']' S=[S ']']; end;
R.H2OPRESLVLS=S;
S=deblank(get(handles.editLevel3,'String'));
if ~isempty(S)&S(1)~='[' S=['[' S]; end; 
if ~isempty(S)&S(end)~=']' S=[S ']']; end;
R.IREMISFREQS=S;
S=deblank(get(handles.editLevel4,'String'));
if ~isempty(S)&S(1)~='[' S=['[' S]; end; 
if ~isempty(S)&S(end)~=']' S=[S ']']; end;
R.MWEMISFREQS=S;

S=get(handles.popupmenuTemporalAverages,'String');j=get(handles.popupmenuTemporalAverages,'Value');
R.TAVG=deblank(S{j});

R.DIRNAME =deblank(get(handles.editDirectoryName,'String'));
if isempty(R.DIRNAME) R.DIRNAME=''; end
R.FNPREFIX=deblank(get(handles.editFileNamePrefix,'String'));

R.SaveWorkspace='n'; 
if get(handles.radiobuttonInWksp,'Value') == 1 R.SaveWorkspace='y'; end
R.SaveFiles='n'; 
if get(handles.radiobuttonToFiles,'Value') == 1 R.SaveFiles='y'; end

LEVEL=1;

S=get(handles.popupmenuMode,'String');j=get(handles.popupmenuMode,'Value');
R.SaveMode=deblank(S{j});

%R.RequestDate=datestr(now,'yyyy-mm-dd HH:MM:SS');
request=R;
handles.request=request;

%---Load the last request and set uicontrol objects.
function pushbuttonLoadLastRequest_Callback(hObject, eventdata, handles)

filename = 'gui_airs_request.mat';
if exist(filename,'file')
    load(filename)
else
    msgbox(['File ',filename,' does not exist.']);
    return
end

R=request;

% select only checked variables
[L1,L2]=size(R.Fields);
for k=1:L1
    eval(['set(handles.checkbox' num2str(k) ',''Value'',0)'])
    if R.Fields{k,2}=='y';
        eval(['set(handles.checkbox' num2str(k) ',''Value'',1)']) 
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
set(handles.editLevel1,'String',R.TEMPPRESLVLS);
set(handles.editLevel2,'String',R.H2OPRESLVLS);
set(handles.editLevel3,'String',R.IREMISFREQS);
set(handles.editLevel4,'String',R.MWEMISFREQS);

%FieldsAnc={'sdev',  'n';
%           'ct',    'n';
%           'LandSeaMask','n';
%           'TotalCounts','n'};
set(handles.checkbox21,'Value',0);
if R.FieldsAnc{1,2}=='y' set(handles.checkbox21,'Value',1); end
set(handles.checkbox22,'Value',0);
if R.FieldsAnc{2,2}=='y' set(handles.checkbox22,'Value',1); end
set(handles.checkbox23,'Value',0);
if R.FieldsAnc{3,2}=='y' set(handles.checkbox23,'Value',1); end

%Pass={'Descending','n';'Ascending','n'};
set(handles.radiobuttonPassDescending,'Value',0);
if R.Pass{1,2}=='y' set(handles.radiobuttonPassDescending,'Value',1); end
set(handles.radiobuttonPassAscending,'Value',0);
if R.Pass{2,2}=='y' set(handles.radiobuttonPassAscending,'Value',1); end

S=get(handles.popupmenuTemporalAverages,'String');
j=find(strcmp(S,R.TAVG));
set(handles.popupmenuTemporalAverages,'Value',j);
S=get(handles.popupmenuTemporalAverages,'String');j=get(handles.popupmenuTemporalAverages,'Value');S=deblank(S{j});
TXT=['Temporal resolution is ' S];
set(handles.textTemporalResolution,'String',TXT,'ForegroundColor',[0 0 1])

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
if get(handles.radiobuttonToFiles,'Value') == 0
    set(handles.editDirectoryName,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
    set(handles.editFileNamePrefix,'BackgroundColor',[0.8 0.8 0.8],'Enable','Off')
elseif get(handles.radiobuttonToFiles,'Value') == 1
    set(handles.editDirectoryName,'BackgroundColor',[1. 1. 1.],'Enable','On')
    set(handles.editFileNamePrefix,'BackgroundColor',[1. 1. 1.],'Enable','On')
end

function CheckAvailableTimeRange(hObject, eventdata, handles)
S=get(handles.popupmenuTemporalAverages,'String');j=get(handles.popupmenuTemporalAverages,'Value');
TAVG=deblank(S{j});
eval(['DSI=handles.DSI.DSI_AIRS_' TAVG ';'])
Time1=DSI.Range.Time1;Time2=DSI.Range.Time2;
DATE1=datestr(Time1,'yyyy-mm-dd');    
DATE2=datestr(Time2,'yyyy-mm-dd');
TIMERANGE=['Available time range: [' DATE1 ' to ' DATE2 ']'];
set(handles.textTimeRange,'String',TIMERANGE); 

