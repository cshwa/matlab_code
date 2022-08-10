function varargout = gui_ocean_toolbox(varargin)
% The main figure m-file for the Matlab OPeNDAP Ocean Toolbox.
%
% gui_ocean_toolbox
% h = gui_ocean_toolbox
%
% INPUT(S):
%
%   [none] -- 
%
% OUTPUT(S):
%
%   h -- (optional) handle to GUI figure.
%
% OPeNDAP Science Team
% Copyright 2007, 2008, 2009
% $Version 3.1.2$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2007/05/01 0.0.x created, ceb
% 2008/01/01 0.0.x removed dataset buttons, replaced with radiobuttons;
%            cleaned up code, ceb
% 2008/03/01 2.0.0 released
% 2008/06/15 2.0.1 adding more guis
%
% Meri Sheremet
% 2009/09/04 3.1.2 update time ranges from dsi_datasetname.mat files
%==========================================================================


% Last Modified by GUIDE v2.5 21-Sep-2009 18:54:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_ocean_toolbox_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_ocean_toolbox_OutputFcn, ...
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


% --- Executes just before gui_ocean_toolbox is made visible.
function gui_ocean_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_ocean_toolbox (see VARARGIN)

load dsi_ocean_toolbox.mat
handles.DSI=DSI;

% Set interface name.
filename = 'ocean_toolbox.m';
txt = ['NASA/REASoN Ocean Data Portal - '];
tmp1 = ['Matlab OPeNDAP Ocean Toolbox'];
tmp2 = '';
txt = [txt,tmp1];
%if (~isdeployed)
%    tmp2 = strrep(op_version('mfile',filename),'Version ','');
%    txt = [txt,' (Ver ',tmp2,')'];
%end
set(hObject,'Name',txt);

% Set logos.
LogoSet(hObject, eventdata, handles)
% Set ToolTipString.
ToolTipSet(hObject, eventdata, handles)

for k=1:11
NUM=num2str(k);
eval(['set(handles.pushbutton' NUM ',''Visible'',''off'',''Enable'',''off'')'])
end
%set(handles.pushbutton11,'Visible','on','Enable','off','String','OceanWind')


DataSets=DSI.DataSetsTable(:,1);
for k=1:length(DataSets)
NUM=num2str(DSI.DataSetsTable{k,2});
S=DSI.DataSetsTable{k,1};
%eval(['set(handles.pushbutton' NUM ',''String'',''' S ''',''Enable'',''off'')'])
eval(['set(handles.pushbutton' NUM ',''String'',''' S ''',''Visible'',''on'',''Enable'',''on'')'])
NUM=num2str(DSI.DataSetsTable{k,2});
S=DSI.DataSetsTable{k,3};
eval(['op_tooltipwrap(handles.pushbutton' NUM ',S,50)'])
S=['Dataset Documentation'];
NUM=num2str(DSI.DataSetsTable{k,2}+100);
eval(['op_tooltipwrap(handles.pushbutton' NUM ',S,50)'])
end

Variables=DSI.VariablesTable(:,1);
for k=1:length(Variables)
NUM=num2str(DSI.VariablesTable{k,2});
S=DSI.VariablesTable{k,1};
eval(['set(handles.checkbox' NUM ',''String'',''' S ''',''Enable'',''on'',''Value'',1)'])
end
set(handles.pushbuttonSelectVars,'String','Clear All')

set(handles.checkbox2,'Visible','on','Enable','off','String','wind stress')

set(handles.checkboxTimeSeries,'Value',1)
set(handles.checkboxClimatology,'Value',1)
set(handles.checkboxSatellite,'Value',1)
set(handles.checkboxModel,'Value',1)

set(handles.text5,'String',...
['Available time range: [' datestr(DSI.Ocean_Toolbox.Range.Time1,'yyyy-mm-dd') ' to ' datestr(DSI.Ocean_Toolbox.Range.Time2,'yyyy-mm-dd') ']*'])

% Choose default command line output for gui_ocean_toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_ocean_toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function ToolTipSet(hObject, eventdata, handles)
% set tooltip for items other than datasets

% txt=['The OPeNDAP Toolbox uses a tool to combine data from multiple days ',...
%      'or from multiple requests into a convenient multidimensional structure. The ',...
%      'following options control how the merging tool places data from the users ',...
%      'workspace into this multidimensional structure.'];
% op_tooltipwrap(handles.pushbutton16,txt,50);

% txt=['The graphical interfaces and their datasets can ',...
%      'be filtered by specifying a set of latitudes, ',...
%      'longitudes, dates of interest or topic. ',...
%      'Use this tool to narrow your datasets of interest.'];
%op_tooltipwrap(handles.pushbutton14,txt,50);

txt=['year'];
op_tooltipwrap(handles.editTime1yyyy,txt,50);
op_tooltipwrap(handles.editTime2yyyy,txt,50);
txt=['month'];
op_tooltipwrap(handles.editTime1mm,txt,50);
op_tooltipwrap(handles.editTime2mm,txt,50);
txt=['day'];
op_tooltipwrap(handles.editTime1dd,txt,50);
op_tooltipwrap(handles.editTime2dd,txt,50);
txt=['starting latitude'];
op_tooltipwrap(handles.editLat1,txt,50);
txt=['ending latitude'];
op_tooltipwrap(handles.editLat2,txt,50);
txt=['starting longitude'];
op_tooltipwrap(handles.editLon1,txt,50);
txt=['ending longitude'];
op_tooltipwrap(handles.editLon2,txt,50);
txt=['Currently support only satellite and model gridded data.'];
op_tooltipwrap(handles.checkboxSatellite,txt,50);
op_tooltipwrap(handles.checkboxModel,txt,50);
txt=['Helpful documenation on how to use the software'];
op_tooltipwrap(handles.pushbutton27,txt,50);
txt=['support@opendap.org'];
op_tooltipwrap(handles.pushbutton28,txt,50);

% --- Outputs from this function are returned to the command line.
function varargout = gui_ocean_toolbox_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function pushbuttonOpendap_Callback(hObject, eventdata, handles)
web('http://www.opendap.org/');
%URL = 'http://www.opendap.org/';
%web(URL,'-browser')

function pushbuttonNASA_Callback(hObject, eventdata, handles)
web('http://www.nasa.gov');
%URL = 'http://www.nasa.gov';
%web(URL,'-browser')

function pushbutton1_Callback(hObject, eventdata, handles)
h = get_airs; 
InputX(hObject, eventdata, handles, h)
function pushbutton2_Callback(hObject, eventdata, handles)
h=get_aviso_altimetry;
InputX(hObject, eventdata, handles, h)
function pushbutton3_Callback(hObject, eventdata, handles)
h=get_ghrsst;
InputX(hObject, eventdata, handles, h)
function pushbutton4_Callback(hObject, eventdata, handles)
h=get_goes;
InputX(hObject, eventdata, handles, h)
function pushbutton5_Callback(hObject, eventdata, handles)
h=get_hycom;
InputX(hObject, eventdata, handles, h)
function pushbutton6_Callback(hObject, eventdata, handles)
h=get_modis;
InputX(hObject, eventdata, handles, h)
function pushbutton7_Callback(hObject, eventdata, handles)
h=get_oaflux;
InputX(hObject, eventdata, handles, h)
function pushbutton8_Callback(hObject, eventdata, handles)
h=get_oceancolor;
InputX(hObject, eventdata, handles, h)
function pushbutton9_Callback(hObject, eventdata, handles)
h=get_oceanwind;
InputX(hObject, eventdata, handles, h)
function pushbutton10_Callback(hObject, eventdata, handles)
h=get_pathfinder1km;
InputX(hObject, eventdata, handles, h)
function pushbutton11_Callback(hObject, eventdata, handles)
h=get_pathfinder4km;
InputX(hObject, eventdata, handles, h)

function pushbuttonReset_Callback(hObject, eventdata, handles)
DSI=handles.DSI;
% Enable all datasets
for k=1:length(DSI.DataSetsTable(:,1))
NUM=num2str(DSI.DataSetsTable{k,2});
eval(['set(handles.pushbutton' NUM ',''Enable'',''on'')'])
end

function pushbutton101_Callback(hObject, eventdata, handles)
k=1;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)
function pushbutton102_Callback(hObject, eventdata, handles)
k=2;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)
function pushbutton103_Callback(hObject, eventdata, handles)
k=3;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)
function pushbutton104_Callback(hObject, eventdata, handles)
k=4;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)
function pushbutton105_Callback(hObject, eventdata, handles)
k=5;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)
function pushbutton106_Callback(hObject, eventdata, handles)
k=6;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)
function pushbutton107_Callback(hObject, eventdata, handles)
k=7;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)
function pushbutton108_Callback(hObject, eventdata, handles)
k=8;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)
function pushbutton109_Callback(hObject, eventdata, handles)
k=9;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)
function pushbutton110_Callback(hObject, eventdata, handles)
k=10;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)
function pushbutton111_Callback(hObject, eventdata, handles)
k=11;
URL=handles.DSI.DataSetsTable{k,4};
web(URL)

function checkbox1_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox2_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox3_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox4_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox5_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox6_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox7_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox8_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox9_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox10_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox11_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox12_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox13_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox14_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox15_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox16_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox17_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkbox18_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)

function pushbuttonSelectVars_Callback(hObject, eventdata, handles)
S=get(handles.pushbuttonSelectVars,'String');
N=length(handles.DSI.VariablesTable(:,1));
if S(1) == 'S' % Select All
    for k=1:N
NUM=num2str(k);
eval(['set(handles.checkbox' NUM ',''Value'',1)'])
    end
set(handles.pushbuttonSelectVars,'String','Clear All')
elseif S(1) == 'C' % Clear All
    for k=1:N
NUM=num2str(k);
eval(['set(handles.checkbox' NUM ',''Value'',0)'])
    end
set(handles.pushbuttonSelectVars,'String','Select All')
end
CheckGUI(hObject, eventdata, handles)

function editTime1yyyy_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function editTime1yyyy_CreateFcn(hObject, eventdata, handles)
function editTime1mm_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function editTime1mm_CreateFcn(hObject, eventdata, handles)
function editTime1dd_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function editTime1dd_CreateFcn(hObject, eventdata, handles)

function editTime2yyyy_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function editTime2yyyy_CreateFcn(hObject, eventdata, handles)
function editTime2mm_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function editTime2mm_CreateFcn(hObject, eventdata, handles)
function editTime2dd_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function editTime2dd_CreateFcn(hObject, eventdata, handles)

function editLat1_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function editLat1_CreateFcn(hObject, eventdata, handles)
function editLat2_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function editLat2_CreateFcn(hObject, eventdata, handles)
function editLon1_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function editLon1_CreateFcn(hObject, eventdata, handles)
function editLon2_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function editLon2_CreateFcn(hObject, eventdata, handles)

function pushbutton27_Callback(hObject, eventdata, handles)
URL = 'http://www.oceanographicdata.org/Toolbox/documentation.html';
web(URL)

function pushbutton28_Callback(hObject, eventdata, handles)
URL = 'http://www.opendap.org/support/index.html';
web(URL)

function checkboxTimeSeries_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkboxClimatology_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkboxSatellite_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)
function checkboxModel_Callback(hObject, eventdata, handles)
CheckGUI(hObject, eventdata, handles)

function LogoSet(hObject, eventdata, handles)

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
CX=(NX+8)/ppcx; % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonOpendap,'CData',img)
pos = get(handles.pushbuttonOpendap,'position'); %units will be in characters
set(handles.pushbuttonOpendap,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_nasa1.jpg'); %('logo_nasamain.jpg');
[NY,NX,N3]=size(img);
CX=(NX+8)/ppcx; % button size in characters including 4 pixel grace margin
CY=(NY+12)/ppcy;
set(handles.pushbuttonNASA,'CData',img)
pos = get(handles.pushbuttonNASA,'position'); %units will be in characters
set(handles.pushbuttonNASA,'position',[pos(1),pos(2),CX,CY]);

[img,MAP]=imread('logo_questionmark.jpg');
[N1,N2,N3]=size(img);
set(handles.pushbutton101,'string','','CData',img,'Units','pixels')
set(handles.pushbutton102,'string','','CData',img,'Units','pixels')
set(handles.pushbutton103,'string','','CData',img,'Units','pixels')
set(handles.pushbutton104,'string','','CData',img,'Units','pixels')
set(handles.pushbutton105,'string','','CData',img,'Units','pixels')
set(handles.pushbutton106,'string','','CData',img,'Units','pixels')
set(handles.pushbutton107,'string','','CData',img,'Units','pixels')
set(handles.pushbutton108,'string','','CData',img,'Units','pixels')
set(handles.pushbutton109,'string','','CData',img,'Units','pixels')
set(handles.pushbutton110,'string','','CData',img,'Units','pixels')
set(handles.pushbutton111,'string','','CData',img,'Units','pixels')

pos = get(handles.pushbutton101,'position');
set(handles.pushbutton101,'position',[pos(1) pos(2) N1 N2])
pos = get(handles.pushbutton102,'position');
set(handles.pushbutton102,'position',[pos(1) pos(2) N1 N2])
pos = get(handles.pushbutton103,'position');
set(handles.pushbutton103,'position',[pos(1) pos(2) N1 N2])
pos = get(handles.pushbutton104,'position');
set(handles.pushbutton104,'position',[pos(1) pos(2) N1 N2])
pos = get(handles.pushbutton105,'position');
set(handles.pushbutton105,'position',[pos(1) pos(2) N1 N2])
pos = get(handles.pushbutton106,'position');
set(handles.pushbutton106,'position',[pos(1) pos(2) N1 N2])
pos = get(handles.pushbutton107,'position');
set(handles.pushbutton107,'position',[pos(1) pos(2) N1 N2])
pos = get(handles.pushbutton108,'position');
set(handles.pushbutton108,'position',[pos(1) pos(2) N1 N2])
pos = get(handles.pushbutton109,'position');
set(handles.pushbutton109,'position',[pos(1) pos(2) N1 N2])
pos = get(handles.pushbutton110,'position');
set(handles.pushbutton110,'position',[pos(1) pos(2) N1 N2])
pos = get(handles.pushbutton111,'position');
set(handles.pushbutton111,'position',[pos(1) pos(2) N1 N2])
%==========================================================================

function CheckGUI(hObject, eventdata, handles)
DSI=handles.DSI;

LCheckTime=1; 
% get strings from GUI menu and convert them to numbers
yyyy=deblank(get(handles.editTime1yyyy,'String'));
mm=deblank(get(handles.editTime1mm,'String')); if length(mm) == 1 mm=['0' mm]; end
dd=deblank(get(handles.editTime1dd,'String')); if length(dd) == 1 dd=['0' dd]; end
if isempty(yyyy)||isempty(mm)||isempty(dd)
    LCheckTime=0; % 
else
R.DATE1=[yyyy mm dd];
R.Time1=datenum(R.DATE1,'yyyymmdd');
end
yyyy=deblank(get(handles.editTime2yyyy,'String'));
mm=deblank(get(handles.editTime2mm,'String')); if length(mm) == 1 mm=['0' mm]; end
dd=deblank(get(handles.editTime2dd,'String')); if length(dd) == 1 dd=['0' dd]; end
if isempty(yyyy)||isempty(mm)||isempty(dd)
    LCheckTime=0; % 
else
R.DATE2=[yyyy mm dd];
R.Time2=datenum(R.DATE2,'yyyymmdd');
end

LCheckLat=1;
R.LAT1=str2num(deblank(get(handles.editLat1,'String')));
R.LAT2=str2num(deblank(get(handles.editLat2,'String')));
if isempty(R.LAT1)||isempty(R.LAT2) 
    LCheckLat=0;
end
LCheckLon=1;
R.LON1=str2num(deblank(get(handles.editLon1,'String')));
R.LON2=str2num(deblank(get(handles.editLon2,'String')));
if isempty(R.LON1)||isempty(R.LON2) 
    LCheckLon=0;
end

DataSets=DSI.DataSetsTable(:,1);
for k=1:length(DataSets)
NUM=num2str(DSI.DataSetsTable{k,2});
DSName=DataSets{k};

% Check Selected Variables
LVAR=0;
    for j=1:length(DSI.VariablesTable(:,1)); % raw in table
        NUMVAR=num2str(j);
        eval(['LCB=get(handles.checkbox' NUMVAR ',''Value'');'])
        if LCB==1         
            VarDataSets=DSI.VariablesTable{j,3};
            if ~isempty(strmatch(DSName,VarDataSets))
                LVAR=1;
            end
        end
    end

% Check Temporal
LTMP=0;
    for j=1:length(DSI.TemporalTable(:,1)); % raw in table
        CB=DSI.TemporalTable{j,2};
        eval(['LCB=get(handles.' CB ',''Value'');'])
        if LCB==1 & ~isempty(strmatch(DSName,DSI.TemporalTable{j,3}))
                LTMP=1;
        end
    end
% Check Source
LSRC=0;
    for j=1:length(DSI.SourceTable(:,1)); % raw in table
        CB=DSI.SourceTable{j,2};
        eval(['LCB=get(handles.' CB ',''Value'');'])
        if LCB==1 & ~isempty(strmatch(DSName,DSI.SourceTable{j,3}))
                LSRC=1;
        end
    end
% Check Time
LTIME=1;
    if LCheckTime==1
        eval(['Time1=DSI.' DSName '.Range.Time1;'])
        eval(['Time2=DSI.' DSName '.Range.Time2;'])
        if Time1>R.Time2 || Time2<R.Time1
            LTIME=0; %DataSet and Request do not overlap
        end
    end
% Check Lat
LLAT=1;
    if LCheckLat==1
        eval(['Lat1=DSI.' DSName '.Range.Latitude1;']);
        eval(['Lat2=DSI.' DSName '.Range.Latitude2;']);
        if Lat1>R.LAT2 || Lat2<R.LAT1
            LLAT=0; %DataSet and Request do not overlap
        end
    end
% Check Lon
LLON=1;
    if LCheckLon==1
        eval(['Lon1=DSI.' DSName '.Range.Longitude1;']);
        eval(['Lon2=DSI.' DSName '.Range.Longitude2;']);
        % because longitudes can be entered -180<180, 0<360, or even >360
        % formula is more complicated
        CD=(Lon1+Lon2)*0.5; % Center of DataSet
        RD=(Lon2-Lon1)*0.5; % Radius of DataSet
        CR=(R.LON1+R.LON2)*0.5; % Center of Request
        RR=(R.LON2-R.LON1)*0.5; % Radius of Request
        % if abs(CD-CR)>180 shift by 360 deg
        CR=CR+floor((CD+180-CR)/360)*360;
        if abs(CR-CD)>RR+RD
            LLON=0; %DataSet and Request do not overlap
        end
    end
    
L=LVAR*LTMP*LSRC*LTIME*LLAT*LLON;
    
    if L==1
    eval(['set(handles.pushbutton' NUM ',''Enable'',''on'')'])
    else
    eval(['set(handles.pushbutton' NUM ',''Enable'',''off'')'])
    end
    
set(handles.checkbox2,'Visible','on','Enable','off','String','wind stress')
end


%==========================================================================
function InputX(hObject, eventdata, handles, h)
hh=guidata(h);
S=get(handles.editTime1yyyy,'String');set(hh.editTime1yyyy,'string',S)
S=get(handles.editTime1mm  ,'String');set(hh.editTime1mm  ,'string',S)
S=get(handles.editTime1dd  ,'String');set(hh.editTime1dd  ,'string',S)
S=get(handles.editTime2yyyy,'String');set(hh.editTime2yyyy,'string',S)
S=get(handles.editTime2mm  ,'String');set(hh.editTime2mm  ,'string',S)
S=get(handles.editTime2dd  ,'String');set(hh.editTime2dd  ,'string',S)
S=get(handles.editLat1     ,'String');set(hh.editLat1     ,'string',S)
S=get(handles.editLat2     ,'String');set(hh.editLat2     ,'string',S)
S=get(handles.editLon1     ,'String');set(hh.editLon1     ,'string',S)
S=get(handles.editLon2     ,'String');set(hh.editLon2     ,'string',S)
%==========================================================================

