function varargout = get_ghrsst(varargin)
% Gets GHRSST data series.
%
% Gets xxx data.
%
%  Usage:
%
%  get_xxx  
%
%  without arguments brings up GUI that allows a user to select variables, 
%  time, latitude, longitude, layer limits and save data to workspace or files. 
%  All selections are saved in a cell array Request.
%
%  get_xxx(Request)
%
%  gets data selected by a prior call which was saved in the cell array Request.
%
%  Files: get_xxx.m      - main function, this file
%         dsi_xxx.m      - describes/generates dataset inventory and interface
%         dsi_xxx.mat    - contains dataset inventory
%         gui_xxx.m      - GUI function
%         gui_xxx.fig    - GUI graphic window
%         gui_xxx_request.mat    - contains the last Request
%  
% OPeNDAP Science Team
% Copyright 2007-2009
% $Version 3.1.2$

%==========================================================================
% Meri Sheremet
%
% REVISION HISTORY:
% 2007/01/01 created, ms
% 2008/03/13 introduced DSI; 
% 2008/12 all toolbox rewritten to be uniform in structure, vs
% 2009/10/03 metastrings introduced to specify URL patterns
%==========================================================================

DSName='GHRSST';
dsname=lower(DSName); %lower case

if ~exist('loaddap')
errordlg('loaddap not found. Make sure it is in your MATLAB path.'); return
end               

dsi_filename = ['dsi_' dsname '.mat'];
%if (~isdeployed)
    pname = mfilename('fullpath');
    pname = pname(1:length(pname) - length(mfilename));
    dsi_filename = [pname,dsi_filename];
%end
dsi_flag = 0;
if exist(dsi_filename,'file')    
  %   Get date of file.
    dum = dir(dsi_filename);    
    if floor(now)>datenum(dum.date) %if today's date is different than that of file
        dsi_flag = 1;        
    end    
else %does not exist
    dsi_flag = 1;    
end
if dsi_flag==1    %
%h=msgbox(['Acquiring metadata from ' DSName ' site, please wait ...']);
h=msgbox(['Acquiring ' DSName ' metadata, please wait ...']);
% next line makes OK button invisible
hh=get(h,'Children'); set(hh(2),'Visible','off'); pause(1);
%dsi_ghrsst; % dsi_xxx
f=ftp('po.gso.uri.edu');cd(f,'/pub/downloads/oceantoolbox');binary(f);mget(f,['dsi_' dsname '.mat']),close(f);
close(h)
end

if nargin<1
% no arguments passed, get parameters from GUI

h = gui_ghrsst;
    if nargout>0
        varargout{1} = h; %handle is output.
    end
elseif nargin > 1
disp('Incorrect number of arguments.'); return  
else
% one argument Request, get xxx data  

request=varargin{1}; 
assignin('base','request',request);
request.RequestDate=datestr(now,'yyyy-mm-dd HH:MM:SS');
R=request;

a=version;a=a(1:3);a=str2num(a);
% for earlier than MATLAB 2008b (Ver 7.7) use third party mexnc package
if strcmp(R.SaveMode,'netcdf') & ~exist('mexnc') & (a < 7.7)
errordlg('mexnc not found. Make sure it is in your MATLAB path.'); return
end   

%R=SetR(R);
R.iTimeIncr=R.DATEINCR;
% load dsi_xxx.mat DSI_xxx_04km_8day
eval(['load dsi_' dsname '.mat DSI_' R.DataSetBranch ';'])
eval(['DSI=DSI_' R.DataSetBranch ';'])
eval(['clear DSI_' R.DataSetBranch ';'])
R.Variables={R.Fields{:},R.Coordinates{:}};
assignin('base','DSI',DSI);

status=op_checkrequest(DSI,R);           if status>0 return; end
D0=op_initd0(DSI,R); % initialize standard opendap structure;
[R,D0,status]=op_coord2ind(DSI,R,D0); if status>0 return; end
status=op_chkamntreqdata(DSI,R); if status>0 return; end
% find the latest saved opendap or netcdf frame number
[kw0,kf0]=op_lastframenum(R);
kf=0;

hwaitbar=waitbar(0,'Acquiring ','Name',[DSName ' Download Progress']);
% browse through time
iTime=R.iTime1:R.iTimeIncr:R.iTime2;
if(isempty(iTime)) 
    errordlg(['No data within specified range of dates/passes.']);
    return
end

%NT=length(R.iTime1:R.iTimeIncr:R.iTime2);
NT=length(iTime);
NV=length(R.Fields);
%for kT=R.iTime1:R.iTimeIncr:R.iTime2
for kT=iTime
kf=kf+1; % count time frames for naming opendap_nnnn and files    
DATE=datestr(DSI.Time(kT),'yyyy-mm-dd HH:MM:SS'); 
Time=DSI.Time(kT);
%TIME=DSI.TIME{kT};
YYYY=datestr(Time,'yyyy');
yearday=(Time-datenum(YYYY,'yyyy')+1);
DDD=sprintf('%03d',floor(yearday));
YYYYMMDD=datestr(Time,'yyyymmdd');
%YYYYMMDD_HHMM=[YYYYMMDD '_' TIME(12:13) TIME(15:16)]; % hour may equal 25

% processing each time frame and constructing data structure D
D=D0;
    
%D.time=DSI.Time(kT); % datenum;  
D.time=(DSI.Time(kT)-datenum(1970,1,1,0,0,0))*86400;% in seconds since 1970-01-01 00:00:00
D.Attributes.time.units='seconds since 1970-01-01 00:00:00';
D.user_friendly_time=datestr(DSI.Time(kT),'yyyy-mm-dd HH:MM:SS');

    %browse through checked Fields
    for kV=1:NV
FIELD=R.Fields{kV};        
xprogress=((kV-1)/NV+(kf-1))/NT; % progress fraction
%waitbar(xprogress, hwaitbar, ['Acquiring ' strrep(FIELD,'_',' ') ' ' DATE],'Name',[DSName ' Download Progress']);
waitbar(xprogress, hwaitbar, ['Acquiring ' ' ' DATE ' ' strrep(FIELD,'_',' ')],'Name',[DSName ' Download Progress']);
% a general request employed by get_DataSet program has the following pattern
% loaddap([URLSITE URLPATH URLFILE '?' URLCVAR URLCTIME URLCDEPTH URLCLAT URLCLON ])
%   constraints URLCTIME URLCDEPTH URLCLAT URLCLON  can be in different order
%   constraints URLCLAT URLCLON URLCDEPTH will be generated by program get

j=0;
for kk=1:length(DSI.Fields)
    if strcmp(DSI.Fields{kk},FIELD) 
        j=kk; 
    end
end
%URL=op_smarturl(URL,FIELD,kT)
URLSITE=DSI.URLSITE; if strfind(URLSITE,'URLSITE=') eval(URLSITE); end 
URLPATH=DSI.URLPATH; if strfind(URLPATH,'URLPATH=') eval(URLPATH); end 
URLFILE=DSI.URLFILE; if strfind(URLFILE,'URLFILE=') eval(URLFILE); end 

URLCVAR=DSI.Fields_NameDS{j};
CTIME=op_smarturl(DSI.URLCTIME,FIELD,kT);
FIELDR=DSI.Fields_NameReturned{j};
URL=[URLSITE URLPATH URLFILE '?' URLCVAR CTIME R.CLAT R.CLON];

disp(URL)
loaddap(['+v','-e'],URL);
%dods_err = 0 means no error.
%dods_err = 1 means error.
%dods_err_msg [variable where error message is stored]
            if dods_err
                Vars=eval(['D.Variables.' FIELD]);
                nd=[];
                for j=1:length(Vars)
                    nd(j)=length(eval(['D.' Vars{j}]));
                end
                FNAN=NaN(nd);    
        disp(dods_err_msg)
        disp(['Missing or corrupt file or site is not accessible. ' FIELD ' is set to NaN.'])
        eval(['D.' FIELD '=FNAN;']);
            else
A=loaddap(['+v','-A'],URL);

% Rename and change dimensions order
% current field is called f
eval(['f=squeeze(' FIELDR ');']);
%FREN=eval(['DSI.FormulaRen.' FIELD]);eval(FREN);
% replace missing values with NaN according to formula stored in DSI.FormulaNaN
eval(eval(['DSI.FormulaNaN.' FIELD]))
%eval(eval(['DSI.FormulaNaN'])); % same for all fields
% convert to physical values according to formula stored in DSI.FormulaCnv
eval(eval(['DSI.FormulaCnv.' FIELD]))
%eval(eval(['DSI.FormulaCnv'])); % same for all fields
% transpose matrix to have f(lat,lon)
%eval(eval(['DSI.FormulaOrderCoord']));


eval(['D.' FIELD '=f;']);
eval(['D.Attributes.' FIELD '.url=URL;']);
eval(['D.Attributes.' FIELD '.OriginalAttributes=A;']);
            end
    end %kV
status=op_saved(D,R,kf,kf0,kw0);     
end %kT

waitbar(1,hwaitbar,'Acquiring ','Name',[DSName ' Download Progress']);
close(hwaitbar)
%msgbox('The data have been acquired.')
disp('The data have been acquired.')

end % nargin

