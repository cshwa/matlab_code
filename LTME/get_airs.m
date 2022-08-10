function varargout = get_airs(varargin)
% Gets Atmospheric Infrared Sounder (AIRS) data.
%
% get_airs
% get_airs(request)
%
% INPUT(S):
%
%   request -- (optional) structure input which contains
%              information necessary to obtain data. The GUI
%              merely provides a graphical method of generating
%              the "request" variable. However, this visual
%              interface may be by-passed by this method.
%
% OUTPUT(S):
%
%   [none]
%
% OPeNDAP Science Team
% Copyright 2007
% $Version 2.1.0$

%==========================================================================
% Meri Sheremet
% July 21 2009
%
% REVISION HISTORY:
% 2007/01/01 created, ms
% 2009/07/21 netcdf compatible datastructures, CF compliant, native mode is added
%==========================================================================

DSName='AIRS';
dsname=lower(DSName); %lower case

if ~exist('loaddap')
errordlg('loaddap not found. Make sure it is in your MATLAB path.'); return
end               

dsi_filename = ['dsi_' dsname '.mat'];
%if (~isdeployed)
pname = strrep(mfilename('fullpath'),mfilename,'');
dsi_filename = [pname,dsi_filename];
%end
dsi_flag = 0;
if exist(dsi_filename,'file')    
%     Get date of file.
    dum = dir(dsi_filename);    
    if floor(now)>datenum(dum.date) %if today's date is different than that of file
        dsi_flag = 1;        
    end    
else %does not exist
    dsi_flag = 1;    
end
if dsi_flag==1 % update dsi once a day (it takes a few seconds)   
%h=msgbox(['Acquiring metadata from ' DSName ' site, please wait ...']);
h=msgbox(['Acquiring ' DSName ' metadata, please wait ...']);
% next line makes OK button invisible
hh=get(h,'Children'); set(hh(2),'Visible','off'); pause(1);
dsi_airs; % dsi_xxx
close(h)
end


if nargin < 1 
h = gui_airs;
if nargout
    varargout{1} = h;
end
elseif nargin > 1
disp('Incorrect number of arguments.'); return  
else
request=varargin{1}; 
request.RequestDate=datestr(now,'yyyy-mm-dd HH:MM:SS');
R=request;

a=version;a=a(1:3);a=str2num(a);
% for earlier than MATLAB 2008b (Ver 7.7) use third party mexnc package
if strcmp(R.SaveMode,'netcdf') & ~exist('mexnc') & (a < 7.7)
errordlg('mexnc not found. Make sure it is in your MATLAB path.'); return
end    


%R=SetR(R);
%R.iTimeIncr=R.DATEINCR;
R.DataSetBranch=['AIRS_' R.TAVG];
% load dsi_xxx.mat DSI_xxx_04km_8day
eval(['load dsi_' dsname '.mat DSI_' R.DataSetBranch ';'])
eval(['DSI=DSI_' R.DataSetBranch ';'])
eval(['clear DSI_' R.DataSetBranch ';'])
%R.Variables={R.Fields{:},R.Coordinates{:}};
%assignin('base','DSI',DSI);assignin('base','R',R);


if isempty(R.LAT1) | isempty(R.LAT2) | isempty(R.LON1) | isempty(R.LON2)
    errordlg(['Latitude/Longitude values are invalid or empty.']);
    return
end

Latitude=DSI.Latitude;   %(89.5:-1:-89.5)';
Longitude=DSI.Longitude; %(-179.5:1:179.5)';
%TempPresLvls, in mbar 
TempPresLvls=DSI.TempPresLvls; %[1000.0;925.0;850.0;700.0;600.0;500.0;400.0;300.0;250.0;200.0;150.0;100.0;70.0;50.0;30.0;20.0;15.0;10.0;7.0;5.0;3.0;2.0;1.5;1.0];
%H2OPresLvls, in mbar
H2OPresLvls=DSI.H2OPresLvls; %TempPresLvls(1:12);
%IREmisFreqs Frequencies for emissivities 
% reported in the AIRS Level-3 product in cm-1
IREmisFreqs=DSI.IREmisFreqs; %[832; 961; 1203; 2616];
% MWEmisFreqs Frequencies for microwave emissivity products 
% reported in AIRS Level-3 in GHz.
MWEmisFreqs=DSI.MWEmisFreqs; %[23.8; 50.3; 89.0];

% adjust entered longitudes modulo 360 to match the dataset range 
LonM=0.5*(Longitude(1)+Longitude(end));
R.LON1=R.LON1+floor((LonM+180-R.LON1)/360)*360;
R.LON2=R.LON2+floor((LonM+180-R.LON2)/360)*360;

% Find indices within the specified range of lons, lats
iLAT2=max(find(Latitude >= R.LAT1))-1;
iLAT1=min(find(Latitude <= R.LAT2))-1;
iLON1=min(find(Longitude >= R.LON1))-1;
iLON2=max(find(Longitude <= R.LON2))-1;

if(isempty(iLAT1) | isempty(iLAT2)) | (iLAT2<iLAT1)
    errordlg(['Requested latitudes are out of range.']);
    return
end
if(isempty(iLON1) | isempty(iLON2)) | (iLON2<iLON1)
    errordlg(['Requested longitudes are out of range.']);
    return
end
if(isempty(R.LATINCR)|isempty(R.LONINCR)|isempty(R.DATEINCR)|R.LATINCR==0|R.LONINCR==0|R.DATEINCR==0) 
    errordlg(['Subsampling values are invalid.']);
    return
end
if(isempty(R.DATE1) | isempty(R.DATE2)) | (str2num(R.DATE2)<str2num(R.DATE1))
    errordlg(['Selected dates are invalid.']);
    return
end

if R.Pass{1,2}=='n' && R.Pass{2,2}=='n' 
    errordlg(['No satellite passes selected.']);
    return
end
if R.SaveWorkspace=='n' & R.SaveFiles=='n'
    errordlg(['No saving mode selected.']);
    return
end
if length(R.DIRNAME) > 0
    if exist(R.DIRNAME)==0
        eval(['mkdir ' R.DIRNAME])
    end    
end

CLAT=['[' num2str(iLAT1) ':' num2str(R.LATINCR) ':' num2str(iLAT2) ']']; % ConstraintLatitude
CLON=['[' num2str(iLON1) ':' num2str(R.LONINCR) ':' num2str(iLON2) ']']; % ConstraintLongitude


% get U structure; DODS_URL,year,month,day %,hour,minute,second
U='';
%U=Get_AIRS_URLs(R.DATE1,R.DATE2,R.DATEINCR,R.TAVG);
U=cat_airs(R);

if isempty(U)
errordlg(['No data within the requested date range.']);
    return 
end
NU=length(U.day);

% depending on the variables selected make the URL constraint
% also count 2D fields
C='';kVR=0;NF2D=0;NV=0;NF=0;
[L1,L2]=size(R.Fields);
    for kV=1:L1
    if R.Fields{kV,2}=='y'
    NV=0;    
    VN=R.Fields{kV,1};
%example, third dimension name TempPrsLvls, the constraint TEMPPRESLVLS 
    CLEV='';
    NF=1;
        if ~isempty(R.Fields{kV,4}) 
            CLEV=eval(['R.' upper(R.Fields{kV,4})]);
            NF=length(str2num(CLEV));
            if NF==0
    errordlg('Incorrect level constraints for 3D variables. Reissue get_airs command and correct the problem.')
            return
            end
        end
        for kP=1:2
        if R.Pass{kP,2}=='y'
            NV=NV+1;
            CPASS=R.Pass{kP,1};CPASS=CPASS(1);
            kVR=kVR+1; jV(kVR)=kV; VR{kVR}=[VN '_' CPASS];
            C=[C VN '_' CPASS CLEV CLAT CLON ','];
            if R.FieldsAnc{1,2}=='y'
            NV=NV+1;    
            kVR=kVR+1; jV(kVR)=kV; VR{kVR}=[VN '_' CPASS '_sdev'];
            C=[C VN '_' CPASS '_sdev' CLEV CLAT CLON ','];
            end
            if R.FieldsAnc{2,2}=='y'
            NV=NV+1;
            kVR=kVR+1; jV(kVR)=kV; VR{kVR}=[VN '_' CPASS '_ct'];
            C=[C VN '_' CPASS '_ct' CLEV CLAT CLON ','];
            end
        end
        end
    NF2D=NF2D+NV*NF;
    end
    end
    NVR=kVR; % number of requested variables
    
% LandSeaMask    
    if R.FieldsAnc{3,2}=='y'
    NF2D=NF2D+1;    
    C=[C 'LandSeaMask' CLAT CLON ','];
    end
C=C(1:end-1); %remove the trailing comma
C=['?' C];
C;
% number of 2D fields
NF2D;

Nv=0;
for k=1:length(R.Fields(:,1))
    if R.Fields{k,2}=='y' Nv=Nv+1; end % number of checked fields
end
if Nv==0
    errordlg(['No variables have been selected.']);
    return
end
% Estimate the amount of requested data
MB=round(NU*NF2D*(iLAT2-iLAT1)/R.LATINCR*(iLON2-iLON1)/R.LONINCR*8/1024/1024);
MBLIMIT=100;
if MB>MBLIMIT
ButtonName=questdlg(['The amount of data requested is ' num2str(MB) ' MB.'], ...
        'AIRS Requested Data Amount','Continue','Cancel','Cancel');
    if strcmp(ButtonName,'Cancel') return; 
    end
end

if strcmp(R.SaveMode,'merged')==0 
% find the latest saved opendap frame number: 
% kw0 in workspace; kf0 in files
[kw0,kf0]=op_lastframenum(R);
end

% Progress bar
hwaitbar=waitbar(0,'Acquiring ','Name','AIRS Download Progress');

readme='';
kf=0;
for k=1:NU
DATE=datestr(datenum(U.year(k),U.month(k),U.day(k)),'yyyymmdd');    
xprogr=(k-1)/NU; % progress fraction
waitbar( xprogr , hwaitbar, ['Acquiring ' DATE],'Name','AIRS Download Progress');

% depending on the variables selected make the URL constraint
        for kP=1:2
            if R.Pass{kP,2}=='y'
clear D % data structure containing a single frame
                kf=kf+1;
C='';kVR=0;NF2D=0;NV=0;NF=0;
[L1,L2]=size(R.Fields);
    for kV=1:L1
    if R.Fields{kV,2}=='y'
    NV=0;    
    VN=R.Fields{kV,1};
%example, third dimension name TempPrsLvls, the constraint TEMPPRESLVLS 
    CLEV='';
    NF=1;
        if ~isempty(R.Fields{kV,4}) 
            CLEV=eval(['R.' upper(R.Fields{kV,4})]);
            NF=length(str2num(CLEV));
            if NF==0
    errordlg('Incorrect level constraints for 3D variables. Reissue get_airs command and correct the problem.')
            return
            end
        end
        
        
            NV=NV+1;
            CPASS=R.Pass{kP,1};CPASS=CPASS(1);
            kVR=kVR+1; jV(kVR)=kV; VR{kVR}=[VN '_' CPASS]; VRO{kVR}=[VN];
            C=[C VN '_' CPASS CLEV CLAT CLON ','];
            if R.FieldsAnc{1,2}=='y'
            NV=NV+1;    
            kVR=kVR+1; jV(kVR)=kV; VR{kVR}=[VN '_' CPASS '_sdev']; VRO{kVR}=[VN '_sdev'];
            C=[C VN '_' CPASS '_sdev' CLEV CLAT CLON ','];
            end
            if R.FieldsAnc{2,2}=='y'
            NV=NV+1;
            kVR=kVR+1; jV(kVR)=kV; VR{kVR}=[VN '_' CPASS '_ct']; VRO{kVR}=[VN '_ct'];
            C=[C VN '_' CPASS '_ct' CLEV CLAT CLON ','];
            end
        
        
    NF2D=NF2D+NV*NF;
    end
    end
    NVR=kVR; % number of requested variables
    
% LandSeaMask    
    if R.FieldsAnc{3,2}=='y'
    NF2D=NF2D+1;    
    C=[C 'LandSeaMask' CLAT CLON ','];
    end
C=C(1:end-1); %remove the trailing comma
C=['?' C];
C;
% number of 2D fields
NF2D;

URL1=U.DODS_URL(k,:);
URL=[URL1 C];
disp(URL)
loaddap(['+v','-e'],URL);%dods_err = 0 means no error.
%dods_err = 1 means error.
%dods_err_msg [variable where error message is stored]
if dods_err
    disp(dods_err_msg)
    msgbox('Could not access AIRS site. Sorry, try again later.');
    close(hwaitbar);
    return
end
A=loaddap(['-A','-e'],URL);


% loaddap returns the following data structures 
% with the names depending on the Pass and Sensor: 
%         descending.Data_2520Fields.TotH2OVap_D_sdev
%          ascending.Data_2520Fields.TotH2OVap_A_sdev
% descending_MW_only.Data_2520Fields.TotH2OVap_MW_D_sdev
%  ascending_MW_only.Data_2520Fields.TotH2OVap_MW_A_sdev
% we need to eliminate that and combine all variables into a single structure        
    for kVR=1:NVR 
VN=VR{kVR};VNO=VRO{kVR}; 
PREFIX='';
if ~isempty(findstr(VN,'_D')) PREFIX='descending.Data_2520Fields.'; end
if ~isempty(findstr(VN,'_A')) PREFIX= 'ascending.Data_2520Fields.'; end
if ~isempty(findstr(VN,'_D')) & ~isempty(findstr(VN,'_MW')) PREFIX='descending_MW_only.Data_2520Fields.'; end
if ~isempty(findstr(VN,'_A')) & ~isempty(findstr(VN,'_MW')) PREFIX= 'ascending_MW_only.Data_2520Fields.'; end
eval(['f=' PREFIX VN ';']);

% replace missing values with NaN
eval(['MV=A.Global_Attributes.' VN '.ml__FillValue;'])
nn=find(f == MV);
f(nn)=NaN;

eval(['D.' VNO '=f;']);
    end
    if R.FieldsAnc{3,2}=='y'
    D.LandMask=location.Data_2520Fields.LandSeaMask;
    end
% Add latitude, longitude, URL and Attributes 
%iLAT1
%iLAT2
eval(['D.latitude =Latitude (' num2str(iLAT1+1) ':' num2str(R.LATINCR) ':' num2str(iLAT2+1) ');'])
eval(['D.longitude=Longitude(' num2str(iLON1+1) ':' num2str(R.LONINCR) ':' num2str(iLON2+1) ');'])

    for kVR=1:NVR 
%Fields={        'Temperature',   'n',   3,    'TempPresLvls'   ,{'latitude','longitude','TempPresLvls'};
        if R.Fields{jV(kVR),3}==3
            VN3D=R.Fields{jV(kVR),4}; % example: VN3D='TempPresLvls';
            VN3DInd=str2num(eval(['R.' upper(VN3D)]))+1;
 %           ['D.' VN3D '=' VN3D '(VN3DInd);']
       eval(['D.' VN3D '=' VN3D '(VN3DInd);'])  
        end
    end

    tpass=0; % approximate LST of the pass
    if kP==1 tpass= 1.5/24; end % D pass 1:30 am
    if kP==2 tpass=13.5/24; end % A pass 1:30 pm
time_dnum=datenum(U.year(k),U.month(k),U.day(k))+tpass;
time_Gregorian=datestr(time_dnum,'yyyy-mm-dd HH:MM:SS');
time_sec1970=86400*(time_dnum-datenum('1970-01-01 00:00:00','yyyy-mm-dd HH:MM:SS'));
D.time=time_sec1970;
D.user_friendly_time=time_Gregorian;

%D.metadata.sst_variables={'latitude','longitude'};
    for kVR=1:NVR     
    VN=VR{kVR};VNO=VRO{kVR};
    ListDim=R.Fields{jV(kVR),5};
        for j=1:length(ListDim)
            Dim=ListDim{j};
            eval(['SDim=length(squeeze(D.' Dim '));']);
            eval(['D.Dimensions.' Dim '=SDim;']);
        end
    end
    D.Dimensions.time=1;
    for kVR=1:NVR     
    VN=VR{kVR};VNO=VRO{kVR};
    eval(['D.Variables.' VNO '=R.Fields{jV(kVR),5};'])
    end
    
    if R.FieldsAnc{3,2}=='y'
    D.Variables.LandMask={'latitude' 'longitude'};
    end
    
    for kVR=1:NVR     
    VN=VR{kVR};VNO=VRO{kVR};
    ListDim=R.Fields{jV(kVR),5};
        for j=1:length(ListDim)
            Dim=ListDim{j};
            eval(['D.Variables.' Dim '=ListDim(j);']);
        end
    end

    for kVR=1:NVR     
    FIELD=VRO{kVR};
    s=eval(['DSI.Attributes.' FIELD '.long_name;']);
    eval(['D.Attributes.' FIELD '.long_name=s;']);
    s=eval(['DSI.Attributes.' FIELD '.units;']);
    eval(['D.Attributes.' FIELD '.units=s;']);
%    eval(['D.Attributes.' FIELD '.url=URL;'])
%    eval(['D.Attributes.' FIELD '.OriginalAttributes=A;'])
   end
    
    if R.FieldsAnc{3,2}=='y'
    D.Attributes.LandMask.long_name=DSI.Attributes.LandMask.long_name;
    D.Attributes.LandMask.units=DSI.Attributes.LandMask.units;
    end

  
    for kVR=1:NVR     
    ListDim=R.Fields{jV(kVR),5};
        for j=1:length(ListDim)
            COORD=ListDim{j};
            s=eval(['DSI.Attributes.' COORD '.long_name';]);
            eval(['D.Attributes.' COORD '.long_name=s;']);
            s=eval(['DSI.Attributes.' COORD '.units;']);
            eval(['D.Attributes.' COORD '.units=s;']);
        end
    end

D.Attributes.url=URL;
D.Attributes.OriginalAttributes=A;

       
%D.Attributes.url=URL;
%D.Attributes.OriginalAttributes=A;
D.Attributes.DataSetName='AIRS';
D.Attributes.DataSetBranch=['AIRS_' R.TAVG];
%D.metadata.temporal='time_series';
% readme_tmp = ['Scale: multiplied "des_stress_Liu_V" by ',...
% 	num2str(attr.des_stress_Liu_V.scale_factor)];
% readme = [readme,'/',readme_tmp];
%readme_tmp=['Comment: CloudFrcVis descending pass is missing in dataset - disabled. '];
%readme = [readme,'/',readme_tmp];
%D.metadata.reference_time=datestr([1970,01,01,00,00,00],'yyyy-mm-dd HH:MM:SS');
readme=['time in seconds since reference_time: 1970-01-01 00:00:00;'  '/' ...
        'bad values replaced with NaN'];
D.Attributes.Readme = readme;  
D.Attributes.RequestDate = R.RequestDate;
D.Attributes.Request=R;

status=op_saved(D,R,kf,kf0,kw0);
%%%%
            end % if kP
       end % kP

end
xprogr=1; % progress fraction
waitbar( xprogr , hwaitbar, ['Acquiring ' DATE],'Name','AIRS Download Progress');
close(hwaitbar)

disp('Data have been acquired.')
end

%**************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U=Get_AIRS_URLs(DATE1,DATE2,DATEINCR,TAVG)

%URLCatBase='http://dods.gso.uri.edu/dods-3.4/nph-dods/catalog/';
URLCatBase='http://acdisc.sci.gsfc.nasa.gov/opendap/catalog/DatapoolCatalog/AIRS/';
URL2='';
if strcmp(TAVG,'daily') URL2='AIRX3STD_005-cat.dat'; UNAME='AIRX3STD_005'; end 
if strcmp(TAVG,'weekly') URL2='AIRX3ST8_005-cat.dat'; UNAME='AIRX3ST8_005'; end 
if strcmp(TAVG,'monthly') URL2='AIRX3STM_005-cat.dat'; UNAME='AIRX3STM_005'; end 

%% constrain date
CDATE=['&date(''' DATE1(1:4) '/' DATE1(5:6) '/' DATE1(7:8) ''',''' DATE2(1:4) '/' DATE2(5:6) '/' DATE2(7:8) ''')'];   
%
%URL=[URLCatBase URL2];
URL=[URLCatBase URL2 '?DODS_URL,year,month,day' CDATE];
loaddap(['+v','-e'],URL);
%dods_err = 0 means no error.
%dods_err = 1 means error.
%dods_err_msg [variable where error message is stored]
if dods_err
    disp(dods_err_msg)
    msgbox('Could not access Catalog server. Sorry, try again later.');
    close(hwaitbar);
    return
end

eval(['U=' UNAME ';'])
% returns the following structure with the list of data file URLs and date-time
% information
%
%AIRSX3STD = 
%
%    DODS_URL: [408x154 char]
%      second: [408x1 double]
%      minute: [408x1 double]
%        hour: [408x1 double]
%         day: [408x1 double]
%       month: [408x1 double]
%        year: [408x1 double]
if exist('U') & length(U)>0

%subsample explicitly    
%JD1=datenum(DATE1,'yyyymmdd');    
%JD2=datenum(DATE2,'yyyymmdd')+1;
%JD=datenum(U.year,U.month,U.day);
%j=find(JD>=JD1 & JD<JD2);
% subsample in time as requested by DATEINCR    
%j=j(1:DATEINCR:length(j));
%U.DODS_URL=U.DODS_URL(j,:);
%U.second  =U.second(j);
%U.minute  =U.minute(j);
%U.hour    =U.hour(j);
%U.day     =U.day(j);
%U.month   =U.month(j);
%U.year    =U.year(j);
%U;

else
U='';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
