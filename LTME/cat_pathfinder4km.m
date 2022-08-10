function [Time,URLSITE,URLPATH,URLFILE]=cat_pathfinder4km(R)

%function [Time,URLSITE,URLPATH,URLFILE]=cat_pathfinder4km_clim(R)
if isfield(R,'DataSetBranch') & strcmp(R.DataSetBranch,'Pathfinder4km_Climatology')
% catalog of files for Pathfinder4km Climatology
% input R.TAVG,R.Passes, optional R.DATE1,R.DATE2,R.DATEINCR

%URLSITE='http://satdat1.gso.uri.edu/opendap/sea_surface_temperature/climatology/pathfinder/Version5.0/';
URLSITE='http://data.nodc.noaa.gov/opendap/pathfinder/Version5.0_Climatologies';
if ~isfield(R,'TAvg') R.TAvg='daily'; end
if ~isfield(R,'Passes') R.Passes={'night','day'}; end

% MATLAB 7.5 has a bug in datenum at year 0001-04-01
% use year 2000 as base instead to deal with climatology
% Note MATLAB datenum uses Julian calendar before October 4, 1582
% and Gregorian calendar after October 15, 1582
% 2000 (and 1600) is a leap year in both Julian and Gregorian calendars
% We need a leap year as a base year in order to access yearday 366
if ~isfield(R,'BaseYear') 
    R.BaseYear=0000; 
elseif ischar(R.BaseYear) 
    R.BaseYear=str2num(R.BaseYear);
end
YYYY=num2str(R.BaseYear,'%04d');
if ~isfield(R,'DATE1') R.DATE1='00000101'; end
if ~isfield(R,'DATE2') R.DATE2='00001231'; end
if ~isfield(R,'DATEINCR') R.DATEINCR=1; end
JD1=datenum([YYYY R.DATE1(5:end)],'yyyymmdd');
JD2=datenum([YYYY R.DATE2(5:end)],'yyyymmdd');
% to allow intervals spanning a New Year, e.g. Dec-Jan
if JD2<JD1 JD2=datenum([num2str(R.BaseYear+1,'%04d') R.DATE2(5:end)],'yyyymmdd'); end
   
PASSPATH='Combined/';PASSFN='combined';hhd=8/24; 
    if isempty(strmatch(R.Passes,'day')) 
%R.Clim_Passes='night';
PASSPATH='Night/';PASSFN='night';hhd=2/24;
    end
    if isempty(strmatch(R.Passes,'night')) 
%R.Clim_Passes='day';
PASSPATH='Day/';PASSFN='day';hhd=14/24;
    end

%error: R2007b
% datenum('00000101','yyyymmdd')=1
% datenum(0,1,1)=0
    
    if ~isempty(strmatch(R.TAvg,{'daily','day','1day'}))
URLPATH=['/Daily/' PASSPATH];ND=1;NT=366;TAVGFN='day';TAVGFNFMT='%03d';
k=JD1;j=k-datenum(str2num(datestr(k,'yyyy')),1,1)+1;JT1=floor((j-1)/ND)+1;
k=JD2;j=k-datenum(str2num(datestr(k,'yyyy')),1,1)+1;JT2=floor((j-1)/ND)+1; 
    elseif ~isempty(strmatch(R.TAvg,{'5day','pentad'}))
URLPATH=['/5day/' PASSPATH];ND=5;NT=73;TAVGFN='pentad';TAVGFNFMT='%02d';
k=JD1;j=k-datenum(str2num(datestr(k,'yyyy')),1,1)+1;JT1=floor((j-1)/ND)+1;
k=JD2;j=k-datenum(str2num(datestr(k,'yyyy')),1,1)+1;JT2=floor((j-1)/ND)+1; 
    elseif ~isempty(strmatch(R.TAvg,{'7day','week','weekly'}))
URLPATH=['/7day/' PASSPATH];ND=7;NT=52;TAVGFN='week';TAVGFNFMT='%02d';
k=JD1;j=k-datenum(str2num(datestr(k,'yyyy')),1,1)+1;JT1=floor((j-1)/ND)+1;
k=JD2;j=k-datenum(str2num(datestr(k,'yyyy')),1,1)+1;JT2=floor((j-1)/ND)+1; 
    elseif ~isempty(strmatch(R.TAvg,{'8day','octad'}))
URLPATH=['/8day/' PASSPATH];ND=8;NT=45;TAVGFN='octad';TAVGFNFMT='%02d';
k=JD1;j=k-datenum(str2num(datestr(k,'yyyy')),1,1)+1;JT1=floor((j-1)/ND)+1;
k=JD2;j=k-datenum(str2num(datestr(k,'yyyy')),1,1)+1;JT2=floor((j-1)/ND)+1; 
    elseif ~isempty(strmatch(R.TAvg,{'monthly','month','30day'}))
URLPATH=['/Monthly/' PASSPATH];ND=30;NT=12;TAVGFN='month';TAVGFNFMT='%02d';
k=JD1;JT1=str2num(datestr(k,'mm'));
k=JD2;JT2=str2num(datestr(k,'mm'));
    elseif ~isempty(strmatch(R.TAvg,{'seasonal','season','90day'}))
URLPATH=['/Seasonal/' PASSPATH];ND=90;NT=4;TAVGFN='season';TAVGFNFMT='%01d';
k=JD1;j=str2num(datestr(k,'mm'));JT1=floor((j-1)/3)+1;
k=JD2;j=str2num(datestr(k,'mm'));JT2=floor((j-1)/3)+1;
    elseif ~isempty(strmatch(R.TAvg,{'annual','yearly','year','365day','12month'}))
URLPATH=['/Annual/' PASSPATH];ND=366;NT=1;TAVGFN='annual';TAVGFNFMT='%00d';
JT1=1;
JT2=1;
    end

if JT2<JT1 JT2=JT2+NT; end
kF=0;
for k=JT1:R.DATEINCR:JT2
kT=mod(k-1,NT)+1;
kF=kF+1;
URLFILE{kF}=[TAVGFN num2str(kT,TAVGFNFMT) '_' PASSFN '.hdf'];
end
URLFILE=URLFILE';

k=(JT1:R.DATEINCR:JT2)';
[tA tB tC]=op_tavgind2datenum(k,R);
Time=tA+hhd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else

% catalog of files for Pathfinder4km
% based on the file naming conventions described in PFV50_UserGuide.pdf
% input R.TAVG,R.Passes, optional R.DATE1,R.DATE2,R.DATEINCR
%
% NODC
% opendap access:
% http://data.nodc.noaa.gov/cgi-bin/nph-dods/pathfinder
% climatology
% http://data.nodc.noaa.gov/cgi-bin/nph-dods/pathfinder/Version5.0_Climatologies

% JPL NASA
% http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/
%data_v5/5day/night/04km/1991/sst/contents.html
% Example of an individual data file URL:
%'http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder
%/data_v5/daily/ascending/04km/2005/sst/2005001.s04d4pfv50-sst.hdf
%?sst[0:1:4095][0:1:8191]'

% old
%URLSITE='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5';

URLSITE='http://data.nodc.noaa.gov/opendap/pathfinder';

%URLPass={'ascending';'descending'};
DSPass={'night';'day'};
DSPassLST={2/24;14/24}; % approx local solar time of a pass 
Fields={'sst';'bsst';'sdev';'num';'qual';'msk1';'msk2'};
Fields2={'sst';'bsst';'sdev';'num';'qual';'mask1';'mask2'};

if ~isfield(R,'TAvg') R.TAvg='daily'; end
TAVG=R.TAvg;
TAvg={'Yearly';'Monthly';'8day';'7day';'5day';'Daily'}; % directory name in the latest dataset
TAVG2=TAvg{strmatch(lower(TAVG),lower(TAvg))};           % directory name in the latest dataset

if isfield(R,'TimeRange')
% dataset available time range
URL1='http://data.nodc.noaa.gov/opendap/pathfinder/Version5.1/Daily/1981/';
YYYY1='1981';
ddd=htmldirlist(URL1,YYYY1);
ddd=min(ddd); %236
DSDate1=datenum(1981,01,ddd);
URL2='http://data.nodc.noaa.gov/opendap/pathfinder/Version5.0_interim_NOAA18/pfrt/Daily/2009/';
YYYY2='2009';
ddd=htmldirlist(URL2,YYYY2);
ddd=max(ddd);
DSDate2=datenum(2009,01,ddd);
else
R.DataSetBranch='Pathfinder4km_daily';
%S=['load dsi_pathfinder4km.mat DSI_' R.DataSetBranch];
eval(['load dsi_pathfinder4km.mat DSI_' R.DataSetBranch])
eval(['DSI=DSI_' R.DataSetBranch ';'])
eval(['clear DSI_' R.DataSetBranch ';'])
DSDate1=DSI.Range.Time1;
DSDate2=DSI.Range.Time2;
end
DSDATE1=datestr(DSDate1,'yyyymmdd');
DSDATE2=datestr(DSDate2,'yyyymmdd');

DSDateRT=datenum('20070101','yyyymmdd'); % real time data after DSDateRT


%timeline
% http://data.nodc.noaa.gov/opendap/pathfinder/
% 1981 236 -- 1985 030            Version5.1 (recommended through 1985 009)
% 1985 004 -- 2006 365            Version5.0
% 2007 001 -- 2008 365            Version5.0_interim
% 2009 001 -- 2009 ???            Version5.0_interim_NOAA18/pfrt/ (these
% are not quality controlled)
%Since NOAA-7 data is less sparse than NOAA-9, and the first 3 days of NOAA-9 are missing, I would say to use Version 5.1 for:
%-First 2 sets of 5-day files (1985001-1985005 and 1985006-1985010)
%-First set of 7-day files (198501)
%-First 8-day file (1985001-1985008)
% And for the monthly use Version 5.0.
% introduce TVer51end to control transitions between version 5.1 and 5.0 for various temporal averages 
TVerOffSet=0;
if strcmp(TAVG,'daily')    jpfv51end=9;end
if strcmp(TAVG,'5day')     jpfv51end=10;end
if strcmp(TAVG,'7day')     jpfv51end=7;end
if strcmp(TAVG,'8day')     jpfv51end=8;end
if strcmp(TAVG,'monthly')  jpfv51end=0;end
if strcmp(TAVG,'yearly')   jpfv51end=0;end


if ~isfield(R,'Passes') R.Passes=DSPass; end
if ~isfield(R,'DATE1') R.DATE1=DSDATE1; end
if ~isfield(R,'DATE2') R.DATE2=DSDATE2; end
if ~isfield(R,'DATEINCR') R.DATEINCR=1; end

Passes=R.Passes;DATE1=R.DATE1;DATE2=R.DATE2;DATEINCR=R.DATEINCR;
Date1=datenum(DATE1,'yyyymmdd');Date2=datenum(DATE2,'yyyymmdd');
%year1=max(str2num(DSDATE1(1:4)),str2num(DATE1(1:4)))
%year2=min(str2num(DSDATE2(1:4)),str2num(DATE2(1:4)))
year1=str2num(DSDATE1(1:4));
year2=str2num(DSDATE2(1:4));

SatPasses={};SatPassesLST={};
for kP=1:length(DSPass)
    if ~isempty(strmatch(DSPass{kP},Passes)) 
        SatPasses{end+1}=DSPass{kP};
        SatPassesLST{end+1}=DSPassLST{kP};
    end
end
% passes are indicated by 1,3
%except between these dates             
J1PASS10=datenum(2003,01,01); 
J2PASS10=datenum(2005,01,155);
if strcmp(TAVG,'daily')    J2PASS10=datenum(2005,01,155);end
if strcmp(TAVG,'5day')     J2PASS10=datenum(2005,01,155);end
if strcmp(TAVG,'7day')     J2PASS10=datenum(2005,01,154);end
if strcmp(TAVG,'8day')     J2PASS10=datenum(2005,01,160);end
if strcmp(TAVG,'monthly')  J2PASS10=datenum(2005,07,0);end
if strcmp(TAVG,'yearly')   J2PASS10=datenum(2005,01,155);end

% 4 2
% tricky for weekly averages
% yearly average for 2005 is not defined. file is not available.

RESO='04'; % 4km
for j=1:length(Fields)
FIELD=Fields{j};TYPE=FIELD;
BitCodes{j}='m'; % for nobs,qual,msk1,msk2
if ~isempty(strmatch(TYPE,{'sst','bsst','sdev'})) BitCodes{j}='s'; end
end


% initialize arrays
Time=[];
for j=1:length(Fields)
FIELD2=Fields2{j};
eval(['URLPATH.' FIELD2 '={};'])
eval(['URLFILE.' FIELD2 '={};'])
end

for year=year1:year2
YYYY=num2str(year);
    NDays=datenum(year,12,31)-datenum(year,01,01)+1; 
    if strcmp(TAVG,'daily')     NT=NDays;         AVGPERIOD='d';  end
    if strcmp(TAVG,'5day')      NT=floor(NDays/5); AVGPERIOD='5';  end
    if strcmp(TAVG,'7day')      NT=floor(NDays/7); AVGPERIOD='w';  end
    if strcmp(TAVG,'8day')      NT=floor(NDays/8); AVGPERIOD='8';  end
    if strcmp(TAVG,'monthly')   NT=12;            AVGPERIOD='m';  end
    if strcmp(TAVG,'yearly')    NT=1;             AVGPERIOD='y';  end
        
    for kt=1:DATEINCR:NT % browse time moments within each year
        
        if strcmp(TAVG,'daily')
j1=kt;            
OBSDATE=[YYYY num2str(j1,'%03d')];
dj=1;
%NOMDATE=datestr(datenum(year,01,j1),'yyyymmdd'); % middle of interval
%NOMDATE=datestr(datenum(year,01,j1),'yyyymmdd'); % start date
NomDate=datenum(year,01,j1); % start date
        elseif strcmp(TAVG,'5day')
dj=5;j1=(kt-1)*dj+1;j2=j1+dj-1;            
OBSDATE=[YYYY num2str(j1,'%03d') '-' YYYY num2str(j2,'%03d')]; 
%NOMDATE=datestr(datenum(year,01,j1+2),'yyyymmdd');
%NOMDATE=datestr(datenum(year,01,j1),'yyyymmdd'); % start date
NomDate=datenum(year,01,j1); % start date
        elseif strcmp(TAVG,'7day')
dj=7;j1=(kt-1)*dj+1;j2=j1+dj-1;            
OBSDATE=[YYYY num2str(kt,'%02d')]; 
%NOMDATE=datestr(datenum(year,01,j1+3),'yyyymmdd');
%NOMDATE=datestr(datenum(year,01,j1),'yyyymmdd'); % start date
NomDate=datenum(year,01,j1); % start date
        elseif strcmp(TAVG,'8day')
dj=8;j1=(kt-1)*dj+1;j2=j1+dj-1;            
OBSDATE=[YYYY num2str(j1,'%03d') '-' YYYY num2str(j2,'%03d')]; 
%NOMDATE=datestr(datenum(year,01,j1+4),'yyyymmdd');
%NOMDATE=datestr(datenum(year,01,j1),'yyyymmdd'); % start date
NomDate=datenum(year,01,j1); % start date
        elseif strcmp(TAVG,'monthly')
dj=30;j1=(kt-1)*dj+1;j2=j1+dj-1;            
OBSDATE=[YYYY num2str(kt,'%02d')];
%NOMDATE=datestr(datenum(year,kt,15),'yyyymmdd');
%NOMDATE=datestr(datenum(year,kt,01),'yyyymmdd'); % start date
NomDate=datenum(year,kt,01); % start date
        elseif strcmp(TAVG,'yearly')
dj=365;
OBSDATE=[YYYY]; 
%NOMDATE=datestr(datenum(year,07,01),'yyyymmdd');
%NOMDATE=datestr(datenum(year,01,01),'yyyymmdd'); % start date
NomDate=datenum(year,01,01); % start date
        end
        
        for kP=1:length(SatPasses)
            PASS=SatPasses{kP};
            PassLST=SatPassesLST{kP};

            
% passes indicated 1,2,3,4            
DAYNIGHT='1'; PassLST=2/24; if strcmp(PASS,'day') DAYNIGHT='3'; PassLST=14/24; end  
if (NomDate>=J1PASS10) & (NomDate<=J2PASS10)
% passes 4 2  (22 and 10)    
DAYNIGHT='4'; PassLST=22/24; if strcmp(PASS,'day') DAYNIGHT='2'; PassLST=10/24; end 
end

Time1=NomDate+PassLST;


%timeline
% http://data.nodc.noaa.gov/opendap/pathfinder/
% 1981 236 -- 1985 030            Version5.1 (recommended through 1985 009)
% 1985 004 -- 2006 365            Version5.0
% 2007 001 -- 2008 365            Version5.0_interim
% 2009 001 -- 2009 ???            Version5.0_interim_NOAA18/pfrt/ (these
% are not quality controlled)
%Since NOAA-7 data is less sparse than NOAA-9, and the first 3 days of NOAA-9 are missing, I would say to use Version 5.1 for:
%-First 2 sets of 5-day files (1985001-1985005 and 1985006-1985010)
%-First set of 7-day files (198501)
%-First 8-day file (1985001-1985008)
% And for the monthly use Version 5.0.

if Time1 <= datenum(1985,01,jpfv51end)
PATH1=['/Version5.1/' TAVG2 '/' YYYY '/'];
VERSION='pfv51';
elseif Time1 >= datenum(1985,01,jpfv51end+1) & Time1 <= datenum(2006,01,365)
PATH1=['/Version5.0/' TAVG2 '/' YYYY '/'];
VERSION='pfv50';
        if strcmp(TAVG,'yearly') & Time1 >= datenum(2005,01,01) & Time1 <= datenum(2005,01,365)
% Yearly averages are not defined due to switch between passes 2,4 and 1,3
% use interim data with caveat
% Version5.0_interim/Yearly/2005/NOAA17/2005.s04y2pfrt-sst.hdf
% Version5.0_interim/Yearly/2005/NOAA17/2005.s04y4pfrt-sst.hdf
        PATH1=['/Version5.0_interim/' TAVG2 '/' YYYY '/NOAA17/'];
        VERSION='pfrt'; 
        end
elseif Time1 >= datenum(2007,01,01) & Time1 <= datenum(2008,01,365)
PATH1=['/Version5.0_interim/' TAVG2 '/' YYYY '/NOAA18/'];
VERSION='pfrt'; 
elseif Time1 >= datenum(2009,01,01) 
PATH1=['/Version5.0_interim_NOAA18/pfrt/' TAVG2 '/' YYYY '/'];
VERSION='pfrt';
end


                if ... 
NomDate >= DSDate1 & NomDate < DSDate2 & NomDate >= Date1 & NomDate <= Date2
Time(end+1)=Time1;

                for j=1:length(Fields)
FIELD=Fields{j};
FIELD2=Fields2{j};
TYPE=FIELD;
BITCODE=BitCodes{j};
BITS=''; if BITCODE=='s' & year  < 2002 & strcmp(VERSION,'pfv50') BITS='-16b'; end

%PATH1=['/' TAVG '/' PASS '/04km/' YYYY '/' FIELD '/'];

FN1=[ OBSDATE '.' BITCODE RESO AVGPERIOD DAYNIGHT VERSION '-' TYPE BITS '.hdf'];

eval(['URLPATH.' FIELD2 '{end+1}=PATH1;'])
eval(['URLFILE.' FIELD2 '{end+1}=FN1;'])
                end
            end
        end %kP
    end
end

% cut and transpose cell arrays
Time=Time';
%TIME=datestr(Time,'yyyymmddHHMMSS');
for j=1:length(Fields)
FIELD2=Fields2{j};
eval(['URLPATH.' FIELD2 '=URLPATH.' FIELD2 ''';'])
eval(['URLFILE.' FIELD2 '=URLFILE.' FIELD2 ''';'])
end

end % if strcmp(R.DataSetBranch,'Pathfinder4km_Climatology')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddd=htmldirlist(URL,YYYY)
% gets the list days in the directory pointed to by URL
%URL='http://data.nodc.noaa.gov/opendap/pathfinder/Version5.0_interim_NOAA18/pfrt/Daily/2009/';
%YYYY='2009';
S=urlread(URL);

%match for 2009ddd. in filenames such as '2009121.s04d3pfrt-sst.hdf'
% after 2009 search for three digits followed by a dot
SS=regexp(S,[YYYY '\d{3,3}\.'],'match');

%cut the the leading and trailing strings 
L1=length(YYYY);
L2=length('.');
S=SS;
for k=1:length(S)
S1=SS{k};    
S{k}=S1(L1+1:end-L2);
end
S(:);

% return ddd (yearday)
ddd=str2num(cell2mat(S(:)));



