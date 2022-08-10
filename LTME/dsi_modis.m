function dsi_modis(varargin)
% updates the DataSetInventory structures
% without arguments for all the following DataSetBranches
if nargin < 1 
DSI.DataSetBranch='MODIS_aqua_04km_daily_all',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_daily_day',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_daily_night',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_8day_all',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_8day_day',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_8day_night',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_monthly_all',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_monthly_day',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_monthly_night',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_annual_all', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_annual_day', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_04km_annual_night', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_daily_all',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_daily_day',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_daily_night',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_8day_all',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_8day_day',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_8day_night',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_monthly_all',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_monthly_day',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_monthly_night',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_annual_all', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_annual_day', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_aqua_09km_annual_night', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_daily_all',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_daily_day',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_daily_night',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_8day_all',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_8day_day',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_8day_night',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_monthly_all',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_monthly_day',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_monthly_night',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_annual_all', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_annual_day', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_04km_annual_night', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_daily_all',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_daily_day',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_daily_night',  DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_8day_all',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_8day_day',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_8day_night',   DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_monthly_all',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_monthly_day',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_monthly_night',DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_annual_all', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_annual_day', DSI=dsi_modis_common(DSI);    
DSI.DataSetBranch='MODIS_terra_09km_annual_night', DSI=dsi_modis_common(DSI);    
elseif nargin > 1
disp('Incorrect number of arguments.'); return  
else
% with one argument for the specified DataSetBranch only     
DSI.DataSetBranch=varargin{1}, DSI=dsi_modis_common(DSI); 
end

function DSI=dsi_modis_common(DSI)
S=DSI.DataSetBranch;
clear DSI
DSI.DataSetBranch=S;
u=findstr(S,'_');
DSI.DataSetName=S(1:u(1)-1); %'MODIS'
DSI.DataSetParameters={'Platform','Algorithm','HRes','TAvg','Pass'};
DSI.Platform=S(u(1)+1:u(2)-1); %'aqua';
DSI.HRes=S(u(2)+1:u(3)-1); %'04km';
    if strcmp(DSI.HRes,'04km')
     DSI.NY=4320;DSI.NX=8640;
    elseif strcmp(DSI.HRes,'09km')
     DSI.NY=2160;DSI.NX=4320;
    else
     disp('incorrect HRes')
     DSI.NY=1;DSI.NX=1;
    end     
DSI.TAvg=S(u(3)+1:u(4)-1); %'8day';
DSI.Pass=S(u(4)+1:end); %'night';
DSI=dsi_modis_common1(DSI);
    if strcmp(DSI.TAvg,'daily')
DSI=dsi_modis_files_daily_rt(DSI);
    else
DSI=dsi_modis_files(DSI);
    end    
DSI=dsi_modis_common2(DSI);

%function DSI=dsi_unisys_variables(DSI)
function DSI=dsi_modis_common1(DSI)
% describe DataSet structure patterns 
% that depend on variable names only 
% and not on time
% examples of data access:
%http
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/modis/
%data/aqua/L3_mapped/sst/daily/04km/2006/365/A2006365.L3m_DAY_SST_4.bz2.html
% aqua,terra;daily,8day,monthly,annual;04km,09km;2002,...,2008;001,...,365;

%loaddap
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/modis/
%data/aqua/L3_mapped/sst/daily/04km/2006/365/
%A2006365.L3m_DAY_SST_4.bz2?l3m_data[0:1:4319][0:1:8639],l3m_qual[0:1:4319][0:1:8639]

% catalog server: http
%http://dods.jpl.nasa.gov/catalogs/modis/L3/daily/04km/daily-04km-SST.dat
%or loaddap
%http://dods.jpl.nasa.gov/opendap/dods_catalogs/modis/L3/daily/04km/daily-04km-SST.dat
% only this catalog
% catalog does not work

FieldsTable={
% opendap name,dataset name,returned name,  long_name               ,gui_name              ,units,    
'sst'       ,'l3m_data','l3m_data','sea_surface_temperature'        ,'sst'                 ,'degrees C';
'sst_qual'  ,'l3m_qual','l3m_qual','sea_surface_temperature_quality','sst_quality'         ,'';
'sst4'      ,'l3m_data','l3m_data','sea_surface_temperature_4micron  ','sst4'              ,'degrees C';
'sst4_qual' ,'l3m_qual','l3m_qual','sea_surface_temperature_4micron_quality','sst4_quality','';
};
CoordinatesTable={
'latitude','latitude','degrees_north';
'longitude','longitude','degrees_east';
'time','time','days since 0000-01-01 00:00:00';
};
DSI.Fields=FieldsTable(:,1);
DSI.Coordinates=CoordinatesTable(:,1);

LaLo={'latitude','longitude'};
DSI.Dimensions={};
DSI.Variables.sst=LaLo;
DSI.Variables.sst_qual=LaLo;
DSI.Variables.sst4=LaLo;
DSI.Variables.sst4_qual=LaLo;

DSI.Variables.latitude={'latitude'};
DSI.Variables.longitude={'longitude'};
DSI.Variables.time={'time'};

for k=1:length(DSI.Fields)
    FIELD=DSI.Fields{k};
eval(['DSI.Attributes.' FIELD '.long_name=FieldsTable{k,4};']);
eval(['DSI.Attributes.' FIELD '.units    =FieldsTable{k,6};']);
end
for k=1:length(DSI.Coordinates)
    COORD=DSI.Coordinates{k};
eval(['DSI.Attributes.' COORD '.long_name=CoordinatesTable{k,2};']);
eval(['DSI.Attributes.' COORD '.units    =CoordinatesTable{k,3};']);
end

% Rename returned variables and change dimensions order
DSI.FormulaRen.sst=      'f=squeeze(l3m_data);';
DSI.FormulaRen.sst_qual= 'f=squeeze(l3m_qual);';
DSI.FormulaRen.sst4 =    'f=squeeze(l3m_data);';
DSI.FormulaRen.sst4_qual='f=squeeze(l3m_qual);';

% Replace missing values with NaN
DSI.FormulaNaN.sst='f(find(f==65535))=NaN;';
DSI.FormulaNaN.sst_qual='f(find(f==255))=NaN;';
DSI.FormulaNaN.sst4 ='f(find(f==65535))=NaN;';
DSI.FormulaNaN.sst4_qual='f(find(f==255))=NaN;';

% Convert to physical units
DSI.FormulaCnv.sst ='f=A.Global_Attributes.HDF_GLOBAL.Slope*f + A.Global_Attributes.HDF_GLOBAL.Intercept;';
DSI.FormulaCnv.sst_qual ='';
DSI.FormulaCnv.sst4='f=A.Global_Attributes.HDF_GLOBAL.Slope*f + A.Global_Attributes.HDF_GLOBAL.Intercept;';
DSI.FormulaCnv.sst4_qual ='';
% url location of the dataset files
DSI.URLSITE='http://dods.jpl.nasa.gov';
for k=1:length(DSI.Fields);
    FIELD=DSI.Fields{k};
eval(['DSI.URLPATH.' FIELD '={};'])
eval(['DSI.URLFILE.' FIELD '={};'])
eval(['DSI.URLCVAR.' FIELD '={};'])
% define empty arrays
eval(['DSI.Time.'  FIELD '=[];'])    % double array: center of time interval
eval(['DSI.TimeA.' FIELD '=[];'])    % double array: beginning of time interval
eval(['DSI.TimeB.' FIELD '=[];'])    % double array: end of time interval
end
DSI.URLCVAR.sst='l3m_data';
DSI.URLCVAR.sst_qual='l3m_qual';
DSI.URLCVAR.sst4='l3m_data';
DSI.URLCVAR.sst4_qual='l3m_qual';

DSI.FormulaCONSTR.sst=      'R.CLAT R.CLON';
DSI.FormulaCONSTR.sst_qual= 'R.CLAT R.CLON';
DSI.FormulaCONSTR.sst4=     'R.CLAT R.CLON';
DSI.FormulaCONSTR.sst4_qual='R.CLAT R.CLON';



function DSI=dsi_modis_files_daily(DSI)
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/modis/data/aqua/L3_mapped/sst/8day/04km/2004/A20040012004008.L3m_8D_NSST_4.bz2?l3m_data[0:1:4319][0:1:8639],l3m_qual[0:1:4319][0:1:8639]
Fields={'sst';'sst4'};
for kFields=1:length(Fields)
FIELD=Fields{kFields};
%for Year=2002:2003 %str2num(datestr(now,'yyyy'))
for Year=2002:str2num(datestr(now,'yyyy'))
    YYYY=num2str(Year);
    NDays=datenum(Year,12,31)-datenum(Year,01,00);
    for Day=1:NDays
    DDD=num2str(Day,'%03d');   
%URLPATH='sst/daily/04km/2006/365/';
    URLPATH=['/opendap/sea_surface_temperature/modis/data/' DSI.Platform '/L3_mapped/' FIELD '/' DSI.TAvg '/' DSI.HRes '/' YYYY '/' DDD '/'];
%A2006365.L3m_DAY_SST_4.bz2.html
% obtain list of files with extention 'bz2'
%FNList=htmldirlist(URL,EXT)
URL=[DSI.URLSITE URLPATH]
FNList=htmldirlist(URL,'bz2')
        for kt=1:length(FNList)
eval(['DSI.URLPATH.' FIELD '{end+1}=URLPATH;'])
eval(['DSI.URLFILE.' FIELD '(end+1)=FNList(kt);'])
%A2006365.L3m_DAY_SST_4.bz2      
FN=FNList{kt};
tA=datenum(str2num(FN(2:5)),01,str2num(FN(6:8)));
tB=tA+1;
    if ~isempty(strfind(FN,'_NSST')) % night images
t=tA;
    elseif ~isempty(strfind(FN,'_SST')) % day images add 12h to nominal time
t=tA+0.5;
    else % incorrect file name
disp(['incorrect file name ' FN])
tA=NaN;tB=NaN;t=NaN;
    end
    L=strcmp(DSI.Pass,'all') | (~isempty(strfind(FN,'_NSST'))&strcmp(DSI.Pass,'night')) | (~isempty(strfind(FN,'_SST'))&strcmp(DSI.Pass,'day'));
    if L
eval(['DSI.TimeA.' FIELD '(end+1)=tA;'])
eval(['DSI.TimeB.' FIELD '(end+1)=tB;'])
eval(['DSI.Time.'  FIELD '(end+1)=t;' ])
    end
        end %kt
    end %Day
end %Year
end %kFields

function DSI=dsi_modis_files(DSI)
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/modis/data/aqua/L3_mapped/sst/8day/04km/2004/A20040012004008.L3m_8D_NSST_4.bz2?l3m_data[0:1:4319][0:1:8639],l3m_qual[0:1:4319][0:1:8639]
Fields={'sst';'sst4'};
for kFields=1:length(Fields)
FIELD=Fields{kFields};
for Year=2002:str2num(datestr(now,'yyyy'))
%for Year=2002:2008 %str2num(datestr(now,'yyyy'))
    YYYY=num2str(Year);
    URLPATH=['/opendap/sea_surface_temperature/modis/data/' DSI.Platform '/L3_mapped/' FIELD '/' DSI.TAvg '/' DSI.HRes '/' YYYY '/'];
% obtain list of files with extention 'bz2'
%FNList=htmldirlist(URL,EXT)
URL=[DSI.URLSITE URLPATH]
FNList=htmldirlist(URL,'bz2')
        for kt=1:length(FNList)
eval(['DSI.URLPATH.' FIELD '{end+1}=URLPATH;'])
eval(['DSI.URLFILE.' FIELD '(end+1)=FNList(kt);'])
%A20040012004008.L3m_8D_NSST_4.bz2     
FN=FNList{kt};
tA=datenum(str2num(FN(2:5)),01,str2num(FN(6:8)));
tB=datenum(str2num(FN(9:12)),01,str2num(FN(13:15)));
    if ~isempty(strfind(FN,'_NSST')) % night images
t=ceil((tA+tB)*0.5);
    elseif ~isempty(strfind(FN,'_SST')) % day images add 12h to nominal time
t=ceil((tA+tB)*0.5)+0.5;
    else % incorrect file name
disp(['incorrect file name ' FN])
tA=NaN;tB=NaN;t=NaN;
    end
L=strcmp(DSI.Pass,'all') | (~isempty(strfind(FN,'_NSST'))&strcmp(DSI.Pass,'night')) | (~isempty(strfind(FN,'_SST'))&strcmp(DSI.Pass,'day'));
    if L
eval(['DSI.TimeA.' FIELD '(end+1)=tA;'])
eval(['DSI.TimeB.' FIELD '(end+1)=tB;'])
eval(['DSI.Time.'  FIELD '(end+1)=t;' ])
    end
        end %kt
end %Year
end %kFields

function DSI=dsi_modis_common2(DSI)
% merge sst and sst4 records

% In priciple, for each time there should be files corresponding to every field. 
% But some may be missing. Merge time arrays in order to find common time array T 
% for which at least one field is available.
% P contains pointers to elements in the original lists.
% Q containt a pointer to the array that contributed the value 
%[T,P,Q]=mergearrays(DSI.Time.Chlor,DSI.Time.SST,DSI.Time.SST_qual,...)
T=struct2cell(DSI.Time);
%assignin('base','DSI',DSI) % new
[T,P,Q]=mergearrays(T{:});

NT=length(T);
%Fields={'sst';'sst4'};
NV=length(DSI.Fields); % modified
for kv=1:NV
FIELD=DSI.Fields{kv};
eval(['F=DSI.URLFILE.' FIELD ';'])
eval(['DSI.URLFILE.' FIELD '=cell(NT,1);'])
    for k=1:NT
        if P(k,kv)>0
eval(['DSI.URLFILE.' FIELD '{k}=F{P(k,kv)};'])
        end
    end
eval(['F=DSI.URLPATH.' FIELD ';'])
eval(['DSI.URLPATH.' FIELD '=cell(NT,1);'])
    for k=1:NT
        if P(k,kv)>0
eval(['DSI.URLPATH.' FIELD '{k}=F{P(k,kv)};'])
        end
    end
end
DSI.Time=T;
TA=struct2cell(DSI.TimeA);
TB=struct2cell(DSI.TimeB);
DSI.TimeA=zeros(NT,1);
DSI.TimeB=zeros(NT,1);
for k=1:NT
DSI.TimeA(k)=TA{Q(k)}(P(k,Q(k)));
DSI.TimeB(k)=TB{Q(k)}(P(k,Q(k)));        
end

% add quality records (same files and times as for primary fields)
DSI.URLPATH.sst_qual=DSI.URLPATH.sst;
DSI.URLFILE.sst_qual=DSI.URLFILE.sst; 
DSI.URLPATH.sst4_qual=DSI.URLPATH.sst4;
DSI.URLFILE.sst4_qual=DSI.URLFILE.sst4;

NY=DSI.NY;NX=DSI.NX;
DSI.Latitude =  90-((1:1:NY)-0.5)'/NY*180;
DSI.Longitude=-180+((1:1:NX)-0.5)'/NX*360;

%DSI.Range.Time1=min(DSI.TimeA);
%DSI.Range.Time2=max(DSI.TimeB);
DSI.Range.Time1=min(DSI.Time);
DSI.Range.Time2=max(DSI.Time);
DSI.Attributes.time.units='days since 0000-01-01 00:00:00, MATLAB datenum';
%DSI.Attributes.time.valid_range=[min(DSI.Time) max(DSI.Time)];

DSI.Range.Latitude1=min(DSI.Latitude);
DSI.Range.Latitude2=max(DSI.Latitude);
DSI.Range.Longitude1=min(DSI.Longitude);
DSI.Range.Longitude2=max(DSI.Longitude);
% Make sure latitude and longitude are increasing arrays.
% In some datasets longitude may jump from 180 to -180 across the date
% line. Fix that by adding 360 degrees modulo 360.
DSI.Attributes.latitude.units='degrees_north';
%DSI.Attributes.latitude.valid_range=[min(DSI.Latitude) max(DSI.Latitude)];
DSI.Attributes.longitude.units='degrees_east';
%DSI.Attributes.longitude.valid_range=[min(DSI.Longitude) max(DSI.Longitude)];


S1=['time in seconds since reference_time;'];
S2=['missing values replaced with NaN;'];
S3=['variables converted to physical units.'];
%S4='';%['time is MATLAB datenum.'];
DSI.Readme=[S1 '/' S2 '/' S3];

DSI.Dimensions.latitude=length(DSI.Latitude);
DSI.Dimensions.longitude=length(DSI.Longitude);
DSI.Dimensions.time=length(DSI.Time);

dsname=lower(DSI.DataSetName);
eval(['DSI_' DSI.DataSetBranch '=DSI;'])
if ~exist(['dsi_' dsname '.mat'],'file') 
eval(['save dsi_' dsname '.mat DSI_' DSI.DataSetBranch ';'])
end
eval(['save dsi_' dsname '.mat -append DSI_' DSI.DataSetBranch ';'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S=htmldirlist(URL,EXT)
% gets the list of files with extension EXT in the directory pointed to by URL

% url location of the dataset files
%DSI.URLSITE='http://dods.jpl.nasa.gov';
%DSI.URLPATH='/opendap/sea_surface_temperature/modis/data/
%aqua/L3_mapped/sst/daily/04km/2006/365/
%A2006365.L3m_DAY_SST_4.bz2.html

S=urlread(URL);

% a fragment of urlread output
%<a href="A2006365.L3m_DAY_NSST_4.bz2.html">A2006365.L3m_DAY_NSST_4.bz2</a>
% match filenames between the leading and trailing strings
S1=[EXT '.html">'];
S2=['</a>'];

%SS=regexp(S,[EXT '\.html">' '.*?' '\.' EXT '</a>'],'match')'
SS=regexp(S,[EXT '\.html">' '\S*?' '\.' EXT '</a>'],'match')';
%SS=regexp(S,[EXT '\.html">' '\S+' '\.' EXT '</a>'],'match')'
%assignin('base','SS',SS) 


%cut off the the leading and trailing strings 
L1=length(S1);
L2=length(S2);
S=SS;
for k=1:length(S)
S1=SS{k};    
S{k}=S1(L1+1:end-L2);
end
S(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,p,q]=mergearrays(varargin)
% [a,p,q]=mergearrays(a1,a2,...,aM) 
%           merges two or more (M) integer arrays 
% input:    a1,a2,...,aM arrays to be merged. Each array aj may have different 
%               lengths.
% output:   a(i) merged array containing all possible values in a1,a2,...,aM
%               (without repetition). 
%               The pointers p(i,j) refer to the index in 
%               the original arrays a1,a2,...aM. The pointer is set to zero p(i,j)=0 
%               if the value is not available in the original array aj.
%               The pointer q(i) inicates which array contributed value first
%
% Example:
% a1    =[010;020;030;050;060;070;080;090;100];
% a2    =[030;040;060;070;080;090;100;120];
% a     =[010;020;030;040;050;060;070;080;090;100;120];
% p(:,1)=[  1;  2;  3;  0;  4;  5;  6;  7;  8;  9;  0];
% p(:,2)=[  0;  0;  1;  2;  0;  3;  4;  5;  6;  7;  8];
% q(:)  =[  1;  1;  1;  2;  1;  1;  1;  1;  1;  1;  2];
%
% if M=1 then a=a1 p=[1:length(a1)]'

A=varargin;
M=length(A); % number of arrays to merge 
NA=cellfun(@length,A); % length of each array
N=max(NA); % the length of the longest array 
Na=sum(NA); % the maximum possible length of array a
a=zeros(Na,1)+inf; % allocate memory, initialize values to inf
p=zeros(Na,M);     % allocate memory, initialize pointers to 0
q=zeros(Na,1);     % allocate memory, initialize pointers to 0
% convert the input cell array into matrix; initialize it to inf (because
% we use min) make one extra row 
B=zeros(N+1,M)+inf;
for k=1:M
B(1:NA(k),k)=A{k};
end
kk=ones(1,M); % pointers to the current element in the original array
ka=0;         % pointer to the element of the output array 

L=1;
while L % at least one of kk points to a valid value in aj 
    x=[];
    for k=1:M
        x=[x B(kk(k),k)];
    end    
    [v,i]=min(x);
% starting from the beginning of the arrays find the minimal value    
   %[v,i]=min([B(kk(1),1) B(kk(2),2) B(kk(3),3)]);  
    if ka==0 | v>a(ka)
% if there is nothing in a yet or the minimal value is larger than the last
% in a
        ka=ka+1;
        a(ka)=v; % then add the new value to a
        p(ka,i)=kk(i); % add the pointer to p
        q(ka)=i;       % add the pointer to q
        kk(i)=kk(i)+1; % and increment the pointer kk
    else 
 % otherwise just add the pointer. The same value is present in more than one array 
        p(ka,i)=kk(i); 
        kk(i)=kk(i)+1;
    end
    L=0; 
    for k=1:M
    L=L | kk(k)<=NA(k); % at least one of kk points to a valid value in aj 
    end
end %while
a=a(1:ka); % keep only the valid values, cut out extra memory 
p=p(1:ka,:);
q=q(1:ka,:);

function DSI=dsi_modis_files_daily_rt(DSI)
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/modis/data/aqua/L3_mapped/sst/8day/04km/2004/A20040012004008.L3m_8D_NSST_4.bz2?l3m_data[0:1:4319][0:1:8639],l3m_qual[0:1:4319][0:1:8639]
Fields={'sst';'sst4'};
for kFields=1:length(Fields)
FIELD=Fields{kFields};
ALG='';if strcmp(FIELD,'sst4') ALG='4'; end
for j=floor(datenum(2002,01,185)):floor(now)
    YYYY=datestr(j,'yyyy');
    Day=j-datenum(str2num(YYYY),01,00); %yearday
    DDD=num2str(Day,'%03d');   
%URLPATH='sst/daily/04km/2006/365/';
    URLPATH=['/opendap/sea_surface_temperature/modis/data/' DSI.Platform '/L3_mapped/' FIELD '/' DSI.TAvg '/' DSI.HRes '/' YYYY '/' DDD '/'];
%A2006365.L3m_DAY_SST_4.bz2.html
% obtain list of files with extention 'bz2'
%FNList=htmldirlist(URL,EXT)
%URL=[DSI.URLSITE URLPATH]
%FNList=htmldirlist(URL,'bz2')
% predict file names
FN1=[upper(DSI.Platform(1)) YYYY DDD '.L3m_DAY_' 'NSST' ALG '_' DSI.HRes(2) '.bz2'];
FN2=[upper(DSI.Platform(1)) YYYY DDD '.L3m_DAY_' 'SST'  ALG '_' DSI.HRes(2) '.bz2'];
FNList={FN1;FN2};
        for kt=1:length(FNList)
eval(['DSI.URLPATH.' FIELD '{end+1}=URLPATH;'])
eval(['DSI.URLFILE.' FIELD '(end+1)=FNList(kt);'])
%A2006365.L3m_DAY_SST_4.bz2      
FN=FNList{kt};
tA=datenum(str2num(FN(2:5)),01,str2num(FN(6:8)));
tB=tA+1;
    if ~isempty(strfind(FN,'_NSST')) % night images
t=tA;
    elseif ~isempty(strfind(FN,'_SST')) % day images add 12h to nominal time
t=tA+0.5;
    else % incorrect file name
disp(['incorrect file name ' FN])
tA=NaN;tB=NaN;t=NaN;
    end
    L=strcmp(DSI.Pass,'all') | (~isempty(strfind(FN,'_NSST'))&strcmp(DSI.Pass,'night')) | (~isempty(strfind(FN,'_SST'))&strcmp(DSI.Pass,'day'));
            if L
eval(['DSI.TimeA.' FIELD '(end+1)=tA;'])
eval(['DSI.TimeB.' FIELD '(end+1)=tB;'])
eval(['DSI.Time.'  FIELD '(end+1)=t;' ])
            end
        end %kt    
end %j
end %kFields

function DSI=dsi_modis_files_rt(DSI) % not ready
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/modis/data/aqua/L3_mapped/sst/8day/04km/2004/A20040012004008.L3m_8D_NSST_4.bz2?l3m_data[0:1:4319][0:1:8639],l3m_qual[0:1:4319][0:1:8639]
Fields={'sst';'sst4'};
for kFields=1:length(Fields)
FIELD=Fields{kFields};
for Year=2002:str2num(datestr(now,'yyyy'))
%for Year=2002:2008 %str2num(datestr(now,'yyyy'))
    YYYY=num2str(Year);
    URLPATH=['/opendap/sea_surface_temperature/modis/data/' DSI.Platform '/L3_mapped/' FIELD '/' DSI.TAvg '/' DSI.HRes '/' YYYY '/'];
% obtain list of files with extention 'bz2'
%FNList=htmldirlist(URL,EXT)
URL=[DSI.URLSITE URLPATH]
FNList=htmldirlist(URL,'bz2')
        for kt=1:length(FNList)
eval(['DSI.URLPATH.' FIELD '{end+1}=URLPATH;'])
eval(['DSI.URLFILE.' FIELD '(end+1)=FNList(kt);'])
%A20040012004008.L3m_8D_NSST_4.bz2     
FN=FNList{kt};
tA=datenum(str2num(FN(2:5)),01,str2num(FN(6:8)));
tB=datenum(str2num(FN(9:12)),01,str2num(FN(13:15)));
    if ~isempty(strfind(FN,'_NSST')) % night images
t=ceil((tA+tB)*0.5);
    elseif ~isempty(strfind(FN,'_SST')) % day images add 12h to nominal time
t=ceil((tA+tB)*0.5)+0.5;
    else % incorrect file name
disp(['incorrect file name ' FN])
tA=NaN;tB=NaN;t=NaN;
    end    
eval(['DSI.TimeA.' FIELD '(end+1)=tA;'])
eval(['DSI.TimeB.' FIELD '(end+1)=tB;'])
eval(['DSI.Time.'  FIELD '(end+1)=t;' ])
        end %kt
end %Year
end %kFields

