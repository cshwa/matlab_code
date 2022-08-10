function U=cat_airs(R)
%function U=Get_AIRS_URLs(DATE1,DATE2,DATEINCR,TAVG)
%U=Get_AIRS_URLs(R.DATE1,R.DATE2,R.DATEINCR,R.TAVG);
% R.DATEINCR==Inf to get first and last entry dates
U='';
if ~isfield(R,'TAVG') 
    disp('R.TAVG not specified')
    return
end

%URLCatBase='http://dods.gso.uri.edu/dods-3.4/nph-dods/catalog/';
URLSite='http://acdisc.sci.gsfc.nasa.gov/opendap/catalog/DatapoolCatalog/AIRS/';
URLPath='';URLFile='';
if strcmp(R.TAVG,'daily')   URLFile='AIRX3STD_005-cat.dat'; UNAME='AIRX3STD_005'; end 
if strcmp(R.TAVG,'weekly')  URLFile='AIRX3ST8_005-cat.dat'; UNAME='AIRX3ST8_005'; end 
if strcmp(R.TAVG,'monthly') URLFile='AIRX3STM_005-cat.dat'; UNAME='AIRX3STM_005'; end 

if isfield(R,'DATEINCR') & R.DATEINCR == Inf
CVAR='?year,month,day';
CDATE='';
else
CVAR='?DODS_URL,year,month,day';
    %% constrain date
CDATE=['&date(''' R.DATE1(1:4) '/' R.DATE1(5:6) '/' R.DATE1(7:8) ''',''' R.DATE2(1:4) '/' R.DATE2(5:6) '/' R.DATE2(7:8) ''')'];   
end
%
URL=[URLSite URLPath URLFile CVAR CDATE];
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

if isfield(R,'DATEINCR') & R.DATEINCR == Inf
    Time1=datenum(U.year(1),U.month(1),U.day(1));
    Time2=datenum(U.year(end),U.month(end),U.day(end));
    clear U; U.TimeMin=Time1;U.TimeMax=Time2;
end

else
U='';
end
