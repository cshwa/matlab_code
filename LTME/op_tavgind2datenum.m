function [tA tB tC]=op_tavgind2datenum(k,R)
% for various temporal averages find the begining, end, and middle 
% of the averaging interval in days since Jan-01 00:00:00 of base year
% Input:
% k is the index (or vector of indecies) of the averaging period.
% First interval is 1, e.g. week 1, week 2, etc.
% If request spans new year week number will exceed 52 
% 
% R.BaseYear, R.TAvg

% base year of the climatology
if ~isfield(R,'BaseYear') y0=0;
else
    if ischar(R.BaseYear) 
        y0=str2num(R.BaseYear);
    else
        y0=R.BaseYear;
    end
end
TAVG=lower(R.TAvg);

tA=NaN;tB=NaN;
if ~isempty(strmatch(TAVG,{'daily','day','1day'}))
    NT=datenum(y0,12,31)-datenum(y0,1,1)+1;
    y=floor((k-1)/NT);j=rem(k-1,NT);
    tA=datenum(y0+y,1,j+1);tB=tA+1;    
elseif ~isempty(strmatch(TAVG,{'5day','pentad'}))
    y=floor((k-1)/73);j=rem(k-1,73);
    tA=datenum(y0+y,1,j*5+1);tB=tA+5;
elseif ~isempty(strmatch(TAVG,{'7day','week','weekly'}))
    y=floor((k-1)/52);j=rem(k-1,52);
    tA=datenum(y0+y,1,j*7+1);tB=tA+7;
elseif ~isempty(strmatch(TAVG,{'8day','octad'}))
    y=floor((k-1)/45);j=rem(k-1,45);
    tA=datenum(y0+y,1,j*8+1);tB=tA+8;
elseif ~isempty(strmatch(TAVG,{'monthly','month','30day'}))
    y=floor((k-1)/12);j=rem(k-1,12);
    tA=datenum(y0+y,j+1,1);tB=datenum(y0+y,j+2,1);
elseif ~isempty(strmatch(TAVG,{'seasonal','season','90day'}))
    y=floor((k-1)/4);j=rem(k-1,4);
    tA=datenum(y0+y,j*3+1,0);tB=datenum(y0+y,j*3+2,0); 
elseif ~isempty(strmatch(TAVG,{'annual','yearly','year','365day','12month'}))
    tA=datenum(k-1,1,1);tB=datenum(k,1,1);
end
tC=0.5*(tA+tB);
