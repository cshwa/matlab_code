function dataset = op_searchtool(search)
% Determines which datasets fall under a certain set of criteria.
%
% dataset = op_searchtool(search)
%
% INPUT(S):
%
%   database_name -- location of file to be read
%   search -- structure defining search criteria. Searches
%                     fields in the OPeNDAP dataset database.
%
% OUTPUT(S):
%
%   dataset -- names of datasets and GUIs which pass search criteria.
%              In the form of a structure.
%
% NOTE: quite complicated and used for main interface to enable/disable
% certain datasets.
%
% OPeNDAP Science Team
% Copyright 2007,2008,2009
% $Version 2.1.5$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2007/06/09 0.0.x created, ceb
% 2008/01/21 1.0.0 removed reference to SoftwareSettings file, preparing for
%                  ocean toolbox with different search tool, ceb
% 2008/02/28 2.1.2 prepare for distribution, ceb
% 2008/03/18 2.1.3 edited "help" text
% 2008/07/05 2.1.4 modified manner in which lat/lon are searched, ceb
% 2008/07/15 2.1.5 modified manner in which lat/lon are searched, ceb
% 2009/05/27       dates in csv files entered as yyyymmdd rather than in
%                  seconds since 1970, ms
%==========================================================================

% DATAFILE INFO.
%database_name = 'Master_Database.csv'; %'OPeNDAP_Ocean_Database.csv';

%load('SoftwareSettings.mat');
%database_name = S.master;

%==========================================================================

% PRELIMINATRIES.
variable = search.variable;
topic = search.topic;
datatype = search.datatype;
datasource = search.datasource;
date1 = search.date1;
date2 = search.date2;
lat1 = search.lat1;
lat2 = search.lat2;
lon1 = search.lon1;
lon2 = search.lon2;
res_time = search.res_time;
res_latlon = search.res_latlon;
keyword = search.keyword;
date_all_switch = search.date_all_switch;
latlon_all_switch = 0; %search.latlon_all_switch;
latlon_any_switch = ~latlon_all_switch;

% SWITCHES.
variable_switch = ~isempty(variable);
topic_switch = 0;
datatype_switch = 1;
datasource_switch = 1;
keyword_switch = ~isempty(keyword);
latitude_switch = ~isempty(lat1) && ~isempty(lat2);
longitude_switch = ~isempty(lon1) && ~isempty(lon2);
%geographic_switch = latitude_switch && longitude_switch
date_switch = ~isempty(date1) && ~isempty(date2);
% date_all_switch = 0;
% date_any_switch = 0;
date_any_switch = ~date_all_switch;

%==========================================================================
if isempty(datatype) || isempty(datasource)
%     dataset.name = [];
%     dataset.function = [];
    dataset = [];
    return%NO PURPOSE IN CONTINUING.
end
if isempty(variable)
    dataset = [];
    return%NO PURPOSE IN CONTINUING.
end
%==========================================================================
if any(ismember(variable,'sea surface temperature')) && ...
        ~any(ismember(variable,'temperature'))
    len = length(variable);
    variable{len+1} = 'temperature';
end
%==========================================================================
% Read the database.
load masterdatabase.mat
% %==========================================================================

% FLAGS.
[nvar,ncol] = size(data);
variable_flag = ones(1,nvar);
topic_flag = ones(1,nvar);
datatype_flag = zeros(1,nvar);
datasource_flag = zeros(1,nvar);
keyword_flag = ones(1,nvar);
latitude_flag = ones(1,nvar);
longitude_flag = ones(1,nvar);
date_flag = ones(1,nvar);
%==========================================================================
% VARIABLE.
if variable_switch
    
variable_flag = zeros(1,nvar);

icol = [3];
for ii = 1:nvar
    for jj = 1:length(icol)
        jj = icol(jj);
        dummy = lower(data{ii,jj});
        %dummy = strtok_apl(dummy,'; ');
        for kk = 1:length(variable)
            idx = strcmpi(dummy,variable{kk});
            if idx
                variable_flag(ii) = 1;
            end
        end%kk
    end%jj
end%ii
end%VARIABLE_SWITCH
%==========================================================================
% TOPIC.
if topic_switch

% % COMBINATIONS.
% %{air} --> 'Atmosphere'
% %{ocean} --> 'Ocean'
% %{air-sea interface} --> 'Air-Sea Interface'
% %{mixed layer} --> 'Mixed Layer'
% %{climate} --> 'Climate'
% %
% %{ocean,air-sea interface} --> {'Ocean','Climate'}
% %{ocean,mixed layer} --> {'Atmosphere','Climate'}
% 
% topic_flag = zeros(1,nvar);
% 
% tmp{1} = {'air'};
% tmp{2} = {'ocean'};
% tmp{3} = {'air_sea_interface'};
% tmp{4} = {'mixed_layer'};
% tmp{5} = {'climate'};
% tmp{6} = {'ocean','air_sea_interface'};
% tmp{7} = {'ocean','mixed_layer'};
% %tmp{8} = {'air','ocean','climate'};
% 
% topic = lower(topic);
% 
% icol = [0];
% for ii = 1:nvar
%     for jj = 1:length(icol)
%         jj = icol(jj);
%         dummy = lower(data{ii,jj});
%         dummy = strtok_apl(dummy,'; ');
%         text_string = '';
%         if isempty(setxor(dummy,tmp{1}))
%             text_string = {'atmosphere'};
%         end
%         if isempty(setxor(dummy,tmp{2}))
%             text_string = {'ocean'};
%         end
%         if isempty(setxor(dummy,tmp{3}))
%             text_string = {'air-sea interface'};
%         end
%         if isempty(setxor(dummy,tmp{4}))
%             text_string = {'mixed layer'};
%         end
%         if isempty(setxor(dummy,tmp{5}))
%             text_string = {'climate'};
%         end
%         if isempty(setxor(dummy,tmp{6}))
%             text_string = {'ocean','air-sea interface'};
%         end
%         if isempty(setxor(dummy,tmp{7}))
%             text_string = {'ocean','mixed layer'};
%         end
%         for kk = 1:length(topic)
%         idx = any(ismember(text_string,topic{kk}));
%         if idx
%             topic_flag(ii) = 1;
%         end
%         end%kk
%     end%jj
% end%ii
end%TOPIC_SWITCH
%==========================================================================
% DATATYPE.
if datatype_switch
datatype = lower(datatype);
icol = [6];
for ii = 1:nvar
    for jj = 1:length(icol)
        jj = icol(jj);
        dummy = lower(data{ii,jj});
        dummy = strrep(dummy,'_',' ');
        %dummy = strtok_apl(dummy,'; ');
        idx = any(ismember(datatype,dummy));
        if idx
            datatype_flag(ii) = 1;
        end
    end%jj
end%ii
end%DATATYPE_SWITCH
%==========================================================================
% DATASOURCE.
if datasource_switch
datatype = lower(datatype);
icol = [7];
for ii = 1:nvar
    for jj = 1:length(icol)
        jj = icol(jj);
        dummy = lower(data{ii,jj});
        dummy = strtok_apl(dummy,'; ');
        for kk = 1:length(datasource)
        idx = any(ismember(dummy,datasource{kk}));
        if idx
            datasource_flag(ii) = 1;
        end
        end%kk
    end%jj
end%ii
end%DATASOURCE_SWITCH
%==========================================================================
% KEYWORD.
if keyword_switch
keyword_flag = zeros(1,nvar);
%keyword = lower(keyword);
%irow = [1:nvar];
icol = [1,2,3,4,5];
%idx = strfind(data(irow,icol),keyword);
for ii = 1:nvar
    %dummy = strfind(eval(hdg{ii}),lower(keyword));
    for jj = 1:length(icol)
        jj = icol(jj);
        dummy = lower(data{ii,jj});
        idx = strfind(dummy,keyword);
        if ~isempty(idx)
            keyword_flag(ii) = 1;
        end
    end%jj
end%ii
end%KEYWORD_SWITCH
%==========================================================================
% LATITUDE. -- ALL & ANY --
if latitude_switch
latitude_flag = zeros(1,nvar);
buffer = 0;

if ischar(lat1) && ischar(lat2) %OLD FORMAT
lat1 = str2num(lat1);
lat2 = str2num(lat2);
end

icol_latstart = [8];
icol_latstop = [9];
for ii = 1:nvar
    lat_start = str2num(data{ii,icol_latstart});
    lat_stop = str2num(data{ii,icol_latstop});
    
    % Add buffer region.
    %none
    lat_flag = 0;
    
    % ALL %MODIFIED 2008/07/15
    if latlon_all_switch
        if lat_start<=lat1 && lat_stop>=lat2
            lat_flag = 1;
        end
    end
    % ANY %MODIFIED 2008/07/15
    if latlon_any_switch
        if lat1>=lat_start && lat2<=lat_stop
            lat_flag = 1;
        elseif (lat1>=lat_start && lat1<=lat_stop) && lat2>=lat_stop
            lat_flag = 1;
        elseif lat1<=lat_start && (lat2>=lat_start && lat2<=lat_stop)
            lat_flag = 1;
        elseif lat1<=lat_start && lat2>=lat_stop
            lat_flag = 1;
        end
    end
    if lat_flag
        latitude_flag(ii) = 1;
    end
end%ii

end%latitude_switch
%==========================================================================
% LONGITUDE. -- ALL & ANY --
if longitude_switch
longitude_flag = zeros(1,nvar);
buffer = 0;

if ischar(lon1) && ischar(lon2) %OLD FORMAT
lon1 = str2num(lon1);
lon2 = str2num(lon2);
end%ischar

icol_lonstart = [10];
icol_lonstop = [11];
for ii = 1:nvar
        
    lon_start = str2num(data{ii,icol_lonstart});
    lon_stop = str2num(data{ii,icol_lonstop});
    
    % HANDLE INTERNATIONAL DATELINE % MODIFIED 2008/07/15
    if (lon_start>0 && lon_stop<0)
        lon_stop = mod(lon_stop,360);
    end
    
    % Add buffer region.
    %none
    lon_flag = 0;
    
    % ALL %MODIFIED 2008/07/15
    if latlon_all_switch
        if lon_start<=lon1 && lon_stop>=lon2
            lon_flag = 1;
        end
    end
    % ANY %MODIFIED 2008/07/15
    if latlon_any_switch
        if lon1>=lon_start && lon2<=lon_stop
            lon_flag = 1;
        elseif (lon1>=lon_start && lon1<=lon_stop) && lon2>=lon_stop
            lon_flag = 1;
        elseif lon1<=lon_start && (lon2>=lon_start && lon2<=lon_stop)
            lon_flag = 1;
        elseif lon1<=lon_start && lon2>=lon_stop
            lon_flag = 1;
        end
    end
    if lon_flag
        longitude_flag(ii) = 1;
    end
end%ii

end%longitude_switch
%==========================================================================
% DATE. -- ALL --
if date_switch && date_all_switch
date_flag = zeros(1,nvar);
buffer = 3*86400;

if ischar(date1) && ischar(date2) %OLD FORMAT
tmp1 = [str2num(date1(1:4)),str2num(date1(5:6)),str2num(date1(7:8))];
tmp2 = [str2num(date2(1:4)),str2num(date2(5:6)),str2num(date2(7:8))];

date1 = op_epoch([tmp1,0,0,0]); %add hours min sec
date2 = op_epoch([tmp2,0,0,0]); %add hours min sec
end%ischar

% Add buffer region.
date1 = date1 - buffer;
date2 = date2 + buffer;

icol_start = [12];
icol_stop = [13];
for ii = 1:nvar
    time_start = data{ii,icol_start};
    time_stop = data{ii,icol_stop};
    time_start = str2num(time_start);
    if ischar(time_stop)
        if strcmpi(time_stop,'present')
            time_stop = op_epoch(datevec(now));
        else
            time_stop = str2num(time_stop);
        end
    end
    %time_start
    %time_stop both in seconds since 1970

    if time_start<=date1 && time_stop>=date2
        date_flag(ii) = 1;
    end
end%ii

end%DATE_SWITCH
%==========================================================================
% DATE. -- ANY --
if date_switch && date_any_switch
date_flag = ones(1,nvar);
buffer = 3*86400;

if ischar(date1) && ischar(date2) %OLD FORMAT
tmp1 = [str2num(date1(1:4)),str2num(date1(5:6)),str2num(date1(7:8))];
tmp2 = [str2num(date2(1:4)),str2num(date2(5:6)),str2num(date2(7:8))];

date1 = op_epoch([tmp1,0,0,0]); %add hours min sec
date2 = op_epoch([tmp2,0,0,0]); %add hours min sec
end%ischar

% Add buffer region.
%date1 = date1 - buffer;
date2 = date2 + buffer;

icol_start = [12];
icol_stop = [13];
for ii = 1:nvar
    time_start = data{ii,icol_start};
    time_stop = data{ii,icol_stop};
%     time_start = str2num(time_start);
    
    %A = data(:,iloc); %mlint shows an error here wrt ":" but no error.
    %B = data(:,iloc+1); %mlint shows an error here wrt ":" but no error.
%    A = strrep(time_start,'present',num2str(round(op_greg2epoch(datevec(now)))));
%    B = strrep(time_stop,'present',num2str(round(op_greg2epoch(datevec(now)))));
%    time_start = str2num(char(A));
%    time_stop = str2num(char(B));

% csv table is supposed to have dates in the form yyyymmdd 
% or present+-ndays
if strmatch(time_start,'present')
    time_start = eval(strrep(time_start,'present',num2str(datenum(now))));
else
    time_start = datenum(time_start,'yyyymmdd');
end
if strmatch(time_stop,'present')
    time_stop  = eval(strrep(time_stop ,'present',num2str(datenum(now))));
else
    time_stop  = datenum(time_stop ,'yyyymmdd');
end
% convert datenum to seconds since 1970
time_start=op_greg2epoch(datevec(time_start));
time_stop =op_greg2epoch(datevec(time_stop ));
    %time_start
    %time_stop both in seconds since 1970

    if time_start>date1 && time_start>date2
        date_flag(ii) = 0;
    end
    if time_stop<date1 && time_stop<date2
        date_flag(ii) = 0;
    end
end%ii

end%DATE_SWITCH
%==========================================================================
% ALL %MODIFIED 2008/07/05
if latlon_all_switch
    geo_switch = latitude_flag & longitude_flag;
end
% ANY %MODIFIED 2008/07/15
if latlon_any_switch
    geo_switch = latitude_flag & longitude_flag;
end
idx = find(variable_flag & ...
      topic_flag & datatype_flag & datasource_flag & keyword_flag ...
      & date_flag & geo_switch);
% [b, i, j] = unique(data(idx,5));
% dataset.name = b';
%%% dataset.function = data(idx(i),6)';
% dataset = data(idx(i),5);

dataset = unique(data(idx,5));
%==========================================================================
return %ENDOFFUNCTION
%==========================================================================
function toks=strtok_apl(str,delim)
%toks=strtok_apl(str,delim)
%
%to extract all tokens from a string within a human lifespan
%replacement to my STRTOKS which iteratively calls the iterative STRTOK
%
%str: a string to extract tokens from
%delim: an array of delimiters (char or byte), defaults to STRTOK defaults
%
%test cases: s={'',' ','a','a ',' a ','a b','a b c',' a b c '}
%				for i=1:length(s),toks=strtok_apl(s{i}),end
%
%arc 8/00
%arc 12/00 add logic for a single delimiter character (caused matrix-vector conversion)

%check for scalar string

if ~ischar(str);
   fprintf(1,'string must be character array\n');
end

if isempty(str)
   toks={};
   return
end

strsiz=size(str);
if sum(strsiz>1)>1
   fprintf(1,'String argument must be scalar \n');
   return
end

%if no delimiters are passed, use the default (as STRTOK)

if nargin<2
   delim=[9:13 32];
else
   if sum(size(delim)>1)>1
      fprintf(1,'Delimiter array must be scalar \n');
      return
   end
end
delsiz=size(delim);

%repmat to make str and delim 2D and same size

if strsiz(1)>1,s=str';else;s=str;end
if delsiz(1)>1,d=delim;else;d=delim';end
ns=length(s);
nd=length(d);

%Tony's trick instead of repmat(s,nd,1)...
sn=double(s);
s=sn(ones(1,nd),:);
d=d(:,ones(ns,1));

%find the non-delimiter characters in s

if length(delim)>1
    good=all(s~=d);	%1 if a good char, 0 if a delimiter
else
    good=s~=d;
end

%need to find the start chars and stop chars of each token
%calc diff(good), +1 transitions are starts, -1 are stops

dif=diff(good);
start=find(dif==1)+1;
stop=find(dif==-1);

%need to set start/end states (is first char a delim or token);
%the first state transition MUST be +1 and the last MUST be -1
%if this is not so, then the string must have begun/ended on a token

if good(1)==1, start=[1 start];end
if good(end)==1 stop=[stop length(good)];end

%extract the tokens

ntoks=length(start);
if ntoks==0
   toks={};
else   
   toks=cell(ntoks,1);
	for i=1:length(start),toks{i}=str(start(i):stop(i));end
end

return
%==========================================================================