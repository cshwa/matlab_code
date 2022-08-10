function op_mergingtool(varargin)
% Merges structures obtained with OPeNDAP GUIs to user-friendly form.
%
% op_mergingtool
% op_mergingtool(prefixin)
% op_mergingtool(prefixin,number_range)
%
% op_mergingtool with no input arguments uses a predefined prefix,
% prefixin = 'opendap'.
%
% INPUT(S):
%
%   prefixin -- (optional) token phrase for the files to be converted, string
%       By default, assumed to be 'opendap'.
%   number_range -- (optional) two-element vector defining limits of conversion.
%       Opendap structures are labeled by number and increment each time a
%       variable is placed in memory. To convert only variables within a
%       certain number range, set number_range = [start_number,stop_number].
%
% OUTPUT(S):
%
%   [none] -- variables are placed in user's base workspace.
%
%   The resulting structure is a multidimensional structure with two types 
%   of data: dependent and independent.
%
%   (1) dependent -- data obtained from GUI and usually the data of interest.
%       This data is multidimensonal and has various dimensions which can
%       be found in the metadata field, usually something like 
%       dimensions = {'latitude','longitude','time','depth'}.
%
%   (2) independent -- data which describe the primary data. This includes
%       vectors .latitude, .longitude, .depth, .time, etc. as corresponding
%       to the dimensions above.
%
%    *  Metadata fields are included in this data structure because it is
%       believed they are helpful in describing the data. Information
%       included in this field are
%
%       name - string used to define output name; used in merging tool
%           temporal - string used to identify type of data, climatology or
%           time_series; used in merging tool
%       readme - general purpose communication string to help the user
%       url - internet protocol used to obtain the data (loaddap(url))
%       request - structure containing information from the request
%           made by the user; this request structure is formed by using the
%           gui and selecting variables, dates, locations, etc. of
%           interest; it can also be formed on the fly by the user, but
%           currently every dataset has a different 'request' format that it
%           uses
%       request_date - string used to identify when the request was
%           made; used in merging tool
%       
%       But, currently, there is some confusion of how to best
%       present this information to the user. Thus, this may change in the
%       future.
%
%   UNITS:
%
%       units for dependent data are specific to the dataset
%
%           HINT: either (1) look on the interface for that dataset and
%           place the cursor over variable of interest or (2) look in the
%           "metadata" field of the merged data structure ("attributes"),
%           and keep in mind that units are usually scaled to whole units
%           in the Toolbox.
%
%       units for independent data are specific to the dataset as well, but
%       are commonly the following:
%
%       latitude -- decimal degrees
%       longitude -- decimal degrees
%       time -- seconds past 1970/01/01 00:00:00
%            -- hint: use op_time.m or op_epoch2greg.m to obtain
%               a six element vector representing the year, month
%               hour, minute, second.
%               Example: round(op_epoch2greg(pathfinder1km.time))
%       depth -- either "levels" (1,2,3,...) or meters
%
% EXAMPLE:
%
% [Let's say an OPeNDAP interface was used to obtain data.]
%
% >> whos
%
%   Name               Size                    Bytes  Class
%   opendap_0001       1x1                    174428  struct array
%   opendap_0002       1x1                    174428  struct array
%   opendap_0003       1x1                    174428  struct array
%   opendap_0004       1x1                    174428  struct array
%   opendap_0005       1x1                    174428  struct array
%   opendap_0006       1x1                    174428  struct array
%   opendap_0007       1x1                    174428  struct array
%   request            1x1                      3064  struct array
%
% Grand total is 150869 elements using 1224060 bytes
%
% The user would merge the data in the following manner
%
% >> op_mergingtool
%
% or
%
% >> op_mergingtool('opendap')
%
% giving
%
% >> whos
%
%   Name                Size                    Bytes  Class
%   pathfinder1km       1x1                    981002  struct array
%   request             1x1                      3064  struct array
%
% Grand total is 121272 elements using 984066 bytes
%
% Or, let's say the user wishes to only merge variables "0003" through
% "0007". Then s/he would
%
% >> op_mergingtool('opendap',[3,7])
% >> whos
%
%   Name                Size                    Bytes  Class
%   opendap_0001        1x1                    174428  struct array
%   opendap_0002        1x1                    174428  struct array
%   pathfinder1km       1x1                    701898  struct array
%   request             1x1                      3064  struct array
%
% Grand total is 129858 elements using 1053818 bytes
%
% OPeNDAP Science Team
% Copyright 2007,2008
% $Version 2.0.4$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2006/12/01 0.0.1 created, ceb
% 2007/05/28 0.0.2 edited, ceb
% 2007/10/17 0.0.3 edited, ceb
% this version of the code implements a procedure in which the
% merging structure is set up beforehand by looping through the variables
% in the workspace to see that the input variables form well-structured
% output variables. if a previously merged structure is found in the user's
% workspace, this structure is parsed out into individual "opendap" struct-
% ures and then checked for merging ability with individual structures.
% 2008/03/01 2.0.3 released (don't know why version number wasn't 2.0.0)ceb
% 2008/03/18 2.0.3 edited documentation a little, ceb
% 2008/07/02 2.0.4 modified so as to retain unmerged variables, ceb
%==========================================================================
%ws = 'caller';
ws = 'base';
%==========================================================================
nan_replac = -99999;
%==========================================================================
% PRE-DEFINED ARGUMENTS.
default_prefx = 'opendap';
%==========================================================================
% DEFINE INDEPENDENT FIELDS.
ii=1;
indep_fields{ii}='latitude'; ii=ii+1;
indep_fields{ii}='longitude'; ii=ii+1;
indep_fields{ii}='time'; ii=ii+1;
indep_fields{ii}='depth'; ii=ii+1;
%==========================================================================
% DEFINE METADATA FIELD.
meta_field='metadata';
%==========================================================================
% DEFINE TOKEN FIELDS (UNION OF INDEPENDENT FIELDS AND METADATA FIELD).
token_fields = indep_fields; ii = length(indep_fields)+1;
token_fields{ii}=meta_field; ii=ii+1;
%==========================================================================
% DEFINE OPENDAP-SPECIFIC FIELDS (WITHIN METADATA).
% THESE WILL BE CELL ARRAYS.
ii=1;
opend_fields{ii}='url'; ii=ii+1;
opend_fields{ii}='attributes'; ii=ii+1;
opend_fields{ii}='name'; ii=ii+1;
opend_fields{ii}='temporal'; ii=ii+1;
opend_fields{ii}='reference_time'; ii=ii+1;
opend_fields{ii}='readme'; ii=ii+1;
opend_fields{ii}='user_friendly_time'; ii=ii+1;
opend_fields{ii}='request_date'; ii=ii+1;
opend_fields{ii}='request'; ii=ii+1;
%==========================================================================
% LOAD PRESET-SETTINGS FROM .MAT FILE.
tmp = mfilename('fullpath');
pname = tmp(1:end - length(mfilename));
fname = 'SoftwareSettings.mat';
load([pname,fname])
prompt_switch = S.merging_settings(4).value; %if 1 prompt user
%ONLY MERGE DATA WITH SAME REQUEST_ID.
merge_switch = S.merging_settings(1).value; %if 1 merge only structures
    %with same request date.
%==========================================================================
% PARSE INPUT(S).
if nargin == 0
    % no inputs, assume prefxin to be a certain form.
    prefxin = default_prefx;
elseif nargin == 1
    % only input is prefxin.
    prefxin = varargin{1};
elseif nargin == 2
    % Case 1 -- (1) input1 is prefxin (2) input2 is prefxout
    % Case 2 -- (1) input1 is prefxin (2) input2 is number_range
    % Case 3 -- (1) input1 is prefxin (2) input2 is date_format
    prefxin = varargin{1};
    variable2 = varargin{2};
    if (isnumeric(variable2))
        % Case 2
        if (length(variable2) == 2)
            range = variable2;
        % Case 3
        else
            date_format = variable2;
        end
    % Case 1
    else
        prefxout = variable2;
    end
elseif nargin == 3
    % Case 1 -- (1) input1 is prefxin (2) input2 is prefxout (3) input3
    % is number_range.
    % Case 2 -- (1) input1 is prefxin (2) input2 is prefxout (3) input3
    % is date_format.
    % Case 3 -- (1) input1 is prefxin (2) input2 is number_range (3)
    % input3 is date_format.
    prefxin = varargin{1};
    variable2 = varargin{2};
    variable3 = varargin{3};
    if (isnumeric(variable3))
        % Case 1
        if (length(variable3) == 2)
            prefxout = variable2;
            range = variable3;
        % Case 2
        else
            prefxout = variable2;
            date_format = variable3;
        end
    % Case 3
    else
        range = variable2;
        date_format = variable3;
    end
elseif nargin == 4
    prefxin = varargin{1};
    prefxout = varargin{2};
    range = varargin{3};
    date_format = varargin{4};
elseif nargin == 5
    % user desires to set the workspace manually ('caller' or 'base').
    prefxin = varargin{1};
    prefxout = varargin{2};
    range = varargin{3};
    date_format = varargin{4};
    ws = varargin{5};
    if strcmp(ws,'caller')
    elseif strcmp(ws,'base')
    else
        error('Input argument "ws" must be "caller" or "base".')
    end
else
    error('Number of input arguments cannot exceed 5.')
end
%==========================================================================
% FORCE FORMAT TO BE "0" -- NO DATE FORMAT OPTIONS AT THIS TIME.
% THIS HAS TO DO WITH THE TIME WHEN WE WERE CONSIDERING OUTPUTTING
% STRUCTURES TO THE USER'S BASE WORKSPACE WITH NAMES LABELED ACCORDING TO
% TIME.
date_format = 0;
%==========================================================================
% SETUP FOR WRITING NEW VARIABLES TO WORKSPACE.
%
% there is an attempt below to handle empty inputs:
% op_mergingtool(prefxin,[],number_range,date_format)
% op_mergingtool(prefxin,prefxout,[],date_format)
% op_mergingtool(prefxin,prefxout,number_range,[])

if (strcmp(prefxin(end),'_'))
    %prefxin = prefxin;
else
    prefxin = strcat(prefxin,'_');
end

if (exist('prefxout','var'))
    if (isempty(prefxout))
        clear prefxout
    elseif (strcmp(prefxout(end),'_'))
        %remove underscore
        prefxout = prefxout(1:end-1);
        name = prefxout;
    else
        %do nothing
        name = prefxout;
    end
end

if (exist('range','var'))
    if (isempty(range))
        range = [0000,9999];
    end
else
    range = [0000,9999];
end

if (exist('date_format','var'))
    if (isempty(date_format))
        date_format = 1;
    end
else
    date_format = 0;
end

%==========================================================================
% DETERMINE VARIABLES WITH PREFIX.
%
% NON-TOKEN FIELDS = PRIMARY DATA.
% TOKEN FIELDS = CONSIDERED SECONDARY DATA.
%==========================================================================
tt = sprintf('%s',''''); %SINGLE QUOTE MARK
%==========================================================================
%==========================================================================
% GET INVENTORY OF ALL VARIABLES.
[variable_in,variable_out,time,depth,temporal,...
    latitude,longitude,nlat,nlon,dependent] = ...
    getinventory(ws,prefxin,range,token_fields,meta_field,nan_replac);
%==========================================================================
% NO VARIABLES.
if isempty(variable_in)
    ws_tmp = ws;
    if strcmpi(ws,'caller')
        ws_tmp = strcat(ws,tt,'s');
    end
    disp(['There are no variables in the ',ws_tmp,' workspace to merge.'])
    %disp(['There are no variables in the ',ws_tmp,' workspace ',...
    %      'which start with ',tt,prefxin,tt,'.'])
    return
end
%==========================================================================
% BELOW, WE CHECK TO SEE IF OUTPUT VARIABLE EXISTS.
% IF IT DOES, WE PROMPT USER FOR (1) MERGE (2) OVERWRITE OR (3) CREATE NEW.
% IF USER SELECTS -- ATTEMPT TO MERGE -- THE PROGRAM PARSES THE POTENTIALLY
% LARGE MATRIX IN MEMORY INTO INDIVIDUAL OPENDAP STRUCTURES AND REPEATS THE
% PROCESS ABOVE OF TAKING INVENTORY OF ALL OPENDAP STRUCTURES.
%==========================================================================
% UNIQUE.
output_name = unique(variable_out);
output_name_prev = output_name;
n_output = length(output_name);
%==========================================================================
% LOOP OVER OUTPUT NAMES.
for ii = 1:n_output
    
varname_out = output_name{ii}; %TEMPORARY OUTPUT NAME.
exist_already = evalin(ws,strcat('exist(',tt,varname_out,tt,',',...
                    tt,'var',tt,')',';'));

% CHECK TO SEE IF OUTPUT VARIABLE EXISTS.
if exist_already
    number = op_getnum(varname_out); %DEFAULT = 0.
    suggestion = strcat(varname_out,'_',...
        num2str(number + 1,'%04d'));
    [choice,varname_out_new] = promptinput(varname_out,suggestion);
    % now varname_out_new is the new output variable name
    %======================================================================
    if strcmpi(choice,'merge') % attempt to merge data ...
        op_parsestruct(varname_out,prefxin,token_fields,opend_fields)
        [variable_in,variable_out,time,depth,temporal,...
        latitude,longitude,nlat,nlon,dependent] = ...
        getinventory(ws,prefxin,range,token_fields,meta_field,nan_replac);
        % UNIQUE.
        output_name = unique(variable_out);
        n_output = length(output_name);
    end
    %======================================================================
    varname_out = varname_out_new;
end
% UPDATE OUTPUT NAME.
output_name{ii} = varname_out;
end %ii -- index unique output name
%==========================================================================
% HANDLE NANs.
inan_flag = zeros(1,length(indep_fields));
for ii = 1:length(indep_fields);
    tmpf = indep_fields{ii};
    if strcmpi(tmpf,'time') || strcmpi(tmpf,'depth')
    eval(['inan','=','isnan(',tmpf,');'])
    if ~isempty(inan)
        inan_flag(ii) = 1;
    end
    eval([tmpf,'(inan)','=','nan_replac',';'])
    end
end
% time(isnan(time)==1) = nan_replac;
% depth(isnan(depth)==1) = nan_replac;

if 0
for ii = 1:length(latitude) % ii -- index latitude
    tmp = latitude{ii};
    tmp(isnan(tmp)==1) = nan_replac;
    latitude{ii} = tmp;
    tmp = longitude{ii};
    tmp(isnan(tmp)==1) = nan_replac;
    longitude{ii} = tmp;
end % ii -- index latitude
end
%==========================================================================
% LOOP OVER OUTPUT NAMES.
OUTPUT = cell(1,n_output);
for ii = 1:n_output
    
    tmptmp_output = output_name_prev{ii};
    tmp_output = output_name{ii}; %TEMPORARY OUTPUT NAME.
    %tmp_output = 'tmp_output'; %TEMPORARY OUTPUT NAME.
    
    % OBTAIN INPUT VARIABLES CORRESPONDING TO THIS OUTPUT VARIABLE.
    idx = find(strcmpi(tmptmp_output,variable_out));
    n_instances = length(idx);
    %NOTE THAT n_instances COULD BE USEFUL
    %LATER ON TO DETERMINE HOW MANY INPUT VARIABLES CORRESPOND TO A SINGLE
    %OUTPUT VARIABLE.
    str = strcat(tmp_output,'.nlat','=','[]',';'); eval(str);
    str = strcat(tmp_output,'.nlon','=','[]',';'); eval(str);
    str = strcat(tmp_output,'.latitude','=','[]',';'); eval(str);
    str = strcat(tmp_output,'.longitude','=','[]',';'); eval(str);
    str = strcat(tmp_output,'.time','=','[]',';'); eval(str);
    str = strcat(tmp_output,'.depth','=','[]',';'); eval(str);
    str = strcat(tmp_output,'.temporal','=','[]',';'); eval(str);
    str = strcat(tmp_output,'.dependent','=','[]',';'); eval(str);
    str = strcat(tmp_output,'.variable_in','=','[]',';'); eval(str);
    str = strcat(tmp_output,'.variable_out','=','[]',';'); eval(str);
    str = strcat(tmp_output,'.merge','=','[]',';'); eval(str);
    for jj = 1:n_instances
        str = strcat(tmp_output,'.nlat(jj)','=','nlat(idx(jj))',';');
        eval(str);
        str = strcat(tmp_output,'.nlon(jj)','=','nlon(idx(jj))',';');
        eval(str);
        str = strcat(tmp_output,'.latitude{jj}','=','latitude{idx(jj)}',';');
        eval(str);
        str = strcat(tmp_output,'.longitude{jj}','=','longitude{idx(jj)}',';');
        eval(str);
        str = strcat(tmp_output,'.time(jj)','=','time(idx(jj))',';');
        eval(str);
        str = strcat(tmp_output,'.depth(jj)','=','depth(idx(jj))',';');
        eval(str);
        str = strcat(tmp_output,'.temporal{jj}','=','lower(temporal{idx(jj)})',';');
        eval(str);
        str = strcat(tmp_output,'.variable_in{jj}','=','variable_in{idx(jj)}',';');
        eval(str);
        str = strcat(tmp_output,'.dependent{jj}','=','dependent{idx(jj)}',';');
        eval(str);
%         str = strcat(tmp_output,'.dependent_unique{jj}','=','tmp_dependent{idx(jj)}',';');
%         eval(str);
        str = [tmp_output,'.temporal{jj}','=','lower(',...
            'strrep(',tmp_output,'.temporal{jj},',tt,' ',tt,',',tt,'_',tt,'));'];
        eval(str);
    end
    eval(strcat(tmp_output,'.variable_out = lower(output_name{ii});'));
    eval(strcat('OUTPUT{ii}','=',tmp_output,';'));
end
%==========================================================================
% CHECK THAT THE DEPENDENT FIELDS (MAIN DATA) AND INDEPENDENT FIELDS
% (LATITUDE,LONGITUDE,TIME,DEPTH) MATCH UP BEFOREHAND.
% LOOP OVER EVERY VARIABLE FIRST TO DETERIMINE IF THERE IS AN ERROR.
%==========================================================================
problem = [];
K = 1;
% LOOP OVER OUTPUT VARIABLES.
for ii = 1:n_output % ii -- index output
    
    varname_in_vec = OUTPUT{ii}.variable_in;
    n_input = length(varname_in_vec);
    merge = ones(1,n_input); %ASSUME ALL VARIABLE SHOULD BE MERGED AT FIRST
    
    % DEFINE SWITCH FOR MERGING INPUT VARIABLES.
    OUTPUT{ii}.merge = merge;
    
    if n_input == 1
        continue %skip this section
    end
    
    % Determine most common grid.
    % The following code chooses a standard lat/lon vector pair
    % as the vectors against which to compare the data.
    nlat_common = mode(OUTPUT{ii}.nlat);
    nlon_common = mode(OUTPUT{ii}.nlon);
    ilat = find(OUTPUT{ii}.nlat == nlat_common);
    ilat = ilat(1); %take first
    ilon = ilat; %set ilon to be ilat so that get pair
    nlat_standard = OUTPUT{ii}.nlat(ilat);
    nlon_standard = OUTPUT{ii}.nlon(ilon);
    lat_standard = OUTPUT{ii}.latitude{ilat};
    lon_standard = OUTPUT{ii}.longitude{ilon};
    
    % Determine most common temporal-type.
    temporal_unique = unique(OUTPUT{ii}.temporal);
    ntemporal = length(temporal_unique);
    npts = zeros(1,ntemporal); %frequency of occur. of each string.
    for kk = 1:ntemporal
    idx = ismember(OUTPUT{ii}.temporal,temporal_unique{kk});
    npts(kk) = length(idx);
    end
    [dum,imax] = max(npts);
    if length(imax>1)
        imax = imax(1);
    elseif isempty(imax)
        imax = 1;
    end
    temporal_standard = temporal_unique{imax}; %most common
    
    % LOOP OVER INPUT VARIABLES.
    for jj = 1:n_input % jj -- index input
    
    tmp = OUTPUT{ii}.variable_in;
    varname_in = tmp{jj};
    varname_out = OUTPUT{ii}.variable_out;
    time = OUTPUT{ii}.time;
    depth = OUTPUT{ii}.depth;
    
    % CHECK FOR CASE WHEN THE LATITUDE AND LONGITUDE FIELDS DO NOT MATCH.
    lat1 = OUTPUT{ii}.latitude{jj};
    lon1 = OUTPUT{ii}.longitude{jj};
    lat_check = 0;
    if length(lat1) == nlat_standard
        if all(lat1 == lat_standard)
            lat_check = 1;
        end
    end
    lon_check = 0;
    if length(lon1) == nlon_standard
        if all(lon1 == lon_standard)
            lon_check = 1;
        end
    end
    if ~lat_check || ~lon_check
        problem(K).varname_in = varname_in;
        problem(K).varname_out = varname_out;
        problem(K).error = 'conflicting grid dimensions';
        problem(K).code = 0;
        K = K+1;
        merge(jj) = 0;
        continue
    end
       
    % CHECK FOR CASE WHEN THE INPUT STRUCTURE DOES NOT HAVE A TIME
    % OR DEPTH THAT IS AN ELEMENT/MEMBER OF THE FINAL OUTPUT STRUCTURE'S
    % TIME OR DEPTH FIELD.
    nn = length(time);
    if jj == 1
        idx = (2:nn);
    elseif jj == nn
        idx = (1:nn-1);
    else
        idx = [1:jj-1,jj+1:nn];
    end

    %IF NEITHER TIME NOR DEPTH MATCH ANY EXISTING TIME/DEPTH VALUES.
    if ~ismember(time(jj),time(idx)) && ~ismember(depth(jj),depth(idx))
        problem(K).varname_in = varname_in;
        problem(K).varname_out = varname_out;
        problem(K).error = 'conflicting time/depth dimensions';
        problem(K).code = 0;
        K = K+1;
        merge(jj) = 0;
        continue
    end
    
    % CHECK FOR CASE WHEN TEMPORAL IS NOT THE SAME.
    if ~strcmpi(OUTPUT{ii}.temporal{jj},temporal_standard)
        %DONT KNOW HOW TO DO THIS YET.
        problem(K).varname_in = varname_in;
        problem(K).varname_out = varname_out;
        problem(K).error = 'conflicting data (temporal) type';
        problem(K).code = 0;
        K = K+1;
        merge(jj) = 0;
        continue
    end
    end % jj -- index input
    
    % DEFINE SWITCH FOR MERGING INPUT VARIABLES.
    OUTPUT{ii}.merge = merge;
    
end % ii -- index output

clear dummy

% DETERMINE WHICH DEPENDENT VARIABLES ARE IN EACH OUTPUT STRUCTURE.
% LOOP OVER OUTPUT VARIABLES.
n_output = length(OUTPUT);
for ii = 1:n_output % ii -- index output
    n_input = length(OUTPUT{ii}.variable_in);
    dependent_unique = [];
    ll = 1;
    for jj = 1:n_input % jj -- index input variable
        if OUTPUT{ii}.merge(jj)
            dependent = OUTPUT{ii}.dependent{jj};
            for kk = 1:length(dependent)
                dependent_unique{ll} = dependent{kk};
                ll = ll+1;
            end
        end %merge?
    end
    OUTPUT{ii}.dependent_unique = unique(dependent_unique);
end

% HANDLE ERRORS.
% THIS MIGHT NEED SOME SPRUCING-UP.
if ~isempty(problem)
    disp(' ')
    disp('The following are a series of errors produced by ')
    disp('the data merging program: ')
    for K = 1:length(problem)
        disp(' ') %disp(['(',num2str(K),')'])
        disp(['output variable name: ',problem(K).varname_out])
        disp(['input variable name: ',problem(K).varname_in])
        disp(['error message: ',problem(K).error])
        disp(' ')
    end % K -- index problem
    disp('Choosing not to merge these data, continuing on')
    disp('with other variables.')
end

% METADATA.
% CONSTRUCT SKELETON VERSION OF OUTPUT VARIABLE.
metadata = [];
for mm = 1:length(opend_fields)
    tmpf = opend_fields{mm};
    eval(strcat('metadata','.',tmpf,'=','[]',';'))
end % mm -- index opednap specific fields

% MERGING
% LOOP OVER OUTPUT VARIABLES.
n_output = length(OUTPUT);
for ii = 1:n_output % ii -- index output
    
    varname_out = OUTPUT{ii}.variable_out;
    
    % OVERWRITE OR CREATE NEW.
    % CONSTRUCT SKELETON VERSION OF OUTPUT VARIABLE.
    jj = 1;
    dummy.metadata = metadata;
    name = lower(varname_out);
    temporal = unique(OUTPUT{ii}.temporal);
    if iscell(temporal) && length(temporal)==1
            temporal = temporal{1};
    end
    eval(['dummy','.','latitude',...
        '= OUTPUT{ii}.latitude{jj};'])
    eval(['dummy','.','longitude',...
        '= OUTPUT{ii}.longitude{jj};'])
    eval(['dummy','.','time',...
        '= sort(unique(OUTPUT{ii}.time));'])
    eval(['dummy','.','depth',...
        '= sort(unique(OUTPUT{ii}.depth));'])
    eval(['dummy','.',meta_field,'.','name',...
        '= name;'])
    eval(['dummy','.',meta_field,'.','temporal',...
        '= temporal;'])
    
    % PUT OUTPUT STRUCTURE IN BASE WORKSPACE.
    assignin(ws,varname_out,dummy)
    
    nlat = length(OUTPUT{ii}.latitude{jj});
    nlon = length(OUTPUT{ii}.longitude{jj});
    tmp_dependent = OUTPUT{ii}.dependent_unique;
    column3_sort = sort(unique(OUTPUT{ii}.time));
    column4_sort = sort(unique(OUTPUT{ii}.depth));
    nt = length(column3_sort);
    nz = length(column4_sort);
    dim = [nlat,nlon,nt,nz]; %dimension of dependent variables
    
    % PREALLOCATE SPACE FOR INCOMING DATA. PUT NANs IN MATRIX.
    assignin(ws,'opendap_temporary',dim);
    for jj = 1:length(tmp_dependent) % jj -- index dependent variables
        tmpf = tmp_dependent{jj};
        evalin(ws,strcat(varname_out,'.',tmpf,' = ',...
            'repmat(nan,','opendap_temporary',')',';'));
    end % jj -- index dependent variables
    evalin(ws,['clear ','opendap_temporary']);
    
    % FILL STRUCTURE.
    n_input = length(OUTPUT{ii}.variable_in);
    for jj = 1:n_input % jj -- index input variable
        
        varname_in = OUTPUT{ii}.variable_in{jj};
        
        % LOOP OVER OPENDAP-SPECIFIC FIELDS. PUT IN CELLS.
        % ASSUMED TO BE WITHIN METADATA.
        for mm = 1:length(opend_fields)
            tmpf = opend_fields{mm};
            if ~any(strcmpi(tmpf,{'name','temporal'}))
            tmp = evalin(ws,strcat(varname_in,'.',meta_field,...
                '.',tmpf,';'));
%             if strcmpi(tmpf,'readme')
%                 tmp = lower(strrep(tmp,' ','_'));
%             end
            tmpvec = evalin(ws,strcat(varname_out,'.',meta_field,...
                '.',tmpf,';'));
            if ~isempty(tmp)
                K = length(tmpvec);
                if iscell(tmp)
                for icell = 1:length(tmp)
                    tmpvec{K+icell} = tmp{icell};
                end
                else
                    tmpvec{K+1} = tmp;
                end
            end
            if jj == n_input % last input variable
                if ~isempty(tmpvec)
                if ~isstruct(tmpvec{1})
                tmpvec = unique(tmpvec);
                end
                end
            end
            assignin(ws,'opendap_temporary',tmpvec);
            evalin(ws,strcat(varname_out,'.',meta_field,'.',tmpf,'=',...
                'opendap_temporary',';'))
            end % not certain fields
        end % mm -- index opednap specific fields
        evalin(ws,['clear ','opendap_temporary']);
        
        tmp_third_dim = OUTPUT{ii}.time(jj);
        tmp_fourth_dim = OUTPUT{ii}.depth(jj);
        
        % MERGE OR NOT?
        if OUTPUT{ii}.merge(jj);
        
        % LOOP OVER DEPENDENT VARIABLES.
        for kk = 1:length(tmp_dependent) % kk -- index dependent variables
            
            tmpf = tmp_dependent{kk};
            
            % IF ISFIELD IN INPUT STRUCTURE.
            if evalin(ws,strcat('isfield','(',varname_in,',',...
                    tt,tmpf,tt,')',';')) %ISFIELD
                
                i3 = find(column3_sort == tmp_third_dim);
                i4 = find(column4_sort == tmp_fourth_dim);
                if length(column4_sort) < 2 %SINGLETON DIMENSION
                    index_str = strcat(':,:,',num2str(i3));
                else
                    index_str = strcat(':,:,',num2str(i3),',',...
                        num2str(i4));
                end
                
                % PLACE DATA IN LOCATION IN NEW MATRIX.
                % THESE FIELDS HAVE ALREADY BEEN SORTED.
                evalin(ws,strcat(varname_out,'.',tmpf,...
                    '(',index_str,')',' = ',...
                    varname_in,'.',tmpf,';'));
                
%                 % ADD "_variables" FIELD WHICH DESCRIBES
%                 % INDEPENDENT VARIABLES FOR EACH DEPENDENT VARIABLE.
%                 evalin(ws,strcat(varname_out,'.',meta_field,'.',...
%                     tmpf,'_variables',...
%                     ' = ',...
%                     varname_in,'.',meta_field,'.',...
%                     tmpf,'_variables',';'));
%                 % NOTE THAT THE ABOVE _variables FIELD MAY BE MODIFIED
%                 % BY THE FACT THAT THE "TIME" AND "DEPTH" FIELDS ARE ADDED
%                 % IN THIS NEW STRUCTURE. (SEE BELOW).
                
            end %ISFIELD?
            
        end % kk -- index dependent variables
        
        end % MERGE OR NOT?
        
        % CLEAR INPUT VARIABLE.
        %evalin(ws,['clear ',varname_in]); %commented 2008/07/02 CEB
                                           %replaced with code below
        tmp = cell(1,length(problem)); %tmp will be a list of "problem"
                                       %variables
        for K = 1:length(problem)
            tmp{K} = problem(K).varname_in;
        end % K -- index problem
        delete_flag = any(strcmp(varname_in,tmp)); %if 1 do not delete
        if ~delete_flag
            evalin(ws,['clear ',varname_in]);
        end
        
    end % jj -- index input variable
    
    % HANDLE NANs.
    removed_time = 0;
    removed_depth = 0;
    for jj = 1:length(indep_fields)
        tmpf = indep_fields{jj};
        if strcmpi(tmpf,'time') || strcmpi(tmpf,'depth')
            if inan_flag(jj) == 1
                % REPLACE -99999 TERM WITH NANs.
                evalin(ws,strcat(varname_out,'.',tmpf,...
                '(',varname_out,'.',tmpf,'==',num2str(nan_replac),')',...
                ' = ','nan',';'));
                
                % REMOVE FIELD IF ALL VALUES WITHIN FIELD ARE NAN.
                % (I.E., NOT PERTINENT TO DATA SET).
                removed = evalin(ws,...
                    ['all(isnan(',varname_out,'.',tmpf,'))']);
                evalin(ws,['if ',num2str(removed),...
                    ',',varname_out,'=',...
                    'rmfield(',varname_out,',',tt,tmpf,tt,'); end'])
                
                if removed && strcmpi(tmpf,'time')
                    removed_time = 1;
                end
                if removed && strcmpi(tmpf,'depth')
                    removed_depth = 1;
                end
                
            end
        end
    end
    if removed_time && ~removed_depth
        variables_text = {'latitude','longitude','depth'};
    elseif ~removed_time && removed_depth
        variables_text = {'latitude','longitude','time'};
    elseif removed_time && removed_depth
        variables_text = {'latitude','longitude','time'};
    else
        variables_text = {'latitude','longitude','time','depth'};
    end
    assignin(ws,'opendap_temporary',variables_text);
    
    % LOOP OVER DEPENDENT VARIABLES.
    for kk = 1:length(tmp_dependent) % kk -- index dependent variables

        tmpf = tmp_dependent{kk};

        % IF ISFIELD IN OUTPUT STRUCTURE.
        if evalin(ws,strcat('isfield','(',varname_out,',',...
                tt,tmpf,tt,')',';')) %ISFIELD

            % ADD "_variables" FIELD WHICH DESCRIBES
            % INDEPENDENT VARIABLES FOR EACH DEPENDENT VARIABLE.
            evalin(ws,strcat(varname_out,'.',meta_field,'.',...
                tmpf,'_variables',...
                ' = ',...
                'opendap_temporary',';'));

        end %ISFIELD?

    end % kk -- index dependent variables
    evalin(ws,['clear ','opendap_temporary']);
    
    % ORDER FIELDS.
    % TOP-LEVEL.
    ff = evalin(ws,strcat('fieldnames(',varname_out,')',';'));
    indep_fields_pr = intersect(ff,indep_fields);
    token_fields_pr = indep_fields_pr; ii = length(indep_fields_pr)+1;
        token_fields_pr{ii}=meta_field; ii=ii+1;
    depend_fields = setdiff(ff,token_fields);
    ll = 1;
    for kk = 1:length(depend_fields)+length(token_fields_pr)
        if kk <= length(depend_fields)
            ff{kk} = depend_fields{kk};
        else
            ff{kk} = token_fields_pr{ll};
            ll = ll+1;
        end
    end
    assignin(ws,'opendap_temporary',ff);
    evalin(ws,strcat(varname_out,' = ',...
                'orderfields(',...
                varname_out,',','opendap_temporary',')',';'));
    evalin(ws,['clear ','opendap_temporary']);
    
    % SECONDARY-LEVEL (METADATA FIELDS).
    ff = evalin(ws,strcat('fieldnames(',varname_out,'.',...
        meta_field,')',';'));
    %opend_fields_pr = intersect(ff,opend_fields);
    variab_fields = setdiff(ff,opend_fields);
    ll = 1;
    for kk = 1:length(variab_fields)+length(opend_fields)
        if kk <= length(variab_fields)
            ff{kk} = variab_fields{kk};
        else
            ff{kk} = opend_fields{ll};
            ll = ll+1;
        end
    end
    assignin(ws,'opendap_temporary',ff);
    evalin(ws,strcat(varname_out,'.',meta_field,' = ',...
                'orderfields(',varname_out,'.',meta_field,...
                ',','opendap_temporary',')',';'));
    evalin(ws,['clear ','opendap_temporary']);
    
    if strcmpi(tmpf,'time') || strcmpi(tmpf,'depth')
    eval(['inan','=','isnan(',tmpf,');'])
    if ~isempty(inan)
        inan_flag(ii) = 1;
    end
    eval([tmpf,'(inan)','=','nan_replac',';'])
    end
        
end % ii -- index output variable
        
%==========================================================================

return %[EOF]

%==========================================================================
function b = strtok_uri(string,dil)
% b = strtok_uri(string,D)
%
% This function is reminiscent of code by Art Croucher..
%
% The function takes a string and removes any occurrence of
% the delimeter D, returning a cell array of the input with
% delimeters removed. D is a character. A is a string. Handles up to
% 1000 occurrences of the delimiter.
%
% OPeNDAP Science Team
% May 2007

% Christian Buckingham
% March 2007

N=1000;
r=string;
b=[];
for ii=1:N;
    [tmp,r]=strtok(r,dil);
    %b=[b,tmp];
    b{ii}=tmp; %Art Croucher, jhu/apl
    if (isempty(r))
        return
    end
end %for

return %[EOF]
%==========================================================================
function varargout = promptinput(varname,varname_example)
% Prompts user for new input name or to overwrite/merge data.
%
% [choice,variable_name] = promptinput(conflicting_name,suggested_name)
%
% INPUT(S):
%
%   conflicting_name -- name of data structure in memory that conflicts
%       with data structure to be written (e.g., 'pathfinder4km')
%   suggested_name -- name of new data structure if user chooses 'create'
%
% OUTPUT(S):
%
%   choice -- 'overwrite','merge', or 'create'
%   variable_name -- user's choice of new variable name
%
% OPeNDAP Science Team
% Copyright 2007
% $Version 1.00$

%==========================================================================
% Christian Buckingham
% October 2007
%
% REVISION HISTORY:
% 2007/10/01 created, ceb
%==========================================================================
merge_switch = 1; %USE A .MAT FILE TO GET USER'S PREFERENCES

choice0 = 'Overwrite';
choice1 = 'Merge';
choice2 = 'Create';
%==========================================================================
tt = sprintf('%s',''''); %SINGLE QUOTE MARK
% ASK QUESTION.
if merge_switch
    choice = questdlg(['A variable with the same name ',tt,...
        varname,tt,' exists in the workspace. Would you like to ',...
        '(1) overwrite existing variable, ',...
        '(2) (attempt to) merge data with existing variable, ',...
        'or ',...
        '(3) create a new variable?'],...
        'Conflict: Variable Name',choice0,choice1,choice2,choice0);
else
    choice = choice1;
end

% EXAMINE USER'S CHOICE.
switch choice
    case choice0
        variable_name = varname;
    case choice1
        variable_name = varname;
    case choice2
        answer = inputdlg('Please enter a variable name, e.g. ',...
            ':',1,{varname_example});
        if ~isempty(answer) && ~isnumeric(answer)
            %NOTHING
        else
            warndlg(['Could not interpret input variable name. ',...
                'Please try again.']);
            %TRY AGAIN.
            answer = inputdlg('Please enter a variable name, e.g. ',...
            ':',1,{varname_example});
        end
        if isempty(answer) && isnumeric(answer)
            error(['Could not interpret input variable name. ',...
                'Exiting ...']);
        end
        variable_name = answer{1};

end%SWITCH

varargout{1} = lower(choice);
varargout{2} = variable_name;

return %[EOF]
%==========================================================================
function [variable_in,variable_out,time,depth,temporal,...
    latitude,longitude,nlat,nlon,dependent] = ...
    getinventory(ws,prefxin,range,token_fields,meta_field,nan_replac)
% Gets inventory of input variables with input prefix prefxin.
%
% INPUT(S):
%
%   ws -- 
%   prefxin --
%   range -- 
%   token_fields --
%   meta_field --
%   nan_replac --
%
% OUTPUT(S):
%
%   variable_in --
%   variable_out --
%   time --
%   depth --
%   temporal --
%   latitude --
%   longitude --
%   nlat --
%   nlon --
%   dependent --
%
% OPeNDAP Science Team
% Copyright 2007
% $Version 1.00$

%==========================================================================
% Christian Buckingham
% October 2007
%
% REVISION HISTORY:
% 2007/10/01 created, ceb
%==========================================================================

tt = sprintf('%s',''''); %SINGLE QUOTE MARK

% GET VARIABLE NAMES.
vars_in = evalin(ws,'whos;');

% DEFINE VARIABLES.
variable_in = [];
variable_out = [];
time = [];
depth = [];
temporal = [];
latitude = [];
longitude = [];
nlat = [];
nlon = [];
dependent = []; % DEPENDENT VARIABLES / NON-TOKEN
kk = 1;
% LOOP OVER ALL VARIABLES IN 'WS' WORKSPACE.
for ii=1:length(vars_in)

% CURRENT VARIABLE NAME.
tmp = vars_in(ii).name;
% IF LENGTH OF VARIABLE NAME EXCEEDS THE LENGTH OF THE PREFIX --> GOOD.
% ELSE --> DON'T CARE ABOUT.
if (length(tmp)>length(prefxin))
if (strcmp(tmp(1:length(prefxin)),prefxin))
    
    % STORE INPUT VARIABLE NAME.
    %vars(jj) = vars_in(ii);
    b = strtok_uri(tmp,'_');
    
    % NEW FEATURE: ONLY EXECUTE OVER A RANGE OF STRUCTURES.
    RRRR = b{end};
    RRRR = str2double(RRRR);
    
    % IF WITHIN RANGE SPECIFIED BY USER ...
    if (RRRR>=range(1) && RRRR<=range(2))
    
    %STORE VARIABLE NAME. 
    varname_in = tmp;
    variable_in{kk} = tmp;
    variable_out{kk} = lower(evalin(ws,strcat(varname_in,'.',...
        meta_field,'.','name',';')));%name{ii};
    temporal{kk} = lower(evalin(ws,strcat(varname_in,'.',...
        meta_field,'.','temporal',';')));%temporal{ii};
    if 1 %independent variables
        if evalin(ws,strcat(...
                'isfield(',varname_in,',',tt,'time',tt,');')...
                );
        time(kk) = evalin(ws,strcat(varname_in,'.',...
            'time',';'));%time(ii);
        else
            time(kk) = nan;
        end
        if evalin(ws,strcat(...
                'isfield(',varname_in,',',tt,'depth',tt,');')...
                );%depth(ii);
        depth(kk) = evalin(ws,strcat(varname_in,'.',...
            'depth',';'));%depth(ii);
        else
            depth(kk) = nan;
        end
    latitude{kk} = evalin(ws,strcat(varname_in,'.',...
        'latitude',';'));%latitude{ii};
    longitude{kk} = evalin(ws,strcat(varname_in,'.',...
        'longitude',';'));%longitude{ii};
    nlat(kk) = length(latitude{kk});
    nlon(kk) = length(longitude{kk});
    end
    
    % DEFINE DEPENDENT FIELDS.
    ll = 1;
    fnames = evalin(ws,strcat('fieldnames(',varname_in,');'));
    tmp_dependent = [];
    for jj = 1:length(fnames)
        tmpf = fnames{jj};
        if ~any(ismember(token_fields,tmpf)) %IF NOT A TOKEN FIELD.
            tmp_dependent{ll} = tmpf;
            ll = ll + 1;
        end
    end % jj -- index field names
    dependent{kk} = tmp_dependent;
       
    kk = kk + 1;
    
    % ELSE CONTINUE WITH CODE.
    else
        %MATLAB's "CONTINUE" FUNCTION INCREMENTS TO THE NEXT ITERATION
        %OF THE LOOP WITHOUT EXECUTING THE REMAINDER OF THE CODE.
        continue
    end

end%IF
end%IF
end%LOOP OVER VARIABLES -- ii

% NOTE -- AT THIS POINT THE FOLLOWING VECTORS HAVE THE SAME LENGTH:
%
% VARIABLE_IN
% VARIABLE_OUT
% TIME
% DEPTH
% LATITUDE
% LONGITUDE
% DEPENDENT
%
% WHERE THE VARIABLE_IN,VARIABLE_OUT,LATITUDE,AND LONGITUDE ARE
% CELL ARRAYS.

return %[EOF]