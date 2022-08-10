function op_parsestruct(varargin)
% Parse merged structure into individual opendap structures.
%
% op_parsestruct(varname)
%
% INPUT(S):
%
%   varname_in -- variable name of structure to be parsed
%   prefx -- string denoting output variable name prefix (e.g., opendap)
%   token_fields -- list of non-data fields (such as time, depth, latitude,
%       longitude, etc.)
%   opendap_fields -- opendap specific fields, assumed to be within
%       "metadata"
%
% OUTPUT(S):
%
%   [none] -- data placed in user's workspace
%             if prefix = 'opendap', then variables
%
%             opendap_0001
%             opendap_0002
%             ...
%
%             will be placed in the user's workspace.
%
% This program is used mostly for parsing out large data structures into
% individual "opendap_xxxx" structures prior to merging. This aids in
% determining if individual structures can be merged or not.
% May be useful in parsing out multidimensional structures (AIRS,HYCOM)
% into smaller ones -- that can then be merged (not yet determined).
%
% OPeNDAP Science Team
% Copyright 2007,2008
% $Version 2.0.3$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2007/10/01 0.0.x created, ceb
% 2007/12/04 0.0.x modified to not write "metadata" to workspace, ceb
% 2008/01/08 1.0.2 modified to write "metadata" to workspace, ceb
% 2008/03/01 2.0.2 released, ceb
% 2008/03/18 2.0.3 modified "help" text, ceb
%==========================================================================

ws = 'base';
%ws = 'caller';

if nargin == 4
    varname_in = varargin{1};
    prefx = varargin{2};
    token_fields = varargin{3};
    opend_fields = varargin{4};
elseif nargin == 3
    varname_in = varargin{1};
    prefx = varargin{2};
    token_fields = varargin{3};
elseif nargin == 2
    varname_in = varargin{1};
    prefx = 'opendap';
    token_fields = varargin{3};
elseif nargin == 1
    varname_in = varargin{1};
    prefx = 'opendap';
end

% DEFINE TOKEN FIELDS.
if ~exist('token_fields','var')
ii=1;
token_fields{ii}='latitude'; ii=ii+1;
token_fields{ii}='longitude'; ii=ii+1;
token_fields{ii}='time'; ii=ii+1;
token_fields{ii}='depth'; ii=ii+1;
token_fields{ii}='name'; ii=ii+1;
token_fields{ii}='url'; ii=ii+1;
token_fields{ii}='attributes'; ii=ii+1;
token_fields{ii}='temporal'; ii=ii+1;
token_fields{ii}='readme'; ii=ii+1;
token_fields{ii}='request_id'; ii=ii+1;
token_fields{ii}='metadata'; ii=ii+1;
token_fields{ii}='request_date'; ii=ii+1;
end

% DEFINE OPENDAP-SPECIFIC FIELDS (WITHIN METADATA).
% THESE WILL BE CELL ARRAYS.
if ~exist('opend_fields','var')
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
end
%==========================================================================
tt = sprintf('%s',''''); %SINGLE QUOTE MARK

% Generate template for token fields.
meta_field = 'metadata';

% time
if evalin(ws,strcat('isfield(',varname_in,',',tt,'time',tt,');'));
time = evalin(ws,strcat(varname_in,'.','time',';'));
else
time = nan;
end
if evalin(ws,strcat('isfield(',varname_in,',',tt,'depth',tt,');'));
depth = evalin(ws,strcat(varname_in,'.','depth',';'));
else
depth = nan;
end
latitude = evalin(ws,strcat(varname_in,'.','latitude',';'));
longitude = evalin(ws,strcat(varname_in,'.','longitude',';'));

template.time = time;
template.depth = depth;
template.latitude = latitude;
template.longitude = longitude;

template2 = template;

% LOOP OVER OPENDAP-SPECIFIC FIELDS. PUT IN CELLS.
% ASSUMED TO BE WITHIN METADATA.
%==========================================================================
% uncomment below to write metadata as well as data -- will increase the
% size of the merged data structure to possibly three times as large ...
%==========================================================================
% for mm = 1:length(opend_fields)
%     tmpf = opend_fields{mm};
%     tmp = evalin(ws,strcat(varname_in,'.',meta_field,'.',tmpf,';'));
%     eval(strcat('template','.',meta_field,'.',tmpf,'=',...
%         'tmp',';'))
% end % mm -- index opednap specific fields
%==========================================================================
% the following lines of code are only here because we commented out the
% above; namely we have removed the ability to write certain (possibly
% large) fields when parsing the data. The reason we have done this is that
% when parsing very large matrices in order to merge individual "opendap"
% data structures with them, we repeat the data so many times that the size
% of the resulting merged structure can be almost 3 times as big as the
% data structure--even when the data is the same. (i.e., it is only the
% metadata that continues to be repeated. The true solution to this problem
% is not to parse out the merged structures at all but examine the
% potential of merging an opendap data structure with an already
% existing merged structure. A separate function would need to be written
% for this purpose.) 
%==========================================================================
for mm = 1:length(opend_fields)
    tmpf = opend_fields{mm};
    tmp = evalin(ws,strcat(varname_in,'.',meta_field,'.',tmpf,';'));
    eval(strcat('template2','.',meta_field,'.',tmpf,'=',...
        'tmp',';'))
end % mm -- index opendap specific fields
%==========================================================================
for mm = 1:length(opend_fields)
    tmpf = opend_fields{mm};
    if strcmpi(tmpf,'name') || strcmpi(tmpf,'temporal')
    tmp = evalin(ws,strcat(varname_in,'.',meta_field,'.',tmpf,';'));
    eval(strcat('template','.',meta_field,'.',tmpf,'=',...
        'tmp',';'))
    else
    tmp = '';
    eval(strcat('template','.',meta_field,'.',tmpf,'=',...
        'tmp',';'))
    end
end % mm -- index opednap specific fields
%==========================================================================


nvar = length(time)*length(depth); %these are how many "opendap_xxxx"
    %variables are required to parse the entire matrix.

pathname = mfilename('fullpath');
pathname = strrep(pathname,mfilename,'');
filename = 'opendap_temporary.mat';

lastreqnum = op_getnum(prefx,ws);
ivar = lastreqnum;
ivar_init = ivar + 1; %first variable written to workspace
                      %(contains metadata information)

% Loop over non-token fields.
fnames = evalin(ws,strcat('fieldnames(',varname_in,');'));
tmp = '';
for ii = 1:length(fnames)
    tmpf = fnames{ii};
    if ~any(ismember(token_fields,tmpf)) %IF NOT A TOKEN FIELD.
        % Parse out data in two-dimensional fields.
        varname_in_exist = evalin(ws,strcat('exist','(',tt,...
            varname_in,tt,',',tt,'var',tt,')',';'));
        %if ~varname_out_exist
        %    evalin(ws,['load ',strcat(pathname,filename)]);
        %end
        data = evalin(ws,strcat(varname_in,'.',tmpf,';'));
        evalin(ws,[varname_in,'.',tmpf,'=','[]',';']);
        %evalin(ws,['save ',strcat(pathname,filename),' ',varname_in]);
        %evalin(ws,['clear ',varname_in]);
        ss = size(data);
        % Loop over time.
        if length(ss)==2, ss = [ss,1,1]; end
        if length(ss)==3, ss = [ss,1]; end
        for jj = 1:ss(3) %time
        for kk = 1:ss(4) %depth
            ivar = ivar + 1;
            tmp_time = time(jj);
            tmp_depth = depth(kk);
            data_slice = data(:,:,jj,kk);
            % if first variable, put in full metadata
            if ivar == ivar_init
                tmp = template2;
            % else use only a limited amount of metadata
            else
                tmp = template;
            end
            tmp.time = tmp_time;
            tmp.depth = tmp_depth;
            eval(strcat('tmp.',tmpf,'=','data_slice',';'));
            varname_out = strcat(prefx,'',num2str(ivar,'%04d'));
            if nvar > 9999
                varname_out = strcat(prefx,'',num2str(ivar,'%05d'));
            end
            assignin(ws,varname_out,tmp);
            clear data_slice
            clear tmp
        end % kk -- index depth
        end % jj -- index time
        clear data
    end%not a token field
end

%delete(strcat(pathname,filename))
evalin(ws,['clear ',varname_in]);
return %[EOF]
%==========================================================================