function varargout = op_version(varargin)
% Obtain version number for OPeNDAP GUIs.
%
% output = op_version('toolbox')
% output = op_version('mfile',filename)
%
% op_version with no output arguments prints the Toolbox and/or associated
% interfaces to the command line.
%
% INPUT(S):
%
%   input -- (optional) -- can be either 'toolbox' or 'mfile'
%
%   filename -- if input = 'mfile', a second input must be specificed
%               as to which GUI is to be queried. For example, 'get_airs'
%               or 'get_seawinds'. If no file is specified
%               op_version returns the version numbers of all GUIs.
%            -- also, filename can be a cell array of filenames.
%
% OUTPUT(S):
%
%   output -- version number (string)
%
% EXAMPLE:
%
% >> op_version
%
% gives
%
% ================================================
% Ocean Toolbox         Version 2.0.1
% ================================================
% Function Name         Version Number
% =============         ================
% get_airs              Version 2.0.0
% get_goes              Version 2.1.0
% get_hycom             Version 2.1.0
% get_oaflux            Version 2.1.0
% get_pathfinder1km     Version 2.0.0
% get_pathfinder4km     Version 2.1.0
% get_seawinds          Version 2.1.1
% get_oceancolor        Version 2.0.0
% (main interface)      Version 2.0.0
%
% Note: this is somewhat of an ad-hoc set of code within the function, but
% it works. May want to upgrade this function to one which obtains a number
% rather than a string. For this reason, this function may change in the
% future.
%
% OPeNDAP Science Team
% Copyright 2007,2008
% $Version 2.0.3$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2007/05/01 0.0.x created, ceb
% 2008/01/22 1.0.0 edited, ceb
% 2008/01/23 1.1.0 added load of "included_interfaces", ceb
% 2008/03/01 2.0.2 released, ceb
% 2008/03/18 2.0.3 modified "help" section
%==========================================================================

% DEFINE TOOLBOX NAME.
toolbox_name = 'Ocean Toolbox';
space = char(' '*ones(1,9));
%space = char(' '*ones(1,9));

% DEFINE GUIS.
if 0
ii = 1;
mname{ii} = 'get_goes'; ii = ii+1;
mname{ii} = 'get_hycom'; ii = ii+1;
mname{ii} = 'get_oaflux'; ii = ii+1;
mname{ii} = 'get_pathfinder1km'; ii = ii+1;
mname{ii} = 'get_pathfinder4km'; ii = ii+1;
mname{ii} = 'get_seawinds'; %ii = ii+1;
else
load NetworkVerification.mat
% mname = included_interfaces;
len = length(interface_directories);
mname = cell(1,len);
for ii = 1:len
    directory = interface_directories{ii};
    mname{ii} = strcat('get_',lower(directory));
end%ii
end

if nargin == 0
    case_num = 3;
elseif nargin == 1
    tmp = varargin{1};
    if strcmpi(tmp,'toolbox')
        % READ VERSION NUMBER OF TOOLBOX.
        case_num = 1;
    elseif strcmpi(tmp,'mfile')
        case_num = 3;
    else
        filename = tmp;
        case_num = 2;
    end
elseif nargin == 2
    tmp = varargin{1};
    filename = varargin{2};
    if strcmpi(tmp,'mfile')
        case_num = 2;
    else
        error('If two inputs, first must be "mfile".')
    end
else
    error('Number of inputs cannot exceed two.')
end

if case_num == 1
    load('PointerToMainFunction.mat','main_function');
    aa = white_space(main_function);
elseif case_num == 2
    % READ VERSION NUMBER OF EVERY MFILE SPECIFIED BY USER.
    if iscell(filename)
        len = length(filename);
        for ii = 1:len
            aa{ii} = white_space(filename{ii});
        end%ii
    else
        aa = white_space(filename);
    end
elseif case_num == 3
    load('PointerToMainFunction.mat','main_function');
    aa{1} = white_space(main_function);
    %space = char(' '*ones(1,9));
    disp('================================================')
    disp([toolbox_name,space,aa{1}])
    disp('================================================')
    filename = mname;
    % READ VERSION NUMBER OF EVERY GUI.
    len = length(filename);
    for ii = 2:len+1
        aa{ii} = white_space(filename{ii-1});
    end%ii%END,LENGTH(MNAME)
end

if nargout == 0 && (nargin == 1 || nargin == 0)
    if iscell(aa)
        %space = char(' '*ones(1,9));
        txt1 = 'Function Name';
        tmptmp = length(toolbox_name) - length(txt1);
        tmp_other = char(' '*ones(1,tmptmp));
        txt2 = 'Version Number';
        disp([txt1,space,tmp_other,txt2])
        disp(['=============',space,tmp_other,'================'])
        len = length(space);
        for ii = 2:length(aa)
            %disp([filename{ii},', Version ',aa{ii}])
            dum1 = filename{ii-1};
            %dum2 = aa{ii};
            tmp1 = length(txt1)+len-length(dum1);
            space = char(' '*ones(1,tmp1+tmptmp));
            disp([filename{ii-1},space,aa{ii}])
        end%ii
        %op_version('mfile','gui_ocean_toolbox')
        tmp = op_version('mfile','gui_ocean_toolbox');
        dum1 = '(main interface)';
        tmp1 = length(txt1)+len-length(dum1);
            space = char(' '*ones(1,tmp1+tmptmp));
        disp([dum1,space,tmp])
    else
        %space = char(' '*ones(1,9));
        disp([toolbox_name,space,aa])
    end
elseif nargout == 0 && nargin == 2
    disp([filename,space,aa])
else
    varargout{1} = aa;
end

return %ENDOFFUNCTION
%==========================================================================
function info = white_space(filename)

% GET ROOT DIRECTORY.
%load('PointerToMainFunction.mat','main_function');
%mfile = main_function; %'OPeNDAP_ToolBox.m';
%tmp = which(mfile);
%tmp = strrep(tmp,mfile,'');
dummy = help(filename);
idx = strfind(dummy,'$');
if length(idx)~=2
    %error(['Could not find "$" sign in HELP part of code: ',filename])
end
if isempty(dummy)
    info = [];
    return
end
info = dummy(idx(1)+1:idx(2)-1);
%info = strrep(info,'Version ','')

return %ENDOFFUNCTION
