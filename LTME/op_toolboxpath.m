function varargout = op_toolboxpath(varargin)
% Adds interface directories and loaddap to matlab path.
%
% op_toolboxpath(opt)
% status = op_toolboxpath(opt)
%
% INPUT(S):
%
%   opt -- scalar, opt == 1 adds "loaddap", opt == 2 adds
%          interface paths
%
% OUTPUT(S):
%
%   status -- (optional) status of paths added or not, 0 if not
%
% (Not sure of the status of this function -- it is duplicated in
% "ocean_toolbox" by typing "ocean_toolbox -nojvm" or "-nogui".)
%
% OPeNDAP Science Team
% Copyright 2007, 2008
% $Version 2.0.1$

%==========================================================================
% Peter Cornillon
%
% REVISION HISTORY:
% 2007/04/XX 0.0.1 created, pcornillon
% 2007/05/07 0.0.2 modified, ceb
% 2008/03/01 2.0.0 released, ceb
% 2008/03/18 2.0.1 modified "help" section, ceb
%==========================================================================

% Get the loaddap directory name

if (nargin==1)
    input1 = varargin{1};
    if (~ischar(input1))
        input1 = num2str(input1);
    end
    AddLoaddap = input1;
    printout = 0;
elseif (nargin==2) %first input number, second input character string
    input1 = varargin{1};
    if (~ischar(input1))
        input1 = num2str(input1);
    end
    if strcmp(varargin{2},'-v') || strcmp(varargin{2},'-nojvm')
        printout = 1;
    else
        printout = 0;
    end
    AddLoaddap = input1;    
else
    error('User must enter either "1" or "2".')
end

% PRELIMINARIES.
quot = sprintf('%s',''''); %SINGLE QUOTE MARK
tt = quot; %RENAME FOR CONCISENESS.
if (ispc), slash='\'; end
if (isunix), slash='/'; end
if (ismac), slash='/'; end

if (ispc)
    location = 'C:\opendap\ml-structs\'; % Default directory
elseif (ismac)
    location = '/usr/local/bin/'; % Default directory
elseif (isunix)
    location = '/usr/local/bin/'; % Default directory
end
AltLoaddapDirectory = strcat(slash,'aaa',slash,'bbb',slash,'...',slash);

if (~exist('input1','var'))
disp(['Do you want to:']);
disp(['(1) Add "loaddap" to your path?']);
disp(['(2) Add OPeNDAP Toolbox mfiles to your path?']);
AddLoaddap = input(['Enter [1] or [2] ([n] to cancel): '],'s');
end

if (strcmp(AddLoaddap,'n'))
    disp('Program cancelled.')
    if nargout
    status = 0;
    varargout{1} = status;
    end
    return
end
    
if (strcmp(AddLoaddap,'1'))
    
if exist('loaddap.m','file')
    %disp('"loaddap" is already in the Matlab path. add "loaddap" canceled.');
    if nargout
    status = 1;
    varargout{1} = status;
    end
    return
end

if exist(strcat(location,'loaddap.m'),'file')
    question_text = ['"loaddap" found in ',location,'. ',...
        'Would you like to add this to your path?'];
    ButtonName = questdlg(question_text,'Add "loaddap" ...','Yes','No','Yes');
    if strcmp(ButtonName,'Yes')
        tmp = 'y';
    else
        tmp = 'n';
    end
%     disp(['"loaddap" found in ',location,'.']);
%     tmp = input(['Would you like to add this to your path? ',...
%         '[y] or [n]: '],'s');
    if strcmpi(tmp,'y')
        eval(['addpath ' location ';'])
        disp(['* Adding ' location ' to your Matlab path.'])
        disp(['The user might consider saving current path settings ',...
              'for future sessions.'])
        if nargout
        status = 1;
        varargout{1} = status;
        end
        return
    else %if strcmp(tmp,'n')
        %nothing
    end
end

% disp('*******************************************************************')
% disp('* loaddap must be in your Matlab path to use the OPeNDAP GUIs.')
% 
% QuestionString = ['* Provide a directory name for loaddap or n for none [/.../ or n]: '];
% AddLoaddap = input(['* Do you want to add ',LoaddapDirectory,' to your path? [y/n]: '],'s');
% if AddLoaddap == 'n'
%     LoaddapDirectory = input(QuestionString,'s');
%     if lower(LoaddapDirectory(1)) == 'n'
%         LoaddapDirectory = [];
%     end
% end

txt1 = '"loaddap" must be in your Matlab path to use the GUIs.';
txt2 = ' Please select the location of "loaddap":';
txt3 = '';
LoaddapDirectory = uigetdir(location,[txt1,txt2,txt3]);

%if ~isempty(LoaddapDirectory)
if LoaddapDirectory
    eval(['addpath ' LoaddapDirectory ';'])
    %%%disp(['* Adding ' LoaddapDirectory ' to your Matlab path.'])
else
    disp('add "loaddap" canceled.')
    if nargout
    status = 0;
    varargout{1} = status;
    end
    return
end

%%%disp('%==================================================================')
else %(strcmp(AddLoaddap,'2'))
%%%disp('%==================================================================')
%%%disp('* Adding OPeNDAP Toolbox mfiles to your Matlab path to use GUIs.')
% Now add all the subdirectories in the ml-opendap-toolbox directory to the 
% Matlab path. First get current directory.

tmp = mfilename('fullpath');
CurrentDirectory = tmp(1:end-length(mfilename));

DD = dir(CurrentDirectory);

for iNames=1:length(DD)
    if DD(iNames).isdir
    if DD(iNames).name(1) ~= '.' && isempty(findstr(DD(iNames).name,'OPeNDAP_ToolBoxPaths.m'))
        eval(['addpath(' tt CurrentDirectory DD(iNames).name tt ');'])
        if (printout)
            disp(['* Adding ' tt CurrentDirectory DD(iNames).name tt ' to your Matlab path.'])
        end
    end
    end
end
addpath(CurrentDirectory); %ADD BASE DIRECTORY.
%%%disp(['* Adding ' tt CurrentDirectory tt ' to your Matlab path.'])

%%%disp('%==================================================================')
end

if nargout
    status = 0;
    varargout{1} = status;
end

return
% 
% % Add directories to the Matlab path necessary to make OPeNDAP GUIs work.
% % See subfunction of a similar name in the main gui program for a similar
% % program.
% %
% % THE PURPOSE OF THIS ROUTINE IS TO ADD THE FOLLOWING DIRECTORIES
% % TO THE MATLAB PATH:
% %
% % (1) loaddap directory name
% % (2) ml-opendap-toolbox directory name
% %
% % IF THE USER WOULD LIKE TO SKIP A STEP, ENTER "n" WHEN PROMPTED.
% % THE "PATH" IS ONLY MODIFIED FOR THE CURRENT SESSION.
% %
% % EXAMPLE:
% %
% % >> op_toolboxpath
% % Do you want to:
% % (1) Add "loaddap" to your path?
% % (2) Add mfiles to your path?
% % Enter [1] or [2]: 2
% %
% % *******************************************************************
% % * Adding mfiles to your Matlab path to use the OPeNDAP GUIs.
% % * Adding 'C:\Documents and Settings\user1\My Documents\work\OPENDAP\SVN_GUIETTE\AIRS' to your Matlab path.
% % * Adding 'C:\Documents and Settings\user1\My Documents\work\OPENDAP\SVN_GUIETTE\CommonFunctions' to your Matlab path.
% % * Adding 'C:\Documents and Settings\user1\My Documents\work\OPENDAP\SVN_GUIETTE\GOES' to your Matlab path.
% % * Adding 'C:\Documents and Settings\user1\My Documents\work\OPENDAP\SVN_GUIETTE\HYCOM' to your Matlab path.
% % * Adding 'C:\Documents and Settings\user1\My Documents\work\OPENDAP\SVN_GUIETTE\MODIS' to your Matlab path.
% % * Adding 'C:\Documents and Settings\user1\My Documents\work\OPENDAP\SVN_GUIETTE\OAFlux' to your Matlab path.
% % * Adding 'C:\Documents and Settings\user1\My Documents\work\OPENDAP\SVN_GUIETTE\Pathfinder4km' to your Matlab path.
% % * Adding 'C:\Documents and Settings\user1\My Documents\work\OPENDAP\SVN_GUIETTE\SeaWinds' to your Matlab path.
% % *******************************************************************
% %
% %   (OPTIONAL)
% %   IF USED WITH INPUTS PROGRAM WILL
% %   1 -- ADD "LOADDAP" TO PATH
% %   2 -- ADD MFILES TO PATH
% %
% % OPeNDAP Science Team
% % Copyright 2007
% % $Version 1.00$
% 
% %==========================================================================
% % Peter Cornillon
% % Christian Buckingham
% % May 2007
% %
% % REVISION HISTORY:
% % 2007/04/XX PCC Created
% % 2007/05/07 CEB Edited
% %==========================================================================
% 
% % Get the loaddap directory name
% 
% if (nargin==1)
%     input1 = varargin{1};
%     if (~ischar(input1))
%         input1 = num2str(input1);
%     end
%     AddLoaddap = input1;
% elseif (nargin>1)
%     error('User must enter either "1" or "2".')
% end
% 
% % PRELIMINARIES.
% quot = sprintf('%s',''''); %SINGLE QUOTE MARK
% tt = quot; %RENAME FOR CONCISENESS.
% if (ispc), slash='\'; end
% if (isunix), slash='/'; end
% if (ismac), slash='/'; end
% 
% if (ispc)
%     location = 'C:\opendap\ml-structs\'; % Default directory
% elseif (ismac)
%     location = '/usr/local/bin/'; % Default directory
% elseif (isunix)
%     location = '/usr/local/bin/'; % Default directory
% end
% AltLoaddapDirectory = strcat(slash,'aaa',slash,'bbb',slash,'...',slash);
% 
% if (~exist('input1','var'))
% disp(['Do you want to:']);
% disp(['(1) Add "loaddap" to your path?']);
% disp(['(2) Add OPeNDAP Toolbox mfiles to your path?']);
% AddLoaddap = input(['Enter [1] or [2] ([n] to cancel): '],'s');
% end
% 
% if (strcmp(AddLoaddap,'n'))
%     disp('Program cancelled.')
%     return
% end
%     
% if (strcmp(AddLoaddap,'1'))
%     
% if exist('loaddap.m','file')
%     disp('"loaddap" is already in the Matlab path. add "loaddap" canceled.');
%     return
% end
% 
% if exist(strcat(location,'loaddap.m'),'file')
%     question_text = ['"loaddap" found in ',location,'. ',...
%         'Would you like to add this to your path?'];
%     ButtonName = questdlg(question_text,'Add "loaddap" ...','Yes','No','Yes');
%     if strcmp(ButtonName,'Yes')
%         tmp = 'y';
%     else
%         tmp = 'n';
%     end
% %     disp(['"loaddap" found in ',location,'.']);
% %     tmp = input(['Would you like to add this to your path? ',...
% %         '[y] or [n]: '],'s');
%     if strcmpi(tmp,'y')
%         eval(['addpath ' location ';'])
%         disp(['* Adding ' location ' to your Matlab path.'])
%         disp(['The user might wish to save the current path settings ',...
%               'for future sessions.'])
%         return
%     else %if strcmp(tmp,'n')
%         %nothing
%     end
% end
% 
% % disp('*******************************************************************')
% % disp('* loaddap must be in your Matlab path to use the OPeNDAP GUIs.')
% % 
% % QuestionString = ['* Provide a directory name for loaddap or n for none [/.../ or n]: '];
% % AddLoaddap = input(['* Do you want to add ',LoaddapDirectory,' to your path? [y/n]: '],'s');
% % if AddLoaddap == 'n'
% %     LoaddapDirectory = input(QuestionString,'s');
% %     if lower(LoaddapDirectory(1)) == 'n'
% %         LoaddapDirectory = [];
% %     end
% % end
% 
% txt1 = '"loaddap" must be in your Matlab path to use the GUIs.';
% txt2 = ' Please select the location of "loaddap":';
% txt3 = '';
% LoaddapDirectory = uigetdir(location,[txt1,txt2,txt3]);
% 
% %if ~isempty(LoaddapDirectory)
% if LoaddapDirectory
%     eval(['addpath ' LoaddapDirectory ';'])
%     disp(['* Adding ' LoaddapDirectory ' to your Matlab path.'])
% else
%     disp('add "loaddap" canceled.')
%     return
% end
% 
% disp('*******************************************************************')
% else %(strcmp(AddLoaddap,'2'))
% disp('*******************************************************************')
% disp('* Adding OPeNDAP Toolbox mfiles to your Matlab path to use GUIs.')
% % Now add all the subdirectories in the ml-opendap-toolbox directory to the 
% % Matlab path. First get current directory.
% 
% CurrentDirectory = strrep(mfilename('fullpath'),mfilename,'');
% 
% DD = dir(CurrentDirectory);
% 
% for iNames=1:length(DD)
%     if DD(iNames).isdir
%     if DD(iNames).name(1) ~= '.' && isempty(findstr(DD(iNames).name,'OPeNDAP_ToolBoxPaths.m'))
%         eval(['addpath(' tt CurrentDirectory DD(iNames).name tt ');'])
%         disp(['* Adding ' tt CurrentDirectory DD(iNames).name tt ' to your Matlab path.'])
%     end
%     end
% end
% addpath(CurrentDirectory); %ADD BASE DIRECTORY.
% disp(['* Adding ' tt CurrentDirectory tt ' to your Matlab path.'])
% 
% disp('*******************************************************************')
% end
% 
% return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function aa = ismac
% 
% dum = computer;
% aa = strcmpi(dum,'mac') || strcmpi(dum,'maci');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
function aa = ismac

dum = computer;
aa = strcmpi(dum,'mac') || strcmpi(dum,'maci');

%==========================================================================