function op_tooltipwrap(hObject,texto,nCols)
% Formats a "tooltipstring" in Matlab's graphical user interfaces.
%
% Author: Joaquim Luis
%
% Note: Matlab needs to have a "time" property that allows the user
% to set the duration for which the tooltipstring is up ... java interfaces
% use a thing called setDelayTime but currently Matlab does not support
% such a feature.
%
% op_tooltipwrap(h,text,ncol)
%
% See documentation below for more information.
%
% INPUT(S):
%
%   h -- handle to uicontrol object
%   text -- text string to format
%   ncol -- number of columns to use in formatting tooltipstring.
%
% OUTPUT(S):
%
%   [none]
%
% OPeNDAP Science Team
% Copyright 2007,2008
% $Version 2.0.1$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2007/01/01 0.0.x edited, ceb
% 2008/03/01 2.0.0 released
% 2008/03/18 2.0.1 modified "help" section
%==========================================================================

%TOOLTIPWRAP Uses the information in TEXTO as the 'TooltipString' for the HOBJECT UI Control.
%
% TEXTO can be a linestring, a cell array of strings or the name of a file containing
% the desired text info.
% NCOLS, is the width in chars of the wrapped TEXTO matrix. If absent it defaults to 20.
% When TEXTO is a filename and NCOLS = Inf, the file contents will be used with no change
% (that is, no text wrapping)
%
%   Example: Run this function with a bad usage and leave the pointer over the OK button.
%           tooltipwrap(1,1)

%   AUTHOR
%       Joaquim Luis    02-July-2006
%       w3.ualg.pt/~jluis/mirone

% Error check
msg = [];
if (nargin < 2 || nargin > 3)
    msg = 'Wrong number of arguments';
elseif (~ishandle(hObject) || ~strcmp(get(hObject,'Type'),'uicontrol') )
    msg = 'First argument is not a valid uicontrol handle.';
elseif ( ~(ischar(texto) || iscell(texto)) )
    msg = 'Second argument must be a line string or a cell array of chars.';
end

if (~isempty(msg))
    h = errordlg(['Tooltipwrap: ' msg],'Error');
    hButt = findobj(h,'Style','pushbutton');
    tooltipwrap(hButt,{msg},10)
    return
end
if (nargin ~= 3),   nCols = 20;     end
    
% Try to open 'texto' as a file. If we can, use the file contents as 'TooltipString'
if (~iscell(texto))
    fid = fopen(texto);
    if (fid > 0)                % Yes, 'texto' is a file name
        str = fscanf(fid,'%c');  fclose(fid);
        if (~isinf(nCols))
            str = textwrap({str},nCols);
        end
    else                        % No, 'texto' is a text string
        str = textwrap({texto},nCols);
    end
else
    str = textwrap(texto,nCols);
end

if (iscell(str))
    str_tip = [];
    for (i = 1:length(str)-1)
        str_tip = sprintf([str_tip str{i} '\n']);
    end
    str_tip = sprintf([str_tip str{end}]);
else
    str_tip = str;
end
set(hObject,'TooltipString',str_tip)
