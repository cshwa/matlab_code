function val=op_getnum(varargin)
% Gets the highest-numbered "opendap_xxxx" variable in the workspace.
%
% num = op_getnum
% num = op_getnum(prefix)
% num = op_getnum(prefix,ws)
%
% INPUT(S):
%
%   prefix -- (optional) prefix of variables for which one wants
%             the number, by default prefix = 'opendap'
%          -- note: the number is assumed to be an integer
%   ws -- desired workspace, 'caller' or ws= 'base', by default
%             ws = 'base'
%
% OUTPUT(S):
%
%   num -- last opendap structure number in base workspace.
%
% OPeNDAP Science Team
% June 2007, 2008
% $Version 2.0.0$

%==========================================================================
% Christian Buckingham
% Meri Sheremet
%
% REVISION HISTORY:
% 2007/06/15 0.0.x created, ceb, msheremet
% 2008/03/01 2.0.0 released
%==========================================================================

% Parse inputs.
if nargin == 0
    prefixin = 'opendap';
    ws = 'base';
elseif nargin == 1
    prefixin = varargin{1};
    ws = 'base';
elseif nargin == 2
    prefixin = varargin{1};
    ws = varargin{2};
else
    error('Number of input arguments must be 1 or 2.')
end

% DETERMINE NUMBER OF VARIABLES WITH PREFIX.
if (strcmp(prefixin(end),'_'))
    %prefixin = prefixin;
else
    prefixin = [prefixin,'_'];
end

% Check to see if structure in workspace.
vars_in = evalin(ws,'whos'); %'base'
len = length(prefixin);
%count = zeros(1,length(vars_in));
jj=1;
for ii=1:length(vars_in)
    tmp = vars_in(ii).name;
    if (length(tmp)>len)
    if (strcmp(tmp(1:len),prefixin))
        %vars(jj) = vars_in(ii);
        %variables{jj} = tmp;
        tmp = tmp(len+1:end);
        tmp = str2double(tmp);
        if isnan(tmp), tmp = 0; end
        count(jj) = tmp; %str2double(tmp);
        jj = jj + 1;
    end
    end
end%ii
nn = jj - 1;
if (nn < 1)
    val = 0;
    return
end
count_asc = sort(count);
val = count_asc(nn); %use num2str(val,'%04d') to convert to text string.
%disp(['Use "num2str(value,"%04d") to convert to text string.'])

return %endoffunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function kreq=LatestReqNum
% % find the latest opendap request number
% a=evalin('base','who');k=strmatch('opendap_',a);
% if isempty(k) k=0; 
% else 
% a=a{k(end)};
% ki=strfind(a,'_');
%     if length(ki)==1 ki(2)=length(a)+1; end
% k=str2num(a(ki(1)+1:ki(2)-1));
% end
% % latest request number
% kreq=k;
% %KREQ=num2str(kreq,'%04d');