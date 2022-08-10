function readme_new = op_format_readme(readme)
% Formats "readme" section of opendap structure.
%
% readme_new = op_format_readme(readme)
%
% INPUT(S):
%
%   readme -- 1 by N element vector containing strings separated
%          by "/" marks.
%
% OUTPUT(S):
%
%   readme_new -- 1 by N element cell array containing strings,
%              each statement in "readme" separated by slashes
%              is placed in its own cell.
%
% If input is a structure, then output is a structure
% with .readme field formatted.
%
% OPeNDAP Science Team
% Copyright 2007, 2008
% $Version 2.0.0$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2007/05/01 0.0.x created, ceb
% 2008/03/01 2.0.0 released, ceb
%==========================================================================

% DEFINE WORKSPACE FOR NO-INPUT CALLING SEQUENCE.
%ws = 'base';
ws = 'caller';

% IF A STRUCTURE, THE NAME OF THE FIELD CONTAINING THE DATA.
ff = 'readme';

% DEFINE DELIMITER.
delim = '/';

% PRELIMINARIES.
struct_flag = 0;
if isstruct(readme)
    struct_flag = 1;
    if (isfield(sec,ff))
        tmp = readme; %RENAME
        readme = readme.readme;
    end
end
%else %NOT A STRUCTURE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PART OF CODE.
readme_new = strtok_apl(readme,delim); %SEE ABOVE FOR DEFN OF "delim".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout || ~struct_flag
    varargout{1} = readme_new;
end

if struct_flag
    eval(['tmp.',ff,' = ','readme_new',';'])
    if nargout
        varargout{1} = tmp;
    else
        assignin(ws,'edited',tmp);
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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