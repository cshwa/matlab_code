function varargout=op_orderfields(varargin)
% Orders the fieldnames of a structure so that they follow a standard format.
%
% S2 = op_orderfields(S1,token_fields)
%
% INPUT(S):
%
%   S1 -- input structure
%   token_fields -- (optional) cell array of fieldnames. If not input
%       the program uses a pre-defined set of token_fields. Enter
%       "type op_orderfields" to see the list.
%
% OUTPUT(S):
%
%   S2 -- (optional) output structure with fieldnames re-ordered.
%
% OPeNDAP Science Team
% Copyright 2007, 2008
% $Version 2.0.0$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2006/12/01 0.0.x created, ceb
% 2007/04/12 0.0.x edited, ceb
% 2008/03/01 2.0.0 released
%==========================================================================

% Define workspace.
%ws = 'caller';

% Inputs.
if (nargin == 1)
    S1 = varargin{1};

% DEFINE INDEPENDENT FIELDS.
ii=1;
indep_fields{ii}='latitude'; ii=ii+1;
indep_fields{ii}='longitude'; ii=ii+1;
indep_fields{ii}='time'; ii=ii+1;
indep_fields{ii}='depth'; ii=ii+1;

% DEFINE METADATA FIELD.
meta_field='metadata';

% DEFINE TOKEN FIELDS (UNION OF INDEPENDENT FIELDS AND METADATA FIELD).
token_fields = indep_fields; ii = length(indep_fields)+1;
token_fields{ii}=meta_field; ii=ii+1;

% DEFINE OPENDAP-SPECIFIC FIELDS (WITHIN METADATA).
% THESE WILL BE CELL ARRAYS.
ii=1;
opend_fields{ii}='url'; ii=ii+1;
opend_fields{ii}='attributes'; ii=ii+1;
opend_fields{ii}='name'; ii=ii+1;
opend_fields{ii}='temporal'; ii=ii+1;
opend_fields{ii}='readme'; ii=ii+1;
opend_fields{ii}='timeGregorian'; ii=ii+1;
opend_fields{ii}='request_date'; ii=ii+1;
opend_fields{ii}='request'; ii=ii+1;

elseif (nargin == 2)
    S1 = varargin{1};
    token_fields = varargin{2};
else
    %error
end

% Field of structure.
tmp1 = fieldnames(S1);
%tmp1 = evalin(ws,['fieldnames(',name,');']);

% Fields not contained in token_fields.
rrr = setdiff(lower(tmp1),token_fields);
[C,IA,IB] = intersect(lower(tmp1),token_fields);
token_fields(IB) = tmp1(IA);

% Sort the above fields.
rrr = sort(rrr);

% Combination of fields in token_fields with above.
N = length(rrr);
M = length(token_fields);

if (M + N ~= length(tmp1))
    error('Problem: input structure does not contain token_fields.');
end
tmp2 = cell(1,N+M);
for ii = 1:M
    tmp2{end - M + ii} = token_fields{ii};
end
for ii = 1:N
    tmp2{ii} = rrr{ii};
end
tmp2 = tmp2';

% Make structure's fields ordered same as token_fields.

if (nargout == 0)
    %evalin(ws,[name,' = orderfields(',name,',',;']);
else
    varargout{1} = orderfields(S1,tmp2);
end

return %endoffunction