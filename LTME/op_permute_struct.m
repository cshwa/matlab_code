function B = op_permute_struct(A,fields,order)
% Permutes the dimensions of a structure.
%
% Acts on fieldnames "fields" and uses the order found in "order".
% Resembles Matlab's "permute".
%
% B=op_permute_struct(A,fields,order);
%
% INPUT(S):
%
%   A -- input structure
%   fields -- cell array of strings, usually obtained by using
%             "fieldnames" function.
%   order -- vector specifying order of desired result
%            e.g., order = [1 2 4 3]
%
% OUTPUT(S):
%
%   B -- output structure
%
% OPeNDAP Science Team
% Copyright 2007, 2008
% $Version 2.0.0$

%==========================================================================
% Christian Buckingham
%
% REVISION HISTORY:
% 2007/05/01 0.0.x created, ceb
% 2008/03/01 2.0.0 released
%==========================================================================

B=A;
%fields=fieldnames(A);
len=length(fields);
for ii=1:len
    tmp_field=fields{ii};
    tmp_data=getfield(A,tmp_field);
    tmp_data=permute(tmp_data,order);
    B=setfield(B,tmp_field,tmp_data);
end

return %endoffunction