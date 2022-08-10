% op_multiplotsst (script) - Plot first 6 SST images to one plot
%
% OPeNDAP Science Team
% Copyright 2008
% $Version 2.0.0$

%==========================================================================
% Peter Cornillon
%
% REVISION HISTORY:
% 2006/xx/xx 0.0.0 created, pcornillon
% 2008/03/01 2.0.0 released
%==========================================================================

VarNames = whos('sst_*');
NumVars = length(VarNames);

figure(1)
for i=1:NumVars
    subplot(3,2,i) 
    eval([ 'PlotImage(1,' VarNames(i).name ',1);'])
end

figure(2)
for i=1:NumVars
    subplot(3,2,i) 
    eval([ 'hist(' VarNames(i).name '.sst(:),100);'])
end


clear VarNames NumVars