function D=op_initd0(DSI,R)
clear D
% create standard blank opendap_nnnn data structure
NV=length(R.Variables);
for kv=1:NV
FIELD=R.Variables{kv}; %
eval(['D.' FIELD '=NaN;']);
end
% R.Variables include lat,lon,time  

D.user_friendly_time='';

D.Dimensions=DSI.Dimensions;
D.Dimensions.time=1;

for kv=1:NV
FIELD=R.Variables{kv}; %
eval(['D.Variables.'  FIELD '=DSI.Variables.'  FIELD ';']);
eval(['D.Attributes.' FIELD '=DSI.Attributes.' FIELD ';']);
end

%Global Attributes
D.Attributes.DataSetName=DSI.DataSetName;
D.Attributes.DataSetBranch=DSI.DataSetBranch;
D.Attributes.Readme=DSI.Readme;
D.Attributes.RequestDate=R.RequestDate;
D.Attributes.Request=R;
