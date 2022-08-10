function D=op_initd(DSI,R,kT)
clear D
% create standard blank opendap_nnnn data structure
NV=length(R.Variables);
for kv=1:NV
FIELD=R.Variables{kv}; %
eval(['D.' FIELD '=NaN;']);
end

% lat,lon,time are supposed to be the same for all fields  

if ndims(DSI.Latitude)==1
D.latitude=DSI.Latitude(R.iLAT1:R.LATINCR:R.iLAT2);
elseif ndims(DSI.Latitude)==2
D.latitude=DSI.Latitude(R.iLAT1:R.LATINCR:R.iLAT2,R.iLON1:R.LONINCR:R.iLON2);
else
    stop('more than 2 latitude dimensions not supported')
end
if ndims(DSI.Longitude)==1
D.longitude=DSI.Longitude(R.iLON1:R.LONINCR:R.iLON2);
elseif ndims(DSI.Longitude)==2
D.longitude=DSI.Longitude(R.iLAT1:R.LATINCR:R.iLAT2,R.iLON1:R.LONINCR:R.iLON2);
else
    stop('more than 2 latitude dimensions not supported')
end

if ~isempty(strmatch('depth',R.Coordinates))
D.depth=DSI.Depth(R.iLEVEL1:R.LEVELINCR:R.iLEVEL2);
end

%D.time=Time(k); % datenum;  
D.time=(DSI.Time(kT)-datenum(1970,1,1,0,0,0))*86400;% in seconds since 1970-01-01 00:00:00
D.user_friendly_time='';

D.Dimensions=DSI.Dimensions;

for kv=1:NV
FIELD=R.Variables{kv}; %
eval(['D.Variables.' FIELD '=DSI.Variables.' FIELD ';']);
eval(['D.Attributes.' FIELD '=DSI.Attributes.' FIELD ';']);
end

%Global Attributes
D.Attributes.DataSetName=DSI.DataSetName;
D.Attributes.DataSetBranch=DSI.DataSetBranch;
D.Attributes.RequestDate=R.RequestDate;
D.Attributes.Request=R;
