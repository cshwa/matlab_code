function[R,D,status]=op_coord2ind(DSI,R,D);
status=0;
if DSI.Latitude(end)>DSI.Latitude(1)
R.iLAT1=min(find(DSI.Latitude >= R.LAT1));
R.iLAT2=max(find(DSI.Latitude <= R.LAT2));
else
% Find indices within the specified range of lats and lons
% Note that array Latitude runs from North to South
R.iLAT1=min(find(DSI.Latitude <= R.LAT2));
R.iLAT2=max(find(DSI.Latitude >= R.LAT1));
end

% adjust entered longitudes modulo 360 to match the dataset range 
LonM=0.5*(DSI.Longitude(1)+DSI.Longitude(end));
R.LON1=R.LON1+floor((LonM+180-R.LON1)/360)*360;
R.LON2=R.LON2+floor((LonM+180-R.LON2)/360)*360;
% the above might not work if Longitude array is not monotonic
% make sure that DSI takes care of that

% Note that Longitude may jump from 180E to -179W in Western Pacific Region
%i=find(DSI.Longitude >= R.LON1);[dummy,j]=min(DSI.Longitude(i));iLON1=i(j);iLON1=iLON1-1;
%i=find(DSI.Longitude <= R.LON2);[dummy,j]=max(DSI.Longitude(i));iLON2=i(j);iLON2=iLON2-1;
R.iLON1=min(find(DSI.Longitude >= R.LON1));
R.iLON2=max(find(DSI.Longitude <= R.LON2));

R.CLAT=['[' num2str(R.iLAT1-1) ':' num2str(R.LATINCR) ':' num2str(R.iLAT2-1) ']']; % Constraint Latitude
R.CLON=['[' num2str(R.iLON1-1) ':' num2str(R.LONINCR) ':' num2str(R.iLON2-1) ']']; % Constraint Longitude

JD1=datenum(R.DATE1,'yyyymmdd');
JD2=datenum(R.DATE2,'yyyymmdd')+1;
%assignin('base','Time',Time)
R.iTime1=min(find(DSI.Time >= JD1));
R.iTime2=max(find(DSI.Time <  JD2));

D.latitude=DSI.Latitude(R.iLAT1:R.LATINCR:R.iLAT2);
D.Dimensions.latitude=length(D.latitude);
D.longitude=DSI.Longitude(R.iLON1:R.LONINCR:R.iLON2);
D.Dimensions.longitude=length(D.longitude);