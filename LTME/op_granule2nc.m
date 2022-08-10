function op_granule2nc(D,FNC)
% converts from opendap_nnnn data granule to
% NetCDF file using MEXNC
% D structure is a data granule
% FNC output NetCDF file name

%MV=1.2677e+030; % replacing NaN is not required

% D consists of fields, coordinates, and metadata
fldnms=fieldnames(D);
k1=find(strcmp(fldnms,'latitude'));
k2=find(strcmp(fldnms,'longitude'));
k3=find(strcmp(fldnms,'time'));
k4=find(strcmp(fldnms,'depth'));
kc=min([k1,k2,k3,k4]);
km=find(strcmp(fldnms,'metadata'));
Fields=fldnms(1:kc-1);
Coords=fldnms(kc:km-1);

[ncid, status]=mexnc('CREATE',FNC,'clobber');

% define dimensions
for k=1:length(Coords)
COORD=Coords{k};
[dimid(k),status]=mexnc('DEF_DIM',ncid,COORD,length(eval(['D.' COORD])));
%[dimid_longitude,   status]=mexnc('DEF_DIM',ncid,'longitude',   length(D.longitude))
%[dimid_depth,       status]=mexnc('DEF_DIM',ncid,'depth',       length(D.depth))
%[dimid_time,        status]=mexnc('DEF_DIM',ncid,'time',        length(D.time))
end

% define coordinate variables 
for k=1:length(Coords)
COORD=Coords{k};
[coordid(k),status]=mexnc('DEF_VAR',ncid,COORD,'DOUBLE',1,dimid(k));
%[varid_longitude,   status]=mexnc('DEF_VAR',ncid,'longitude','DOUBLE',1,dimid(2))
%[varid_depth,       status]=mexnc('DEF_VAR',ncid,'depth','DOUBLE',1,dimid(3))
%[varid_time,        status]=mexnc('DEF_VAR',ncid,'time','DOUBLE',1,dimid(4))
end

%define field variables
for k=1:length(Fields)
FIELD=Fields{k};
Field_Coords=eval(['D.metadata.' FIELD '_variables']);
    for kk=1:length(Field_Coords)
    List(kk)=strmatch(Field_Coords(kk),Coords);
    end
[fieldid(k),status]=mexnc('DEF_VAR',ncid,FIELD,'DOUBLE',length(List),dimid(List));
end

% add attributes
for k=1:length(Fields)
FIELD=Fields{k};
    if isstruct(D.metadata.url)==1
        if isfield(D.metadata.url,FIELD)
S=eval(['D.metadata.url.' FIELD]);
        else
            S=''; % example: sst_masked does not have url as sst
        end
    else
        S=D.metadata.url;
    end
status=mexnc('PUT_ATT_TEXT',ncid,fieldid(k),'url','CHAR',length(S),S);
%status=mexnc('PUT_ATT_TEXT',ncid,varid_longitude,'units','char',length('degrees'),'degrees')
%status=mexnc('PUT_ATT_TEXT',ncid,varid_depth,'units','char',length('meters'),'meters')
%status=mexnc('PUT_ATT_TEXT',ncid,varid_time,'units','char',length('seconds since 1970-01-01 00:00:00'),'seconds since 1970-01-01 00:00:00')
%status=mexnc('PUT_ATT_TEXT',ncid,varid_salinity,'units','char',length('PSU'),'PSU')
% replacing NaN is not required
%status=mexnc('PUT_ATT_DOUBLE',ncid,fieldid(k),'missing_value','DOUBLE',1,MV);
end

% add time attributes
k=strmatch('time',Coords);
    if isfield(D.metadata,'temporal')
S=D.metadata.temporal;
status=mexnc('PUT_ATT_TEXT',ncid,coordid(k),'temporal','CHAR',length(S),S);
    end
    if isfield(D.metadata,'reference_time')
S=D.metadata.reference_time;
status=mexnc('PUT_ATT_TEXT',ncid,coordid(k),'reference_time','CHAR',length(S),S);
    end
    if isfield(D.metadata,'user_friendly_time')
S=D.metadata.user_friendly_time;
status=mexnc('PUT_ATT_TEXT',ncid,coordid(k),'user_friendly_time','CHAR',length(S),S);
    end

% Global Attributes
S=D.metadata.name;
status=mexnc('PUT_ATT_TEXT',ncid,-1,'name','char',length(S),S);
S=D.metadata.request_date;
status=mexnc('PUT_ATT_TEXT',ncid,-1,'request_date','char',length(S),S);

% end define mode
status=mexnc('ENDDEF',ncid);

% add data
for k=1:length(Coords)
COORD=Coords{k};
A=eval(['D.' COORD]);
status=mexnc('PUT_VAR_DOUBLE',ncid,coordid(k),A);
end

for k=1:length(Fields)
FIELD=Fields{k};
% dimensions in reverse order
A=eval(['D.' FIELD]);
%ndims(A);
%size(A);
A=permute(A,[ndims(A):-1:1]);
%size(A);
%A(find(isnan(A)))=MV; % NaN to Missing Values
status=mexnc('PUT_VAR_DOUBLE',ncid,fieldid(k),A);
end

status=mexnc('CLOSE',ncid);
