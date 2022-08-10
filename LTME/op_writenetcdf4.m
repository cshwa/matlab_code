function op_writenetcdf(D,FNC)
% converts from opendap_nnnn data frame to
% NetCDF file using MATLAB internal netcdf functions
% D structure is a data frame
% FNC output NetCDF file name

% 2009-07-13 converted from mexnc to MATLAB internal netcdf functions

%MV=1.2677e+030; % replacing NaN is not required

%[ncid, status]=mexnc('CREATE',FNC,'clobber');
ncid = netcdf.create(FNC,'clobber');

% D consists of (fields, coordinates), dimensions, variables and attributes

% define dimensions
Dimensions=fieldnames(D.Dimensions);
for k=1:length(Dimensions)
DIM=Dimensions{k};
%[dimid(k),status]=mexnc('DEF_DIM',ncid,DIM,eval(['D.Dimensions.' DIM]));
dimid(k) = netcdf.defDim(ncid,DIM,eval(['D.Dimensions.' DIM]));
%[dimid_longitude,   status]=mexnc('DEF_DIM',ncid,'longitude',   length(D.longitude))
%[dimid_depth,       status]=mexnc('DEF_DIM',ncid,'depth',       length(D.depth))
%[dimid_time,        status]=mexnc('DEF_DIM',ncid,'time',        length(D.time))
end

Variables=fieldnames(D.Variables);
%define variables
for k=1:length(Variables)
VAR=Variables{k};
VarDims=eval(['D.Variables.' VAR]);
    List=0;
    for kk=1:length(VarDims)
    List(kk)=strmatch(VarDims(kk),Dimensions);
    end
%[varid(k),status]=mexnc('DEF_VAR',ncid,VAR,'DOUBLE',length(List),dimid(List));
varid(k) = netcdf.defVar(ncid,VAR,'double',dimid(List));
end

% add attributes
for k=1:length(Variables)
VAR=Variables{k};
Atts=eval(['fieldnames(D.Attributes.' VAR ');']);
    for kk=1:length(Atts)
        ATT=Atts{kk};
        Att=eval(['D.Attributes.' VAR '.' ATT ';']);
        if ischar(Att);     
%status=mexnc('PUT_ATT_TEXT',ncid,varid(k),ATT,'CHAR',length(Att),Att);
netcdf.putAtt(ncid,varid(k),ATT,Att);
        end
        if isnumeric(Att);
%status=mexnc('PUT_ATT_DOUBLE',ncid,varid(k),ATT,'DOUBLE',length(Att),Att);
netcdf.putAtt(ncid,varid(k),ATT,Att);
        end
    end
end


% Global Attributes
ATT='DataSetName';Att=eval(['D.Attributes.' ATT ';']);
%status=mexnc('PUT_ATT_TEXT',ncid,-1,ATT,'char',length(Att),Att);
netcdf.putAtt(ncid,-1,ATT,Att);
ATT='Readme';Att=eval(['D.Attributes.' ATT ';']);
%status=mexnc('PUT_ATT_TEXT',ncid,-1,ATT,'char',length(Att),Att);
netcdf.putAtt(ncid,-1,ATT,Att);
ATT='DataSetBranch';Att=eval(['D.Attributes.' ATT ';']);
%status=mexnc('PUT_ATT_TEXT',ncid,-1,ATT,'char',length(Att),Att);
netcdf.putAtt(ncid,-1,ATT,Att);
ATT='RequestDate';Att=eval(['D.Attributes.' ATT ';']);
%status=mexnc('PUT_ATT_TEXT',ncid,-1,ATT,'char',length(Att),Att);
netcdf.putAtt(ncid,-1,ATT,Att);

% end define mode
%status=mexnc('ENDDEF',ncid);
netcdf.endDef(ncid);

% add data

for k=1:length(Variables)
VAR=Variables{k};
% dimensions in reverse order
A=eval(['D.' VAR]);
%ndims(A);
%size(A);
A=permute(A,[ndims(A):-1:1]);
%size(A);
%A(find(isnan(A)))=MV; % NaN to Missing Values
%status=mexnc('PUT_VAR_DOUBLE',ncid,varid(k),A);
netcdf.putVar(ncid,varid(k),A);
end

%status=mexnc('CLOSE',ncid);
netcdf.close(ncid);