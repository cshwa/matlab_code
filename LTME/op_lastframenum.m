function [kw,kf]=op_lastframenum(R)
% find the latest saved opendap frame number 
%KF=num2str(kf,'%04d');
% the latest workspace frame number of form opendap_nnnn
a=evalin('base','who');k=strmatch('opendap_',a);
if isempty(k) k=0; 
else 
a=a{k(end)};
ki=strfind(a,'_');
    if length(ki)==1 ki(2)=length(a)+1; end
k=str2num(a(ki(1)+1:ki(2)-1));
end
if isempty(k) k=0;
end
kw=k;
% the latest saved file of form FNPREFIXnnnn.mat or .nc
EXT='mat'; 
if strcmp(R.SaveMode,'netcdf')==1
    EXT='nc';
end
if isempty(R.DIRNAME) R.DIRNAME='.'; end
%[R.DIRNAME '\' R.FNPREFIX '*.mat']
a=dir([R.DIRNAME '\' R.FNPREFIX '*.' EXT]);
a={a.name};
if isempty(a) k=0;
else
a=a{end};
k=str2num(a(length(R.FNPREFIX)+1:end-length(EXT)-1));
end
if isempty(k) k=0;
end
kf=k;