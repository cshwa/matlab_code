function status=op_saved(D,R,kf,kf0,kw0)
status=0;
if R.SaveWorkspace=='y' & strcmp(R.SaveMode,'opendap')==1
KF=num2str(kf+kw0,'%04d');
WSN=['opendap_' KF];
assignin('base',WSN,D)
end
if R.SaveFiles == 'y' & strcmp(R.SaveMode,'opendap')==1
% save extracted data to files
KF=num2str(kf+kf0,'%04d');
FN3=[R.FNPREFIX 'opendap_' KF '.mat'];
    if length(R.DIRNAME) > 0
    FN3=[R.DIRNAME '/' FN3];
    else
    FN3=[pwd '/' FN3];
    end
FileName=FN3;
disp(['    saving data to file: ' FN3])
save(FN3,'D')
end
if R.SaveFiles == 'y' & strcmp(R.SaveMode,'netcdf')==1
% save extracted data to files
KF=num2str(kf+kf0,'%04d');
FN3=[R.FNPREFIX 'opendap_' KF '.nc'];
    if length(R.DIRNAME) > 0
    FN3=[R.DIRNAME '/' FN3];
    else
    FN3=[pwd '/' FN3];
    end
FileName=FN3;
disp(['    saving data to file: ' FN3])
a=version;a=a(1:3);a=str2num(a);
% for earlier than MATLAB 2008b (Ver 7.7) use third party mexnc package
    if a < 7.7
    op_writenetcdf3(D,FN3);
    else % otherwise use internal netcdf support 
    op_writenetcdf4(D,FN3);
    end
end

DSN=D.Attributes.DataSetName;
DATE=datestr(datenum(D.user_friendly_time,'yyyy-mm-dd HH:MM:SS'),'yyyymmddHHMMSS');
if R.SaveWorkspace=='y' & strcmp(R.SaveMode,'native')==1
WSN=[DSN '_' DATE];
assignin('base',WSN,D)
end
if R.SaveFiles == 'y' & strcmp(R.SaveMode,'native')==1
% save extracted data to files
FN3=[R.FNPREFIX DSN '_' DATE '.mat'];
    if length(R.DIRNAME) > 0
    FN3=[R.DIRNAME '/' FN3];
    else
    FN3=[pwd '/' FN3];
    end
FileName=FN3;
disp(['    saving data to file: ' FN3])
save(FN3,'D')

end


