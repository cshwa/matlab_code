function status=op_checkrequest(DSI,R)
status=0;
if length(R.Fields)==0
    errordlg(['No variables selected.']);status=1;return;
end
if(isempty(R.DATE1) | isempty(R.DATE2)) | (str2num(R.DATE2)<str2num(R.DATE1))
    errordlg(['Selected dates are invalid.']);status=1;return
end
if(isempty(R.LAT1) | isempty(R.LAT2)) | (R.LAT2<R.LAT1)
    errordlg(['Selected latitudes are invalid.']);status=1;return
end
if(isempty(R.LON1) | isempty(R.LON2)) % LON2 can be less than LON1
    errordlg(['Selected longitudes are invalid.']);status=1;return
end
if(isempty(R.DATEINCR) | isempty(R.LATINCR) | isempty(R.LONINCR) |R.LATINCR==0|R.LONINCR==0|R.DATEINCR==0) 
    errordlg(['Subsampling values are invalid.']);status=1;return
end
if R.SaveWorkspace=='n' & R.SaveFiles=='n'
errordlg('No saving data in workspace or to files selected.'); status=1;return
end

