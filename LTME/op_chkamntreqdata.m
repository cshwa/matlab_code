function status=op_chkamntreqdata(DSI,R)
status=0;
% Estimate the amount of requested data
Nt=length(R.iTime1:R.iTimeIncr:R.iTime2);
Ny=length(R.iLAT1:R.LATINCR:R.iLAT2);
Nx=length(R.iLON1:R.LONINCR:R.iLON2);
if isfield(R,'iLEVEL1') & isfield(R,'LEVELINCR') & isfield(R,'iLEVEL2')
Nz=length(R.iLEVEL1:R.LEVELINCR:R.iLEVEL2);
else
    Nz=1;
end
Nv=0;
for k=1:length(R.Fields)
    FIELD=R.Fields{k};
    Coords=eval(['DSI.Variables.' FIELD]);   
    if length(Coords)== 3
        Nv=Nv+Nz;
    elseif length(Coords)== 2
        Nv=Nv+1;
    end
end    
MB=round(Nt*Ny*Nx*Nv*8/1024/1024);
MBLIMIT=100;
    if MB>MBLIMIT
ButtonName=questdlg(['The amount of data requested is ' num2str(MB) ' MB.'],...
'Requested Data Amount','Continue','Cancel','Cancel');
        if strcmp(ButtonName,'Cancel') 
            status=1;
            return; 
        end
    end

if Nt==0
    errordlg(['Selected time range contains no data.']); status=1;return; 
end
if Ny==0
    errordlg(['Selected latitude range contains no data.']);status=1;return; 
end
if Nx==0
    errordlg(['Selected longitude range contains no data or corresponds to an eastern and western portion of the requested field(s) with a gap between them. We suggest that you break your request into two separate requests.']);status=1;return;end 
if Nv==0
    errordlg(['No variables selected.']);status=1;return;
end 

