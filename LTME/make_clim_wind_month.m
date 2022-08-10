clc;clear all;close all

grdfile='D:\장기생태\yellow_sea\roms_grd_auto_rdrg2_new8_smooth.nc';
data_info = ncinfo(grdfile, 'lon_rho'); 

styear=2006;
endyear=2010;

outfileu='auto_ERA5_clim01-10_Uwind_monthly.nc';
outfilev='auto_ERA5_clim01-10_Vwind_monthly.nc';

ncid = netcdf.create(outfileu,'CLOBBER');

eta_rho_dimid = netcdf.defDim(ncid,'eta_rho',data_info.Size(2));
xi_rho_dimid = netcdf.defDim(ncid, 'xi_rho', data_info.Size(1));
eta_u_dimid = netcdf.defDim(ncid,'eta_u',data_info.Size(2));
xi_u_dimid = netcdf.defDim(ncid, 'xi_u', data_info.Size(1)-1);
eta_v_dimid = netcdf.defDim(ncid,'eta_v',data_info.Size(2)-1);
xi_v_dimid = netcdf.defDim(ncid, 'xi_v', data_info.Size(1));
eta_psi_dimid = netcdf.defDim(ncid,'eta_psi',data_info.Size(2)-1);
xi_psi_dimid = netcdf.defDim(ncid, 'xi_psi', data_info.Size(1)-1);
time_dimid = netcdf.defDim(ncid,'Uwind_time', 4);

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' ROMS Surface forcing file ');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'title', ' Bulk Formular Forcing file ');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'source', ' ECMWF ERA5 Clim');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by Y.J.Tak');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

timevarid=netcdf.defVar(ncid, 'Uwind_time', 'NC_DOUBLE', time_dimid);
netcdf.putAtt(ncid,timevarid,'long_name','cyclic month');
netcdf.putAtt(ncid,timevarid,'units','MONTHS');
netcdf.putAtt(ncid,timevarid,'cycle_length', 4);

dvarid=netcdf.defVar(ncid,'Uwind', 'NC_FLOAT', [xi_rho_dimid eta_rho_dimid time_dimid]);  %% [x y t]
netcdf.putAtt(ncid,dvarid,'long_name','ECMWF ERA5 10 Meter zonal wind velocity');
netcdf.putAtt(ncid,dvarid,'units','meter second-1');
netcdf.putAtt(ncid,dvarid,'time','Uwind_time');

netcdf.endDef(ncid);


ncid2 = netcdf.create(outfilev,'CLOBBER');

eta_rho_dimid = netcdf.defDim(ncid2,'eta_rho',data_info.Size(2));
xi_rho_dimid = netcdf.defDim(ncid2, 'xi_rho', data_info.Size(1));
eta_u_dimid = netcdf.defDim(ncid2,'eta_u',data_info.Size(2));
xi_u_dimid = netcdf.defDim(ncid2, 'xi_u', data_info.Size(1)-1);
eta_v_dimid = netcdf.defDim(ncid2,'eta_v',data_info.Size(2)-1);
xi_v_dimid = netcdf.defDim(ncid2, 'xi_v', data_info.Size(1));
eta_psi_dimid = netcdf.defDim(ncid2,'eta_psi',data_info.Size(2)-1);
xi_psi_dimid = netcdf.defDim(ncid2, 'xi_psi', data_info.Size(1)-1);
time_dimid = netcdf.defDim(ncid2, 'Vwind_time', 4);

netcdf.putAtt(ncid2,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' ROMS Surface forcing file ');
netcdf.putAtt(ncid2,netcdf.getConstant('NC_GLOBAL'), ...
    'title', ' Bulk Formular Forcing file ');
netcdf.putAtt(ncid2,netcdf.getConstant('NC_GLOBAL'), ...
    'source', ' ECMWF ERA5 Clim');
netcdf.putAtt(ncid2,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by Y.J.Tak');
netcdf.putAtt(ncid2,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

timevarid2=netcdf.defVar(ncid2, 'Vwind_time', 'NC_DOUBLE', time_dimid);
netcdf.putAtt(ncid2,timevarid2,'long_name','cyclic month');
netcdf.putAtt(ncid2,timevarid2,'units','MONTHS');
netcdf.putAtt(ncid2,timevarid2,'cycle_length', 4);

dvarid2=netcdf.defVar(ncid2,'Vwind', 'NC_FLOAT', [xi_rho_dimid eta_rho_dimid time_dimid]);  %% [x y t]
netcdf.putAtt(ncid2,dvarid2,'long_name','ECMWF ERA5 10 Meter meridional wind velocity');
netcdf.putAtt(ncid2,dvarid2,'units','meter second-1');
netcdf.putAtt(ncid2,dvarid2,'time','Vwind_time');

netcdf.endDef(ncid2);



%meanwd=atand(mean(sind(wd))/mean(cosd(wd)));
%meanwd=mod(meanwd,360)
 t1 = datetime(1981,01,01,00,00,00); % >> 21,00,00 (GMT)
 t2 = datetime(1981,12,31,23,00,00); %>> 22,00,00  (GMT)
 t = t1:hours(6):t2
 time = datenum(t)
 time_c = datestr(time,'mm-dd hh'); 


 
 matcht{1}={'03-01 00'};
 matcht{2}={'05-31 18'};
 
 matcht{3}={'06-01 00'};
 matcht{4}={'08-31 18'};
 
 matcht{5}={'09-01 00'};
 matcht{6}={'11-30 18'};
 
 matcht{7}={'01-01 00'};
 matcht{8}={'02-28 18'};
 matcht{9}={'12-01 00'};
 matcht{10}={'12-31 18'};
 
 for i = 1:length(matcht)
idx=strcmp(matcht{i},time_c);
idxx(i)= find(idx ~= 0);
 end
 

wtime=[6/24:6/24:365];
yyc=0;
% wm=zeros(386,488);wd=zeros(386,488);
f_u=zeros(386,488,4);f_v=zeros(386,488,4);
for tt=1:1:length(wtime)
        yyc=yyc+1;
        filenameu=['auto_ERA5_clim01-10_Uwind.nc'];
        filenamev=['auto_ERA5_clim01-10_Vwind.nc'];
        Uwind=ncread(filenameu,'Uwind',[1 1 tt],[Inf Inf 1]);
        Vwind=ncread(filenamev,'Vwind',[1 1 tt],[Inf Inf 1]);
        [twm,twd]=uv2compass(Uwind,Vwind);
        wm(:,:,yyc)=twm;wd(:,:,yyc)=twd;
end

ti=0
for i = 1:2:6
    ti=ti+1;
    f_tu=mean(wm(:,:,idxx(i):idxx(i+1)),3).*mean(sin(wd(:,:,idxx(i):idxx(i+1))*pi/180),3);
    f_tv=mean(wm(:,:,idxx(i:i+1)),3).*mean(cos(wd(:,:,idxx(i:i+1))*pi/180),3);
    f_u(:,:,ti)=f_tu;f_v(:,:,ti)=f_tv;yyc=0;
end
%winter (DJF)
f_tu=mean(wm(:,:,[idxx(7):idxx(8) idxx(9):idxx(10)]),3).*mean(sin(wd(:,:,[idxx(7):idxx(8) idxx(9):idxx(10)])*pi/180),3);
f_tv=mean(wm(:,:,[idxx(7):idxx(8) idxx(9):idxx(10)]),3).*mean(cos(wd(:,:,[idxx(7):idxx(8) idxx(9):idxx(10)])*pi/180),3);
f_u(:,:,4)=f_tu;f_v(:,:,4)=f_tv;
 
 
 

% netcdf.putVar(ncid, timevarid, wtime);
% netcdf.putVar(ncid2, timevarid2, wtime);

netcdf.putVar(ncid, dvarid,f_u); %%[x y t], Z'(x,y)
netcdf.putVar(ncid2, dvarid2,f_v); %%[x y t], Z'(x,y)

netcdf.close(ncid);
netcdf.close(ncid2);


%% test
u=ncread('auto_ERA5_clim01-10_Uwind_monthly.nc','Uwind');
v=ncread('auto_ERA5_clim01-10_Vwind_monthly.nc','Vwind');

quiver(f_u(1:10:end,1:10:end,1),f_v(1:10:end,1:10:end,1))
quiver(u(1:10:end,1:10:end,2),v(1:10:end,1:10:end,2))
quiver(u(1:10:end,1:10:end,3),v(1:10:end,1:10:end,3))
quiver(u(1:10:end,1:10:end,4),v(1:10:end,1:10:end,4))

