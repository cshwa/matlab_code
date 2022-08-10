clc;clear all;close all

grdfile='D:\장기생태\yellow_sea\roms_grd_auto_rdrg2_new8_smooth.nc';
data_info = ncinfo(grdfile, 'lon_rho'); 

styear=2001;
endyear=2010;

outfileu='auto_ERA5_clim01-10_Uwind.nc';
outfilev='auto_ERA5_clim01-10_Vwind.nc';

ncid = netcdf.create(outfileu,'CLOBBER');

eta_rho_dimid = netcdf.defDim(ncid,'eta_rho',data_info.Size(2));
xi_rho_dimid = netcdf.defDim(ncid, 'xi_rho', data_info.Size(1));
eta_u_dimid = netcdf.defDim(ncid,'eta_u',data_info.Size(2));
xi_u_dimid = netcdf.defDim(ncid, 'xi_u', data_info.Size(1)-1);
eta_v_dimid = netcdf.defDim(ncid,'eta_v',data_info.Size(2)-1);
xi_v_dimid = netcdf.defDim(ncid, 'xi_v', data_info.Size(1));
eta_psi_dimid = netcdf.defDim(ncid,'eta_psi',data_info.Size(2)-1);
xi_psi_dimid = netcdf.defDim(ncid, 'xi_psi', data_info.Size(1)-1);
time_dimid = netcdf.defDim(ncid,'Uwind_time', 1460);

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
netcdf.putAtt(ncid,timevarid,'long_name','cyclic day');
netcdf.putAtt(ncid,timevarid,'units','DAYS');
netcdf.putAtt(ncid,timevarid,'cycle_length', 365);

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
time_dimid = netcdf.defDim(ncid2, 'Vwind_time', 1460);

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
netcdf.putAtt(ncid2,timevarid2,'long_name','cyclic day');
netcdf.putAtt(ncid2,timevarid2,'units','DAYS');
netcdf.putAtt(ncid2,timevarid2,'cycle_length', 365);

dvarid2=netcdf.defVar(ncid2,'Vwind', 'NC_FLOAT', [xi_rho_dimid eta_rho_dimid time_dimid]);  %% [x y t]
netcdf.putAtt(ncid2,dvarid2,'long_name','ECMWF ERA5 10 Meter meridional wind velocity');
netcdf.putAtt(ncid2,dvarid2,'units','meter second-1');
netcdf.putAtt(ncid2,dvarid2,'time','Vwind_time');

netcdf.endDef(ncid2);



%meanwd=atand(mean(sind(wd))/mean(cosd(wd)));
%meanwd=mod(meanwd,360)
wtime=[6/24:6/24:365];
yyc=0;
wm=zeros(386,488,10);wd=zeros(386,488,10);
f_u=zeros(386,488,1460);f_v=zeros(386,488,1460);
for tt=1:1:length(wtime)
    for yy=styear:endyear
        yyc=yyc+1;
        filenameu=['auto_ERA5_',num2str(yy),'_Uwind.nc'];
        filenamev=['auto_ERA5_',num2str(yy),'_Vwind.nc'];
        Uwind=ncread(filenameu,'Uwind',[1 1 tt],[Inf Inf 1]);
        Vwind=ncread(filenamev,'Vwind',[1 1 tt],[Inf Inf 1]);
        [twm,twd]=uv2compass(Uwind,Vwind);
        wm(:,:,yyc)=twm;wd(:,:,yyc)=twd;
    end
    f_tu=mean(wm,3).*mean(sin(wd*pi/180),3);
    f_tv=mean(wm,3).*mean(cos(wd*pi/180),3);
    f_u(:,:,tt)=f_tu;f_v(:,:,tt)=f_tv;yyc=0;
end


netcdf.putVar(ncid, timevarid, wtime);
netcdf.putVar(ncid2, timevarid2, wtime);

netcdf.putVar(ncid, dvarid,f_u); %%[x y t], Z'(x,y)
netcdf.putVar(ncid2, dvarid2,f_v); %%[x y t], Z'(x,y)

netcdf.close(ncid);
netcdf.close(ncid2);


