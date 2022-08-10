function surface_forcing_ECMWF(grdfile, year, varname, tinterval,data_daily)
% % This function is based on MATLAB 2017a
    root_dir = '/data1/temp/ECMWF_interim/';
%     filename = strcat('ECMWF_Interim_',varname,'_',num2str(year,'%04i'),'.nc');
%     atmfile  = strcat(root_dir,varname,'/',filename);
switch varname
        case 'airT'
            varname_input = 't2m';
            varname_output = 'Tair';
            long_name = 'ECMWF 2 meter Temperature';
            units = 'Celsius';           
        case 'dewt'
            varname_input = 'd2m';
            varname_output = 'Dair';
            long_name = 'ECMWF 2 meter Dewpoint temperature';
            units = 'Celsius';   
        case 'msl'
            varname_input = varname;
            varname_output = 'Pair';
            long_name = 'ECMWF Mean Sea Level pressure';
            units = 'mbar';   
        case 'ssrd'
            varname_input = varname;
            varname_output = 'swrad';
            long_name = 'ECMWF Short Wave RAdiation Downwards';
            units = 'W/m^2';
        case 'u10'
            varname_input = varname;
            varname_output = 'Uwind';
            long_name = 'ECMWF 10 Meter zonal wind velocity';
            units = 'm/s';
        case 'v10'
            varname_input = varname;
            varname_output = 'Vwind';
            long_name = 'ECMWF 10 Meter meridional wind velocity';
            units = 'm/s';
        otherwise
            '?'
            return
end
    
    varname_time = strcat(varname_output,'_time');
    outfile = strcat(num2str(year,'%04i'),'_',varname_output,'.nc');
    disp(['forcing file is the ',outfile])
    
    totalday  = yeardays(year); 
    time=[0.5:1:totalday-0.5];
    
    lon_rho = ncread(grdfile,'lon_rho');
    lat_rho = ncread(grdfile,'lat_rho');
    data_info = ncinfo(grdfile, 'lon_rho'); 
    
    lon_grd_max = lon_rho(data_info.Size(1),data_info.Size(2));
    lat_grd_max = lat_rho(data_info.Size(1),data_info.Size(2));
    lon_grd_min = lon_rho(1,1);
    lat_grd_min = lat_rho(1,1);

ncid = netcdf.create(outfile,'CLOBBER');

    eta_rho_dimid = netcdf.defDim(ncid,'eta_rho',data_info.Size(2));
    xi_rho_dimid = netcdf.defDim(ncid, 'xi_rho', data_info.Size(1));
    eta_u_dimid = netcdf.defDim(ncid,'eta_u',data_info.Size(2));
    xi_u_dimid = netcdf.defDim(ncid, 'xi_u', data_info.Size(1)-1);
    eta_v_dimid = netcdf.defDim(ncid,'eta_v',data_info.Size(2)-1);
    xi_v_dimid = netcdf.defDim(ncid, 'xi_v', data_info.Size(1));
    eta_psi_dimid = netcdf.defDim(ncid,'eta_psi',data_info.Size(2)-1);
    xi_psi_dimid = netcdf.defDim(ncid, 'xi_psi', data_info.Size(1)-1);
    time_dimid = netcdf.defDim(ncid, varname_time, 0);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' ROMS Surface forcing file ');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', ' Bulk Formular Forcing file ');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', ' ECMWF ERA Interim ');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by Y.Y.Kim');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    timevarid=netcdf.defVar(ncid, varname_time, 'NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid,timevarid,'long_name','cyclic day');
    netcdf.putAtt(ncid,timevarid,'units','DAYS');
    netcdf.putAtt(ncid,timevarid,'cycle_length', totalday);

    dvarid=netcdf.defVar(ncid,varname_output, 'NC_FLOAT', [xi_rho_dimid eta_rho_dimid time_dimid]);  %% [x y t]
    netcdf.putAtt(ncid,dvarid,'long_name',long_name);
    netcdf.putAtt(ncid,dvarid,'units',units);
    netcdf.putAtt(ncid,dvarid,'time',varname_time);

    netcdf.endDef(ncid);

    for i=1:totalday
        netcdf.putVar(ncid, timevarid, i-1, 1, time(i));
        Z=squeeze(data_daily(:,:,i)); %%[x y]
        netcdf.putVar(ncid, dvarid, [0 0 i-1], [data_info.Size(1) data_info.Size(2) 1], Z); %%[x y t], Z'(x,y)
    end

    netcdf.close(ncid);
    status=1;