clear all;clc;close all;

%read reference setting file
filepath ='D:/data/avhrr/monthly/';
name1 = 'avhrr_monthly';
file1 = strcat(filepath, 'avhrr_monthly1982_01.nc');
nc = netcdf(file1);
row_lat = nc{'lat'}(:); %unit: degrees_E  
row_lon = nc{'long'}(:);%unit: degrees_N 
row_temp= nc{'temp'}(:);%unit: Celsius
clear nc;
%creating nc file : yearly temp
for i = 2012:2013
    temp_sum = zeros(size(row_temp));
    name2 = num2str(i);
    for j = 1:12
        if j < 10 
            name3 = strcat('0', num2str(j));
        else
            name3 = num2str(j);
        end
        file_name = strcat(filepath,name1, name2, '_',name3,'.nc');
        nc = netcdf(file_name);
        temp_sum = temp_sum + nc{'temp'}(:);
        clear nc;
    end
    temp_mean = temp_sum / 12;
    
    %--------------------creating .nc code-------------------------
    nx = 720 ; ny = 1440;
    filenc=strcat('avhrr_', name2,'.nc');
    ncid = netcdf.create(filenc,'clobber');
     % Crear dimensiones
    dimid_eta_rho = netcdf.defDim(ncid,'eta_rho',nx);
    dimid_xi_rho = netcdf.defDim(ncid,'xi_rho',ny); 

     % Crear variables atributos
    varid_lon = netcdf.defVar(ncid,'long','double',[dimid_xi_rho, dimid_eta_rho]);
    netcdf.putAtt(ncid,varid_lon,'long_name','Longitude of T points')
    netcdf.putAtt(ncid,varid_lon,'units','degrees_E')
     % 
    varid_lat = netcdf.defVar(ncid,'lat','double',[dimid_xi_rho, dimid_eta_rho]);
    netcdf.putAtt(ncid,varid_lat,'long_name','Latitude of T points')
    netcdf.putAtt(ncid,varid_lat,'units','degrees_N')
     % 
    varid_temp = netcdf.defVar(ncid,'temp','double',[dimid_xi_rho, dimid_eta_rho]);
    netcdf.putAtt(ncid,varid_temp,'long_name','Yearly Temperature')
    netcdf.putAtt(ncid,varid_temp,'units','Celsius')
    netcdf.endDef(ncid)

     % % % Agregar datos de coordenadas
    netcdf.putVar(ncid,varid_lon,row_lon');
    netcdf.putVar(ncid,varid_lat,row_lat');
    netcdf.putVar(ncid,varid_temp,temp_mean');
    netcdf.close(ncid) 
end




