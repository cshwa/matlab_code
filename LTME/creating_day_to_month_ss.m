% origin by Young-Ki Kim reveced 2017.06.14
% edited by E.B. Cho 2017.06.14

clear all;clc;close all;

%read reference setting file
yr = 2013;
filepath1 = ['E:\11.사업\장기생태_3단계\Data\위성자료\daily\' num2str(yr) '\'];%daily folder
filepath2 ='E:\11.사업\장기생태_3단계\Data\위성자료\monthly\';   %monthly folder
name1 = ['avhrr-only-v2.' num2str(yr)];
name11 = ['avhrr_monthly' num2str(yr)];
file1 = strcat(filepath2, 'avhrr_monthly1982_01.nc');
nc = netcdf(file1);
row_lat = nc{'lat'}(:); %unit: degrees_E  
row_lon = nc{'long'}(:);%unit: degrees_N 
row_temp= nc{'temp'}(:);
clear nc;

%creating nc file : yearly temp
for i = 1:12
    temp_sum = zeros(size(row_temp));
    %define name1(month name)
    if i < 10 
        name2 = strcat('0', num2str(i));
    else
        name2 = num2str(i);
    end
    %define mx
    if i == 2;mx = 28;
    elseif i==4 || 6 || 9 || 11;mx=30;
    elseif i==1 || 3 || 5 || 7 || 8 ||10||12; mx=31;
    end
    for j = 1:mx
        if j <10
            name3 = strcat('0', num2str(j));
        else
            name3 = num2str(j);
        end
        file_name = strcat(filepath1,name1, name2,name3,'.nc');
        nc = netcdf(file_name);
        %make sst data
        sst = nc{'sst'}(1,1,:,:);
        sst(sst==-999) = nan;  %missing process
        sst = 0.01 * sst + 0 ; %scale,add factor 
        temp_sum = temp_sum + sst;
        clear nc;
    end
    temp_mean = temp_sum / mx;
    
    %--------------------creating .nc code-------------------------
    nx = 720 ; ny = 1440;
%     filenc=strcat(name11,'_',name2, '.nc');
 filenc=strcat(filepath2,name11,'_',name2, '.nc');
    ncid = netcdf.create(filenc,'NC_WRITE');
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




