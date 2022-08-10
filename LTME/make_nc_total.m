
cd '/data1/cshwa/tmp/history'
close all;clear all; clc;

minyear=1976;
maxyear=2005;
list = dir('*hur_interped_nor_*');

hur_sd_1 = load(list(2).name,'hur_sd');
hur_sd_2 = load(list(3).name,'hur_sd');
hur_sd_3 = load(list(4).name,'hur_sd');
hur_sd_4 = load(list(5).name,'hur_sd');
hur_sd_1= hur_sd_1.hur_sd;
hur_sd_2= hur_sd_2.hur_sd;
hur_sd_3= hur_sd_3.hur_sd;
hur_sd_4= hur_sd_4.hur_sd;
load('standard_grid.mat');

[size11 size12 size13]= size(hur_sd_1);
[size21 size22 size23]= size(hur_sd_2);
[size31 size32 size33]= size(hur_sd_3);
[size41 size42 size43]= size(hur_sd_4);

%%% concatenate array along specified dimension
hur_sd = cat(3, hur_sd_1, hur_sd_2, hur_sd_3, hur_sd_4); %combine matrix into one

t_sd=size(hur_sd,3)

timetime=minyear:1/365:maxyear+1-1/365;

tstep=0;
leap_yr = [1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004];
 for j = minyear:maxyear
    tstep =tstep +1;
    fstep=(tstep-1)*365+1:1:tstep*365;
    ftime=timetime(fstep); % slice every year in time
    var_temp=hur_sd(:,:,fstep); % slice every year in hur
    
    
    if find(leap_yr == j)
    t_temp(1:58) = ftime(1:58);
    t_temp(59) = ftime(58);
    t_temp(60:366) = ftime(59:end);
    
    temp_v(:,:,1:58) = var_temp(:,:,1:58); % feb 28 
    temp_v(:,:,59) = var_temp(:,:,58); %make feb 29 (same with 28)
    temp_v(:,:,60:366)=var_temp(:,:,59:end);
    hur_temp=temp_v;
    
    ncid = netcdf.create(strcat('./total/NorESM-M_hur_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',366);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'relative humidity at the 2meter height');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','lon');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','lat');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    hurvarid = netcdf.defVar(ncid, 'hur', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,hurvarid,'standard_name','hur');
    netcdf.putAtt(ncid,hurvarid,'long_name','Relative humidity');
    netcdf.putAtt(ncid,hurvarid,'units','Percent');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, hurvarid, [0 0 0], [size11 size12 366], hur_temp);
    netcdf.putVar(ncid, timevarid, 0, length(t_temp), t_temp);
    netcdf.close(ncid);   
        
    
    else
        
    hur_temp=var_temp;
    ncid = netcdf.create(strcat('./total/NorESM-M_hur_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',365);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'relative humidity at the 2meter height');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','lon');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','lat');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    hurvarid = netcdf.defVar(ncid, 'hur', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,hurvarid,'standard_name','hur');
    netcdf.putAtt(ncid,hurvarid,'long_name','Relative humidity');
    netcdf.putAtt(ncid,hurvarid,'units','Percent');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
    netcdf.putVar(ncid, hurvarid, [0 0 0], [size11 size12 365], hur_temp);

    netcdf.close(ncid);
    end
 end
 
 
clc; clear all; close all;

minyear=1976;
maxyear=2005;
list = dir('*psl_interped_nor_*');

psl_sd_1 = load(list(2).name,'psl_sd');
psl_sd_2 = load(list(3).name,'psl_sd');

psl_sd_1= psl_sd_1.psl_sd;
psl_sd_2= psl_sd_2.psl_sd;
load('standard_grid.mat');

[size11 size12 size13]= size(psl_sd_1);
[size21 size22 size23]= size(psl_sd_2);

%%% concatenate array along specified dimension
psl_sd = cat(3, psl_sd_1, psl_sd_2); %combine matrix into one

t_sd=size(psl_sd,3)

timetime=minyear:1/365:maxyear+1-1/365;

tstep=0;
leap_yr = [1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004];
for j = minyear:maxyear
    tstep =tstep +1;
    fstep=(tstep-1)*365+1:1:tstep*365;
    ftime=timetime(fstep); % slice every year in time
    var_temp=psl_sd(:,:,fstep); % slice every year in psl

if find(leap_yr == j)
    t_temp(1:58) = ftime(1:58);
    t_temp(59) = ftime(58);
    t_temp(60:366) = ftime(59:end);
    
    temp_v(:,:,1:58) = var_temp(:,:,1:58); % feb 28 
    temp_v(:,:,59) = var_temp(:,:,58); %make feb 29 (same with 28)
    temp_v(:,:,60:366)=var_temp(:,:,59:end);
    psl_temp=temp_v;
    
    ncid = netcdf.create(strcat('./total/NorESM-M_psl_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',366);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'psl');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','latitude');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    pslvarid = netcdf.defVar(ncid, 'psl', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,pslvarid,'standard_name','psl');
    netcdf.putAtt(ncid,pslvarid,'long_name','Air pressure at the sea level');
    netcdf.putAtt(ncid,pslvarid,'units','Pa');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, pslvarid, [0 0 0], [size11 size12 366], psl_temp);

    netcdf.close(ncid);
else
    psl_temp=var_temp;
    ncid = netcdf.create(strcat('./total/NorESM-M_psl_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',365);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'psl');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','latitude');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    pslvarid = netcdf.defVar(ncid, 'psl', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,pslvarid,'standard_name','psl');
    netcdf.putAtt(ncid,pslvarid,'long_name','Air pressure at the sea level');
    netcdf.putAtt(ncid,pslvarid,'units','Pa');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, pslvarid, [0 0 0], [size11 size12 365], psl_temp);

    netcdf.close(ncid);
end
end

 close all; clear all; clc;

minyear=1976;
maxyear=2005;
list = dir('*rsds_interped_nor_*');

rsds_sd_1 = load(list(2).name,'rsds_sd');
rsds_sd_2 = load(list(3).name,'rsds_sd');

rsds_sd_1= rsds_sd_1.rsds_sd;
rsds_sd_2= rsds_sd_2.rsds_sd;
load('standard_grid.mat');

[size11 size12 size13]= size(rsds_sd_1);
[size21 size22 size23]= size(rsds_sd_2);

%%% concatenate array along specified dimension
rsds_sd = cat(3, rsds_sd_1, rsds_sd_2); %combine matrix into one

t_sd=size(rsds_sd,3)

timetime=minyear:1/365:maxyear+1-1/365;

tstep=0;
% for j = minyear:maxyear
leap_yr = [1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004];
for j = minyear:maxyear
    tstep =tstep +1;
    fstep=(tstep-1)*365+1:1:tstep*365;
    ftime=timetime(fstep); % slice every year in time
    var_temp=rsds_sd(:,:,fstep); % slice every year in rsds

if find(leap_yr == j)
    t_temp(1:58) = ftime(1:58);
    t_temp(59) = ftime(58);
    t_temp(60:366) = ftime(59:end);
    
    temp_v(:,:,1:58) = var_temp(:,:,1:58); % feb 28 
    temp_v(:,:,59) = var_temp(:,:,58); %make feb 29 (same with 28)
    temp_v(:,:,60:366)=var_temp(:,:,59:end);
    rsds_temp=temp_v;
    
    ncid = netcdf.create(strcat('./total/NorESM-M_rsds_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid,'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',366);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'solar radiation');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','latitude');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    rsdsvarid = netcdf.defVar(ncid, 'rsds', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,rsdsvarid,'standard_name','rsds');
    netcdf.putAtt(ncid,rsdsvarid,'long_name','Surface_downwelling_shortwave_flux_in_air');
    netcdf.putAtt(ncid,rsdsvarid,'units','W/m^2');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, rsdsvarid, [0 0 0], [size11 size12 366], rsds_temp);

    netcdf.close(ncid);
    
else
    rsds_temp=var_temp;
    ncid = netcdf.create(strcat('./total/NorESM-M_rsds_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid,'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',365);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'solar radiation');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','latitude');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    rsdsvarid = netcdf.defVar(ncid, 'rsds', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,rsdsvarid,'standard_name','rsds');
    netcdf.putAtt(ncid,rsdsvarid,'long_name','Surface_downwelling_shortwave_flux_in_air');
    netcdf.putAtt(ncid,rsdsvarid,'units','W/m^2');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, rsdsvarid, [0 0 0], [size11 size12 365], rsds_temp);

    netcdf.close(ncid);
end
end

clc; clear all; close all;

minyear=1976;
maxyear=2005;
list = dir('*tas_interped_nor_*');

tas_sd_1 = load(list(2).name,'tas_sd');
tas_sd_2 = load(list(3).name,'tas_sd');

tas_sd_1= tas_sd_1.tas_sd;
tas_sd_2= tas_sd_2.tas_sd;
load('standard_grid.mat');

[size11 size12 size13]= size(tas_sd_1);
[size21 size22 size23]= size(tas_sd_2);

%%% concatenate array along specified dimension
tas_sd = cat(3, tas_sd_1, tas_sd_2); %combine matrix into one

t_sd=size(tas_sd,3)

timetime=minyear:1/365:maxyear+1-1/365;

tstep=0;
leap_yr = [1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004];
for j = minyear:maxyear
    tstep =tstep +1;
    fstep=(tstep-1)*365+1:1:tstep*365;
    ftime=timetime(fstep); % slice every year in time
    var_temp=tas_sd(:,:,fstep); % slice every year in tas

if find(leap_yr == j)
    t_temp(1:58) = ftime(1:58);
    t_temp(59) = ftime(58);
    t_temp(60:366) = ftime(59:end);
    
    temp_v(:,:,1:58) = var_temp(:,:,1:58); % feb 28 
    temp_v(:,:,59) = var_temp(:,:,58); %make feb 29 (same with 28)
    temp_v(:,:,60:366)=var_temp(:,:,59:end);
    tas_temp=temp_v;
    
    ncid = netcdf.create(strcat('./total/NorESM-M_tas_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',366);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'Air temperature at the 2meter height');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','lon');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','lat');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    tasvarid = netcdf.defVar(ncid, 'tas', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,tasvarid,'standard_name','tas');
    netcdf.putAtt(ncid,tasvarid,'long_name','Air temperature at the surface');
    netcdf.putAtt(ncid,tasvarid,'units','W/m^2');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, tasvarid, [0 0 0], [size11 size12 366], tas_temp);

    netcdf.close(ncid);
    
else
    tas_temp=var_temp;
    ncid = netcdf.create(strcat('./total/NorESM-M_tas_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',365);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'Air temperature at the 2meter height');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','lon');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','lat');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    tasvarid = netcdf.defVar(ncid, 'tas', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,tasvarid,'standard_name','tas');
    netcdf.putAtt(ncid,tasvarid,'long_name','Air temperature at the surface');
    netcdf.putAtt(ncid,tasvarid,'units','W/m^2');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, tasvarid, [0 0 0], [size11 size12 365], tas_temp);

    netcdf.close(ncid);
end
end

clc; clear all; close all;

minyear=1976;
maxyear=2005;
list = dir('*uo_interped_nor_*');

uo_sd_1 = load(list(2).name,'uo_sd');
uo_sd_2 = load(list(3).name,'uo_sd');
uo_sd_3 = load(list(4).name,'uo_sd');
uo_sd_4 = load(list(5).name,'uo_sd');
uo_sd_1= uo_sd_1.uo_sd;
uo_sd_2= uo_sd_2.uo_sd;
uo_sd_3= uo_sd_3.uo_sd;
uo_sd_4= uo_sd_4.uo_sd;
load('standard_grid.mat');

[size11 size12 size13]= size(uo_sd_1);
[size21 size22 size23]= size(uo_sd_2);
[size31 size32 size33]= size(uo_sd_3);
[size41 size42 size43]= size(uo_sd_4);

%%% concatenate array along specified dimension
uo_sd = cat(3, uo_sd_1, uo_sd_2, uo_sd_3, uo_sd_4); %combine matrix into one

t_sd=size(uo_sd,3)

timetime=minyear:1/365:maxyear+1-1/365;

tstep=0;
leap_yr = [1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004];
for j = minyear:maxyear
    tstep =tstep +1;
    fstep=(tstep-1)*365+1:1:tstep*365;
    ftime=timetime(fstep); % slice every year in time
    var_temp=uo_sd(:,:,fstep); % slice every year in uo

if find(leap_yr == j)
    t_temp(1:58) = ftime(1:58);
    t_temp(59) = ftime(58);
    t_temp(60:366) = ftime(59:end);
    
    temp_v(:,:,1:58) = var_temp(:,:,1:58); % feb 28 
    temp_v(:,:,59) = var_temp(:,:,58); %make feb 29 (same with 28)
    temp_v(:,:,60:366)=var_temp(:,:,59:end);
    uo_temp=temp_v;
    
    ncid = netcdf.create(strcat('./total/NorESM-M_uo_wind_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',366);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'Eastward_wind at the 10meter above the surface');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','lon');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','lat');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    uovarid = netcdf.defVar(ncid, 'uo', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,uovarid,'standard_name','uo');
    netcdf.putAtt(ncid,uovarid,'long_name','Eastward_wind at the 10meter above the surface');
    netcdf.putAtt(ncid,uovarid,'units','m/s');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, uovarid, [0 0 0], [size11 size12 366], uo_temp);

    netcdf.close(ncid);
else
    uo_temp=var_temp;
    ncid = netcdf.create(strcat('./total/NorESM-M_uo_wind_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',365);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'Eastward_wind at the 10meter above the surface');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','lon');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','lat');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    uovarid = netcdf.defVar(ncid, 'uo', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,uovarid,'standard_name','uo');
    netcdf.putAtt(ncid,uovarid,'long_name','Eastward_wind at the 10meter above the surface');
    netcdf.putAtt(ncid,uovarid,'units','m/s');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, uovarid, [0 0 0], [size11 size12 365], uo_temp);

    netcdf.close(ncid);
end
end

close all; clear all; clc;

minyear=1976;
maxyear=2005;
list = dir('*vo_interped_nor_*');

vo_sd_1 = load(list(2).name,'vo_sd');
vo_sd_2 = load(list(3).name,'vo_sd');
vo_sd_3 = load(list(4).name,'vo_sd');
vo_sd_4 = load(list(5).name,'vo_sd');
vo_sd_1= vo_sd_1.vo_sd;
vo_sd_2= vo_sd_2.vo_sd;
vo_sd_3= vo_sd_3.vo_sd;
vo_sd_4= vo_sd_4.vo_sd;
load('standard_grid.mat');

[size11 size12 size13]= size(vo_sd_1);
[size21 size22 size23]= size(vo_sd_2);
[size31 size32 size33]= size(vo_sd_3);
[size41 size42 size43]= size(vo_sd_4);

%%% concatenate array along specified dimension
vo_sd = cat(3, vo_sd_1, vo_sd_2, vo_sd_3, vo_sd_4); %combine matrix into one

t_sd=size(vo_sd,3)

timetime=minyear:1/365:maxyear+1-1/365;

tstep =0;
leap_yr = [1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004];
for j = minyear:maxyear
    tstep =tstep +1;
    fstep=(tstep-1)*365+1:1:tstep*365;
    ftime=timetime(fstep); % slice every year in time
    var_temp=vo_sd(:,:,fstep); % slice every year in vo

if find(leap_yr == j)
    t_temp(1:58) = ftime(1:58);
    t_temp(59) = ftime(58);
    t_temp(60:366) = ftime(59:end);
    
    temp_v(:,:,1:58) = var_temp(:,:,1:58); % feb 28 
    temp_v(:,:,59) = var_temp(:,:,58); %make feb 29 (same with 28)
    temp_v(:,:,60:366)=var_temp(:,:,59:end);
    vo_temp=temp_v;
    
        ncid = netcdf.create(strcat('./total/NorESM-M_vo_wind_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',366);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'Northward_wind at the 10meter above the surface');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','lon');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','lat');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    vovarid = netcdf.defVar(ncid, 'vo', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,vovarid,'standard_name','vo');
    netcdf.putAtt(ncid,vovarid,'long_name','Northward_wind at the 10meter above the surface');
    netcdf.putAtt(ncid,vovarid,'units','m/s');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, vovarid, [0 0 0], [size11 size12 366], vo_temp);

    netcdf.close(ncid);
    
else
    vo_temp=var_temp;
    ncid = netcdf.create(strcat('./total/NorESM-M_vo_wind_',num2str(j),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    time_dimid = netcdf.defDim(ncid, 'time',365);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'Northward_wind at the 10meter above the surface');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','lon');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','lat');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    vovarid = netcdf.defVar(ncid, 'vo', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,vovarid,'standard_name','vo');
    netcdf.putAtt(ncid,vovarid,'long_name','Northward_wind at the 10meter above the surface');
    netcdf.putAtt(ncid,vovarid,'units','m/s');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
    netcdf.putVar(ncid, vovarid, [0 0 0], [size11 size12 365], vo_temp);

    netcdf.close(ncid);
end
end


cd '/data1/cshwa/tmp/history/ocean'
clc; clear all; close all;

minyear=1974;
maxyear=2005;
list = dir('*uo_interped_nor_*')

uo_sd_1 = load(list(2).name,'uo_sd');
uo_sd_2 = load(list(3).name,'uo_sd');
uo_sd_3 = load(list(4).name,'uo_sd');
uo_sd_4 = load(list(5).name,'uo_sd');
uo_sd_5 = load(list(6).name,'uo_sd');
uo_sd_6 = load(list(7).name,'uo_sd');
uo_sd_7 = load(list(8).name,'uo_sd');
uo_sd_8 = load(list(9).name,'uo_sd');
uo_sd_1= uo_sd_1.uo_sd;
uo_sd_2= uo_sd_2.uo_sd;
uo_sd_3= uo_sd_3.uo_sd;
uo_sd_4= uo_sd_4.uo_sd;
uo_sd_5= uo_sd_5.uo_sd;
uo_sd_6= uo_sd_6.uo_sd;
uo_sd_7= uo_sd_7.uo_sd;
uo_sd_8= uo_sd_8.uo_sd;
load('standard_grid_ocean.mat');

[size11 size12 size13 size14]= size(uo_sd_1);

%%% concatenate array along specified dimension
uo_sd = cat(4, uo_sd_1, uo_sd_2, uo_sd_3, uo_sd_4,uo_sd_5,uo_sd_6,uo_sd_7,uo_sd_8); %combine matrix into one

t_sd=size(uo_sd,4)
size(uo_sd)

timetime=minyear:1/365:maxyear+1-1/365;

tstep=0;
% for j = minyear:maxyear
for j = minyear:maxyear
tstep =tstep +1;
fstep=(tstep-1)*12+1:1:tstep*12;
ftime=timetime(fstep); % slice every year in time
uo_temp=uo_sd(:,:,:,fstep); % slice every year in uo

ncid = netcdf.create(strcat('../total/NorESM-M_uo_',num2str(j),'.nc'),'CLOBBER');
lon_dimid = netcdf.defDim(ncid,'lon',size11);
lat_dimid = netcdf.defDim(ncid, 'lat', size12);
dep_dimid = netcdf.defDim(ncid, 'depth', size13);
time_dimid = netcdf.defDim(ncid, 'time',12);  

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' NorESM1-M interpolated data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'title', 'sea_water_x_velocity');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'uource', 'CMIP5 data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by S.H Chae');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
netcdf.putAtt(ncid,lonvarid,'long_name','lon');
netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
netcdf.putAtt(ncid,latvarid,'standard_name','lat');
netcdf.putAtt(ncid,latvarid,'long_name','lat');
netcdf.putAtt(ncid,latvarid,'units','degrees_north');

depvarid = netcdf.defVar(ncid, 'depth', 'NC_float', dep_dimid);
netcdf.putAtt(ncid,depvarid,'standard_name','ocean depth');
netcdf.putAtt(ncid,depvarid,'long_name','depth_below_ocean_surface');
netcdf.putAtt(ncid,depvarid,'units','meter');

timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
netcdf.putAtt(ncid,timevarid,'standard_name','time');
netcdf.putAtt(ncid,timevarid,'long_name','time');
netcdf.putAtt(ncid,timevarid,'units','months from 1974-jan');

uovarid = netcdf.defVar(ncid, 'uo', 'NC_float', [lon_dimid lat_dimid dep_dimid time_dimid]);
netcdf.putAtt(ncid,uovarid,'standard_name','uo');
netcdf.putAtt(ncid,uovarid,'long_name','sea_water_x_velocity');
netcdf.putAtt(ncid,uovarid,'units','m/s');

netcdf.endDef(ncid);

netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1,1))), squeeze(lon_ec_f(:,1,1)));
netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:,1))), squeeze(lat_ec_f(1,:,1)));
netcdf.putVar(ncid, depvarid, 0, length(squeeze(stan_dep(1,1,:))), squeeze(stan_dep(1,1,:)));
netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
netcdf.putVar(ncid, uovarid, [0 0 0 0], [size11 size12 size13 12], uo_temp);

netcdf.close(ncid);
end

clc; clear all; close all;

minyear=1974;
maxyear=2005;
list = dir('*vo_interped_nor_*')

vo_sd_1 = load(list(2).name,'vo_sd');
vo_sd_2 = load(list(3).name,'vo_sd');
vo_sd_3 = load(list(4).name,'vo_sd');
vo_sd_4 = load(list(5).name,'vo_sd');
vo_sd_5 = load(list(6).name,'vo_sd');
vo_sd_6 = load(list(7).name,'vo_sd');
vo_sd_7 = load(list(8).name,'vo_sd');
vo_sd_8 = load(list(9).name,'vo_sd');
vo_sd_1= vo_sd_1.vo_sd;
vo_sd_2= vo_sd_2.vo_sd;
vo_sd_3= vo_sd_3.vo_sd;
vo_sd_4= vo_sd_4.vo_sd;
vo_sd_5= vo_sd_5.vo_sd;
vo_sd_6= vo_sd_6.vo_sd;
vo_sd_7= vo_sd_7.vo_sd;
vo_sd_8= vo_sd_8.vo_sd;
load('standard_grid_ocean.mat');

[size11 size12 size13 size14]= size(vo_sd_1);

%%% concatenate array along specified dimension
vo_sd = cat(4, vo_sd_1, vo_sd_2, vo_sd_3, vo_sd_4,vo_sd_5,vo_sd_6,vo_sd_7,vo_sd_8); %combine matrix into one

t_sd=size(vo_sd,4)
size(vo_sd)

timetime=minyear:1/365:maxyear+1-1/365;

tstep=0;
% for j = minyear:maxyear
for j = minyear:maxyear
tstep =tstep +1;
fstep=(tstep-1)*12+1:1:tstep*12;
ftime=timetime(fstep); % slice every year in time
vo_temp=vo_sd(:,:,:,fstep); % slice every year in vo

ncid = netcdf.create(strcat('../total/NorESM-M_vo_',num2str(j),'.nc'),'CLOBBER');
lon_dimid = netcdf.defDim(ncid,'lon',size11);
lat_dimid = netcdf.defDim(ncid, 'lat', size12);
dep_dimid = netcdf.defDim(ncid, 'depth', size13);
time_dimid = netcdf.defDim(ncid, 'time',12);  

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' NorESM1-M interpolated data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'title', 'sea_water_y_velocity');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'vource', 'CMIP5 data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by S.H Chae');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
netcdf.putAtt(ncid,lonvarid,'long_name','lon');
netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
netcdf.putAtt(ncid,latvarid,'standard_name','lat');
netcdf.putAtt(ncid,latvarid,'long_name','lat');
netcdf.putAtt(ncid,latvarid,'units','degrees_north');

depvarid = netcdf.defVar(ncid, 'depth', 'NC_float', dep_dimid);
netcdf.putAtt(ncid,depvarid,'standard_name','ocean depth');
netcdf.putAtt(ncid,depvarid,'long_name','depth_below_ocean_surface');
netcdf.putAtt(ncid,depvarid,'units','meter');

timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
netcdf.putAtt(ncid,timevarid,'standard_name','time');
netcdf.putAtt(ncid,timevarid,'long_name','time');
netcdf.putAtt(ncid,timevarid,'units','months from 1974-jan');

vovarid = netcdf.defVar(ncid, 'vo', 'NC_float', [lon_dimid lat_dimid dep_dimid time_dimid]);
netcdf.putAtt(ncid,vovarid,'standard_name','vo');
netcdf.putAtt(ncid,vovarid,'long_name','sea_water_y_velocity');
netcdf.putAtt(ncid,vovarid,'units','m/s');

netcdf.endDef(ncid);

netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1,1))), squeeze(lon_ec_f(:,1,1)));
netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:,1))), squeeze(lat_ec_f(1,:,1)));
netcdf.putVar(ncid, depvarid, 0, length(squeeze(stan_dep(1,1,:))), squeeze(stan_dep(1,1,:)));
netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
netcdf.putVar(ncid, vovarid, [0 0 0 0], [size11 size12 size13 12], vo_temp);

netcdf.close(ncid);
end

clc; clear all; close all;

minyear=1974;
maxyear=2005;
list = dir('*so_interped_nor_*');

so_sd_1 = load(list(2).name,'so_sd');
so_sd_2 = load(list(3).name,'so_sd');
so_sd_3 = load(list(4).name,'so_sd');
so_sd_4 = load(list(5).name,'so_sd');
so_sd_5 = load(list(6).name,'so_sd');
so_sd_6 = load(list(7).name,'so_sd');
so_sd_7 = load(list(8).name,'so_sd');
so_sd_8 = load(list(9).name,'so_sd');
so_sd_1= so_sd_1.so_sd;
so_sd_2= so_sd_2.so_sd;
so_sd_3= so_sd_3.so_sd;
so_sd_4= so_sd_4.so_sd;
so_sd_5= so_sd_5.so_sd;
so_sd_6= so_sd_6.so_sd;
so_sd_7= so_sd_7.so_sd;
so_sd_8= so_sd_8.so_sd;
load('standard_grid_ocean.mat');

[size11 size12 size13 size14]= size(so_sd_1);

%%% concatenate array along specified dimension
so_sd = cat(4, so_sd_1, so_sd_2, so_sd_3, so_sd_4,so_sd_5,so_sd_6,so_sd_7,so_sd_8); %combine matrix into one

t_sd=size(so_sd,4)
size(so_sd)

timetime=minyear:1/365:maxyear+1-1/365;

tstep=0;
% for j = minyear:maxyear
for j = minyear:maxyear
tstep =tstep +1;
fstep=(tstep-1)*12+1:1:tstep*12;
ftime=timetime(fstep); % slice every year in time
so_temp=so_sd(:,:,:,fstep); % slice every year in so

ncid = netcdf.create(strcat('../total/NorESM-M_so_',num2str(j),'.nc'),'CLOBBER');
lon_dimid = netcdf.defDim(ncid,'lon',size11);
lat_dimid = netcdf.defDim(ncid, 'lat', size12);
dep_dimid = netcdf.defDim(ncid, 'depth', size13);
time_dimid = netcdf.defDim(ncid, 'time',12);  

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' NorESM1-M interpolated data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'title', 'sea_water_salinity');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'source', 'CMIP5 data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by S.H Chae');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
netcdf.putAtt(ncid,lonvarid,'long_name','lon');
netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
netcdf.putAtt(ncid,latvarid,'standard_name','lat');
netcdf.putAtt(ncid,latvarid,'long_name','lat');
netcdf.putAtt(ncid,latvarid,'units','degrees_north');

depvarid = netcdf.defVar(ncid, 'depth', 'NC_float', dep_dimid);
netcdf.putAtt(ncid,depvarid,'standard_name','ocean depth');
netcdf.putAtt(ncid,depvarid,'long_name','depth_below_ocean_surface');
netcdf.putAtt(ncid,depvarid,'units','meter');

timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
netcdf.putAtt(ncid,timevarid,'standard_name','time');
netcdf.putAtt(ncid,timevarid,'long_name','time');
netcdf.putAtt(ncid,timevarid,'units','months from 1974-jan');

sovarid = netcdf.defVar(ncid, 'salt', 'NC_float', [lon_dimid lat_dimid dep_dimid time_dimid]);
netcdf.putAtt(ncid,sovarid,'standard_name','salt');
netcdf.putAtt(ncid,sovarid,'long_name','sea_water_salinity');
netcdf.putAtt(ncid,sovarid,'units','PSU');

netcdf.endDef(ncid);

netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1,1))), squeeze(lon_ec_f(:,1,1)));
netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:,1))), squeeze(lat_ec_f(1,:,1)));
netcdf.putVar(ncid, depvarid, 0, length(squeeze(stan_dep(1,1,:))), squeeze(stan_dep(1,1,:)));
netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
netcdf.putVar(ncid, sovarid, [0 0 0 0], [size11 size12 size13 12], so_temp);

netcdf.close(ncid);
end

clc; clear all; close all;

minyear=1974;
maxyear=2005;
list = dir('*thetao_interped_nor_*');

thetao_sd_1 = load(list(2).name,'thetao_sd');
thetao_sd_2 = load(list(3).name,'thetao_sd');
thetao_sd_3 = load(list(4).name,'thetao_sd');
thetao_sd_4 = load(list(5).name,'thetao_sd');
thetao_sd_5 = load(list(6).name,'thetao_sd');
thetao_sd_6 = load(list(7).name,'thetao_sd');
thetao_sd_7 = load(list(8).name,'thetao_sd');
thetao_sd_8 = load(list(9).name,'thetao_sd');
thetao_sd_1= thetao_sd_1.thetao_sd;
thetao_sd_2= thetao_sd_2.thetao_sd;
thetao_sd_3= thetao_sd_3.thetao_sd;
thetao_sd_4= thetao_sd_4.thetao_sd;
thetao_sd_5= thetao_sd_5.thetao_sd;
thetao_sd_6= thetao_sd_6.thetao_sd;
thetao_sd_7= thetao_sd_7.thetao_sd;
thetao_sd_8= thetao_sd_8.thetao_sd;
load('standard_grid_ocean.mat');

[size11 size12 size13 size14]= size(thetao_sd_1);

%%% concatenate array along specified dimension
thetao_sd = cat(4, thetao_sd_1, thetao_sd_2, thetao_sd_3, thetao_sd_4,thetao_sd_5,thetao_sd_6,thetao_sd_7,thetao_sd_8); %combine matrix into one

t_sd=size(thetao_sd,4)
size(thetao_sd)

timetime=minyear:1/365:maxyear+1-1/365;

tstep=0;
% for j = minyear:maxyear
for j = minyear:maxyear
tstep =tstep +1;
fstep=(tstep-1)*12+1:1:tstep*12;
ftime=timetime(fstep); % slice every year in time
thetao_temp=thetao_sd(:,:,:,fstep); % slice every year in thetao

ncid = netcdf.create(strcat('../total/NorESM-M_temp_',num2str(j),'.nc'),'CLOBBER');
lon_dimid = netcdf.defDim(ncid,'lon',size11);
lat_dimid = netcdf.defDim(ncid, 'latitude', size12);
dep_dimid = netcdf.defDim(ncid, 'depth', size13);
time_dimid = netcdf.defDim(ncid, 'time',12);  

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' NorESM1-M interpolated data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'title', 'sea_water_potential_temperature');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'thetaource', 'CMIP5 data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by S.H Chae');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
netcdf.putAtt(ncid,lonvarid,'long_name','lon');
netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

latvarid = netcdf.defVar(ncid, 'latitude', 'NC_float', lat_dimid);
netcdf.putAtt(ncid,latvarid,'standard_name','latitude');
netcdf.putAtt(ncid,latvarid,'long_name','latitude');
netcdf.putAtt(ncid,latvarid,'units','degrees_north');

depvarid = netcdf.defVar(ncid, 'depth', 'NC_float', dep_dimid);
netcdf.putAtt(ncid,depvarid,'standard_name','ocean depth');
netcdf.putAtt(ncid,depvarid,'long_name','depth_below_ocean_surface');
netcdf.putAtt(ncid,depvarid,'units','meter');

timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
netcdf.putAtt(ncid,timevarid,'standard_name','time');
netcdf.putAtt(ncid,timevarid,'long_name','time');
netcdf.putAtt(ncid,timevarid,'units','months from 1974-jan');

thetaovarid = netcdf.defVar(ncid, 'temp', 'NC_float', [lon_dimid lat_dimid dep_dimid time_dimid]);
netcdf.putAtt(ncid,thetaovarid,'standard_name','temp');
netcdf.putAtt(ncid,thetaovarid,'long_name','sea_water_potential_temperature');
netcdf.putAtt(ncid,thetaovarid,'units','K');

netcdf.endDef(ncid);

netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1,1))), squeeze(lon_ec_f(:,1,1)));
netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:,1))), squeeze(lat_ec_f(1,:,1)));
netcdf.putVar(ncid, depvarid, 0, length(squeeze(stan_dep(1,1,:))), squeeze(stan_dep(1,1,:)));
netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
netcdf.putVar(ncid, thetaovarid, [0 0 0 0], [size11 size12 size13 12], thetao_temp);

netcdf.close(ncid);
end


close all; clear all; clc;

minyear=1850;
maxyear=2005;
list = dir('zos_interped_nor*');

zos_sd_1 = load(list(1).name,'zos_sd');
zos_sd_1= zos_sd_1.zos_sd;

load('standard_grid.mat');

[size11 size12 size13] = size(zos_sd_1);

%%% concatenate array along specified dimension
zos_sd =zos_sd_1;

t_sd = size(zos_sd,4)
size(zos_sd)

timetime=minyear:1/12:maxyear+1-1/12;

tstep=0;
for j = minyear:maxyear
tstep = tstep + 1;
fstep=(tstep-1)*12+1:1:tstep*12;
ftime=timetime(fstep); % slice every year in time
zos_temp=zos_sd(:,:,fstep); % slice every year in zos

ncid = netcdf.create(strcat('../total/NorESM-M_zos_',num2str(j),'.nc'),'CLOBBER');
lon_dimid = netcdf.defDim(ncid,'lon',size11);
lat_dimid = netcdf.defDim(ncid, 'lat', size12);
time_dimid = netcdf.defDim(ncid, 'time',12);  

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' NorESM1-M interpolated data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'title', 'sea_surface_height_above geoid');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'zosurce', 'CMIP5 data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by S.H Chae');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
netcdf.putAtt(ncid,lonvarid,'long_name','lon');
netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
netcdf.putAtt(ncid,latvarid,'standard_name','lat');
netcdf.putAtt(ncid,latvarid,'long_name','lat');
netcdf.putAtt(ncid,latvarid,'units','degrees_north');

timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
netcdf.putAtt(ncid,timevarid,'standard_name','time');
netcdf.putAtt(ncid,timevarid,'long_name','time');
netcdf.putAtt(ncid,timevarid,'units','months from 1850-jan');

zosvarid = netcdf.defVar(ncid, 'zos', 'NC_float', [lon_dimid lat_dimid time_dimid]);
netcdf.putAtt(ncid,zosvarid,'standard_name','zos');
netcdf.putAtt(ncid,zosvarid,'long_name','sea_surface_height_above geoid');
netcdf.putAtt(ncid,zosvarid,'units','meter');

netcdf.endDef(ncid);

netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1))), squeeze(lon_ec_f(:,1)));
netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:))), squeeze(lat_ec_f(1,:)));
netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
netcdf.putVar(ncid, zosvarid, [0 0 0], [size11 size12 12], zos_temp);
netcdf.close(ncid);
end