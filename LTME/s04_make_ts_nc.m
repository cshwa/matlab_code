clc;clear all;close all
% 
% Xi=[129.410:0.0001:129.435];lengx=length(Xi);
% Yi=[36.020:0.0001:36.035];lengy=length(Yi);
% Zi=[0:0.5:16];lengz=length(Zi);

load data\sumjin_constant_ts.mat
% load E:\[10]model\ROMS\03_initial_TS\data\sumjin_small_ts.mat.


% f_NO3=f_salt;f_NO3(1:4,:)=16.6;f_NO3(5,:,:)=(16.6+10.4)/2;f_NO3(6:15,:)=10.4;
% f_NH4=f_salt;f_NH4(1:4,:)=189;f_NH4(5,:,:)=(189+127)/2;f_NH4(6:15,:)=127;

%% temp
nw=netcdf('temp_GY.nc','clobber');

nw('X') = length(Xi);
nw('Y') = length(Yi);
% nw('T') = 12;
nw('Z') = length(z_dep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=Xi;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=Yi;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=z_dep';

nw{'temperature'} = ncfloat('Z','Y','X');
nw{'temperature'}.long_name = ncchar('temperature');
nw{'temperature'}.units = ncchar('deg.C');
nw{'temperature'}.missing_value=ncfloat(-99.9999008178711);
nw{'temperature'}(:)=f_temp;

close(nw);

%% salt
nw=netcdf('salt_GY.nc','clobber');

nw('X') = length(Xi);
nw('Y') = length(Yi);
% nw('T') = 12;
nw('Z') = length(z_dep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=Xi;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=Yi;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=z_dep;

nw{'salinity'} = ncfloat('Z','Y','X');
nw{'salinity'}.long_name = ncchar('salinity');
nw{'salinity'}.units = ncchar('salinity');
nw{'salinity'}.missing_value=ncfloat(-99.9999008178711);
nw{'salinity'}(:)=f_salt;

close(nw);



%% NO3

%{
nw=netcdf('NO3_pohang.nc','clobber');

nw('X') = length(Xi);
nw('Y') = length(Yi);
% nw('T') = 12;
nw('Z') = length(z_dep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=Xi;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=Yi;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=z_dep;

nw{'NO3'} = ncfloat('Z','Y','X');
nw{'NO3'}.long_name = ncchar('NO3');
nw{'NO3'}.units = ncchar('micro mole N/L');
nw{'NO3'}.missing_value=ncfloat(-99.9999008178711);
nw{'NO3'}(:)=f_NO3;

close(nw);

%% NH4
nw=netcdf('NH4_pohang.nc','clobber');

nw('X') = length(Xi);
nw('Y') = length(Yi);
% nw('T') = 12;
nw('Z') = length(z_dep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=Xi;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=Yi;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=z_dep;

nw{'NH4'} = ncfloat('Z','Y','X');
nw{'NH4'}.long_name = ncchar('NH4');
nw{'NH4'}.units = ncchar('micro mole N/L');
nw{'NH4'}.missing_value=ncfloat(-99.9999008178711);
nw{'NH4'}(:)=f_NH4;

close(nw);


%}

