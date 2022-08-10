% % % This script is based on MATLAB R2017a
% % 21-Jun-2018

clear all;close all; clc;
warning off;

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
%     dropboxpath='C:\Users\KYY\Dropbox';
%     addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
%     addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_ktotalday']));
    addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path

elseif (strcmp(system_name,'GLNXA64'))
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_ktotalday']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Run']));
    addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
end
% start


grdfile='/home/yjtak/roms/auto_NWP/input/auto01_grid.nc';

rawdata_dir='/data2/temp/ECMWF_ERA5/';
varname={
%     '2m_temperature', 'mean_sea_level_pressure', ...
%     '2m_dewpoint_temperature', ...
%     '10m_u_component_of_wind', '10m_v_component_of_wind',...
%  };
     'total_precipitation','surface_solar_radiation_downwards'
  };
% varname={'surface_solar_radiation_downwards'};
year=[1991:2000];
tinterval = 24; % unit : hours

testname = 'autoNWP_ERA5';
output_dir=['/data1/yjtak/auto_fennel/input/ERA5/'];
    
% ncread('/data1/temp/ECMWF_interim/ssrd/ECMWF_Interim_ssrd_2016.nc',ssrd);
for i=1:length(year)
    for j=1:length(varname)
        status = ROMS_surface_forcing_ERA5(grdfile, year(i), varname{j}, tinterval,output_dir,testname);
    end
    status = ROMS_surface_forcing_Qair(grdfile, year(i), tinterval,output_dir,testname);
end