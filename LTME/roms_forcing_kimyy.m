% % % This script is based on MATLAB R2017a

clear all;close all; clc;
warning off;
linux=1; windows=0;
if (windows ==1)
    % % for windows
%     dropboxpath='C:\Users\KYY\Dropbox';
%     addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
%     addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_ktotalday']));
elseif (linux==1)
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_ktotalday']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Run']));
end
% start

grdfile='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/input/test37/2001/roms_grid_combine2_test37.nc';

rawdata_dir='/data1/temp/ECMWF_interim/';
varname={'airT', 'msl', 'dewt', 'ssrd', 'u10', 'v10'};
year=[1980];
tinterval = 1;
output_dir='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/forcing_matlab/';

% ncread('/data1/temp/ECMWF_interim/ssrd/ECMWF_Interim_ssrd_2016.nc',ssrd);
for i=1:length(year)
    for j=1:length(varname)
%         status = ROMS_surface_forcing_ECMWF(grdfile, year(i), varname{j}, tinterval);
    end
    status = ROMS_surface_forcing_Qair(grdfile, year(i), tinterval);
end