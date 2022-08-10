clear all;close all;clc
addpath(genpath('/home/yjtak/roms/matlab'));
addpath(genpath('/home/yjtak/Dropbox/source/matlab/Common/m_map'));

gn = grd('auto8');
Vtransform=2;
Vstretching=4;
% (3) roms data file from a larger domain
head='test06_monthly_';
foot='.nc';

% sourcestr = [ 'From 1/4 deg. : to 10km resolution '];
% check file number   ex) avg_1606.nc ~ avg_1666.nc
% start_file = 0050;  %% use//  file_type = 'sequence_number';
start_file = 1;   %% use//  file_type = 'year_mon_number';
end_file = 12;
bndy_time=[15:30:365];
%file_type = 'sequence_number'; %% sequence_number // 0001.nc
file_type = 'year_mon_number'; %% year_mon_number // YYYY_MM.nc


% (5) output file, boundary file for the nested grid
% out_file = 'roms_yw_init1996_layer20.nc';
% mkdir(out_head);
out_head = ['./'];

out_file = [out_head,'roms_bndy_auto_NWP_1_10_test6_clim06_10',foot];
in_file_h = [out_head,'roms_bndy_auto_NWP_1_10_test6_'];

disp([' create an empty boundary file = ' out_file])
noclobber = 0;
base_date = [2006 1 1 0 0 0];
time_variable = 'ocean_time';
cycle_length = 360;
time_variable_units = basedate2str(base_date);
roms.grd = gn;
grd_file = gn.grd_file;
details = [ 'Quick boundary file from NWP 1/10 deg ' ...
    'by script ' which(mfilename) ];
donuts = 0;
nc_create_roms_bndy_nest_new  % create an nc file, write variables and close it.
disp([' creating an empty boundary file done ..... '])
disp('  ')


st_year=2006;
end_year=2010;
ic=0;
for year=st_year:end_year
    ic=ic+1;
    in_file=[in_file_h,num2str(year),foot];
    nc=netcdf(in_file,'r');
    for varlist = { 'temp','salt','u','v'}
        varname = char(varlist);
        eval([varname '_w(ic,:,:,:) = nc{''',varname,'_west''}(:,:,:);']);
        eval([varname '_e(ic,:,:,:) = nc{''',varname,'_east''}(:,:,:);']);
        eval([varname '_n(ic,:,:,:) = nc{''',varname,'_north''}(:,:,:);']);
        eval([varname '_s(ic,:,:,:) = nc{''',varname,'_south''}(:,:,:);']);
    end
    
    for varlist = { 'zeta','ubar','vbar'}
        varname = char(varlist);
        eval([varname '_w(ic,:,:) = nc{''',varname,'_west''}(:,:);']);
        eval([varname '_e(ic,:,:) = nc{''',varname,'_east''}(:,:);']);
        eval([varname '_n(ic,:,:) = nc{''',varname,'_north''}(:,:);']);
        eval([varname '_s(ic,:,:) = nc{''',varname,'_south''}(:,:);']);
    end
    
    
end

for varlist ={'temp','salt','u','v'}
    varname=char(varlist);
    for dir ={'e','w','s','n'}
        dirname=char(dir);
        eval([varname,'_',dirname,'_clim = squeeze(mean(',varname,'_',dirname,'));']);
        eval([varname,'_',dirname,'_permute = permute(',varname,'_',dirname,'_clim,[3,2,1]);']);
        eval(['netcdf.putVar(ncid,',varname,'_',dirname,'_ID,',varname,'_',dirname,'_permute);']);
    end

end

for varlist ={'ubar','vbar','zeta'}
    varname=char(varlist);
    for dir ={'e','w','s','n'}
        dirname=char(dir);
        eval([varname,'_',dirname,'_clim = squeeze(mean(',varname,'_',dirname,'));']);
        eval([varname,'_',dirname,'_permute = permute(',varname,'_',dirname,'_clim,[2,1]);']);
        eval(['netcdf.putVar(ncid,',varname,'_',dirname,'_ID,',varname,'_',dirname,'_permute);']);
    end

end

