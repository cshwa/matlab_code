clc;clear all;close all
% 
% Xi=[129.410:0.0001:129.435];lengx=length(Xi);
% Yi=[36.020:0.0001:36.035];lengy=length(Yi);
% Zi=[0:0.5:16];lengz=length(Zi);

%%  txt     txt
%%  lat     lon
%%  1   2    3       4           5      6          7          8      9
%% Temp	DO	Cond	Salinity	pH	Turbidity	Chlorophyll	 COD	TOC    
%% 10      11      12     13       14    15    16   17      18
%%NO2-N   NO3-N   NH4-N	  T-N     PO4-P	 TP    Ca	Mg      Fe

[raw1,txt1]=xlsread('H:\E\포항모델\포항측자료\서울대.xlsx','Sheet1','');

% unit conversion
temp_obs = raw1(1:24,1);
temp_obs = reshape(temp_obs,3,8);
DO_obs = raw1(1:24,2)*0.7*44.661;    % mg/L -> millimole_oxygen meter-3;
DO_obs = reshape(DO_obs,3,8);
salt_obs = raw1(1:24,4);
salt_obs = reshape(salt_obs,3,8);
chla_obs = raw1(1:24,7)/1000;
chla_obs = reshape(chla_obs,3,8);
no3_obs = raw1(1:24,11)*1000/14;  % mg/L -> millimole_N meter-3;
no3_obs = reshape(no3_obs,3,8);
NH4_obs = raw1(1:24,12)*1000/14;  % mg/L -> millimole_N meter-3;
NH4_obs = reshape(NH4_obs,3,8);
% % tPO4
% raw1(6,14) = 0;
% tPO4_obs_pre = raw1(1:24,14); tPO4_obs_pre(tPO4_obs_pre==0) = min(tPO4_obs_pre(tPO4_obs_pre~=0));
% tPO4_obs = tPO4_obs_pre*1000/30.9;  % mg/L -> millimole_N meter-3 용존총인자료 넣음.;
tPO4_obs = ones(3,8).*(0.01*(1000/30.9));

% deg, minute, sec unit to decimal degree

[raw2,txt2]=xlsread('H:\E\포항모델\포항측자료\서울대_변환.xlsx','Sheet4','');
lon_obs = raw2(1:8,2);
lat_obs = raw2(1:8,1);

% lon_obs= reshape(lon_obs,2,4);
% lat_obs= reshape(lat_obs,2,4);
% temp1 = reshape(squeeze(temp_obs(1,:)),2,4);

[lat_me lon_me]= meshgrid(lat_obs,lon_obs);
% lon=ncread('grid_pohang_csh_fine.nc','lon_rho');
% lat=ncread('grid_pohang_csh_fine.nc','lat_rho');
lon = (129.4200 : 0.0001 : 129.4350);
lat = (36.0200:0.0001:36.0330);
[lat lon]= meshgrid(lat,lon);

%% temp

sur_temp= griddata(lon_obs, lat_obs, squeeze(temp_obs(1,:)),lon_me,lat_me);
sur_temp_in = griddata(lon_obs, lat_obs, sur_temp,lon,lat);
% Now find the nan's
nanLocations = isnan(sur_temp_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(sur_temp_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(sur_temp_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(sur_temp_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = sur_temp_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	sur_temp_in(x,y,z) = goodValue;
end
% u should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(sur_temp_in);
numberOfNans = sum(nanLocations(:));

mid_temp = griddata(lon_obs, lat_obs, squeeze(temp_obs(2,:)),lon_me,lat_me);
mid_temp_in = griddata(lon_obs, lat_obs, mid_temp,lon,lat);
% Now find the nan's
nanLocations = isnan(mid_temp_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(mid_temp_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(mid_temp_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(mid_temp_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = mid_temp_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	mid_temp_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(mid_temp_in);
numberOfNans = sum(nanLocations(:));

bot_temp = griddata(lon_obs, lat_obs, squeeze(temp_obs(3,:)),lon_me,lat_me);
bot_temp_in = griddata(lon_obs, lat_obs, bot_temp,lon,lat);
% Now find the nan's
nanLocations = isnan(bot_temp_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(bot_temp_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(bot_temp_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(bot_temp_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = bot_temp_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	bot_temp_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(bot_temp_in);
numberOfNans = sum(nanLocations(:));


%% salt

sur_salt= griddata(lon_obs, lat_obs, squeeze(salt_obs(1,:)),lon_me,lat_me);
sur_salt_in = griddata(lon_obs, lat_obs, sur_salt,lon,lat);
% Now find the nan's
nanLocations = isnan(sur_salt_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(sur_salt_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(sur_salt_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(sur_salt_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = sur_salt_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	sur_salt_in(x,y,z) = goodValue;
end
% u should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(sur_salt_in);
numberOfNans = sum(nanLocations(:));

mid_salt = griddata(lon_obs, lat_obs, squeeze(salt_obs(2,:)),lon_me,lat_me);
mid_salt_in = griddata(lon_obs, lat_obs, mid_salt,lon,lat);
% Now find the nan's
nanLocations = isnan(mid_salt_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(mid_salt_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(mid_salt_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(mid_salt_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = mid_salt_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	mid_salt_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(mid_salt_in);
numberOfNans = sum(nanLocations(:));

bot_salt = griddata(lon_obs, lat_obs, squeeze(salt_obs(3,:)),lon_me,lat_me);
bot_salt_in = griddata(lon_obs, lat_obs, bot_salt,lon,lat);
% Now find the nan's
nanLocations = isnan(bot_salt_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(bot_salt_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(bot_salt_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(bot_salt_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = bot_salt_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	bot_salt_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(bot_salt_in);
numberOfNans = sum(nanLocations(:));

%% no3

sur_no3= griddata(lon_obs, lat_obs, squeeze(no3_obs(1,:)),lon_me,lat_me);
sur_no3_in = griddata(lon_obs, lat_obs, sur_no3,lon,lat);
% Now find the nan's
nanLocations = isnan(sur_no3_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(sur_no3_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(sur_no3_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(sur_no3_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = sur_no3_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	sur_no3_in(x,y,z) = goodValue;
end
% u should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(sur_no3_in);
numberOfNans = sum(nanLocations(:));

mid_no3 = griddata(lon_obs, lat_obs, squeeze(no3_obs(2,:)),lon_me,lat_me);
mid_no3_in = griddata(lon_obs, lat_obs, mid_no3,lon,lat);
% Now find the nan's
nanLocations = isnan(mid_no3_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(mid_no3_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(mid_no3_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(mid_no3_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = mid_no3_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	mid_no3_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(mid_no3_in);
numberOfNans = sum(nanLocations(:));

bot_no3 = griddata(lon_obs, lat_obs, squeeze(no3_obs(3,:)),lon_me,lat_me);
bot_no3_in = griddata(lon_obs, lat_obs, bot_no3,lon,lat);
% Now find the nan's
nanLocations = isnan(bot_no3_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(bot_no3_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(bot_no3_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(bot_no3_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = bot_no3_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	bot_no3_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(bot_no3_in);
numberOfNans = sum(nanLocations(:));

%% chla

sur_chla= griddata(lon_obs, lat_obs, squeeze(chla_obs(1,:)),lon_me,lat_me);
sur_chla_in = griddata(lon_obs, lat_obs, sur_chla,lon,lat);
% Now find the nan's
nanLocations = isnan(sur_chla_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(sur_chla_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(sur_chla_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(sur_chla_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = sur_chla_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	sur_chla_in(x,y,z) = goodValue;
end
% u should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(sur_chla_in);
numberOfNans = sum(nanLocations(:));

mid_chla = griddata(lon_obs, lat_obs, squeeze(chla_obs(2,:)),lon_me,lat_me);
mid_chla_in = griddata(lon_obs, lat_obs, mid_chla,lon,lat);
% Now find the nan's
nanLocations = isnan(mid_chla_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(mid_chla_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(mid_chla_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(mid_chla_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = mid_chla_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	mid_chla_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(mid_chla_in);
numberOfNans = sum(nanLocations(:));

bot_chla = griddata(lon_obs, lat_obs, squeeze(chla_obs(3,:)),lon_me,lat_me);
bot_chla_in = griddata(lon_obs, lat_obs, bot_chla,lon,lat);
% Now find the nan's
nanLocations = isnan(bot_chla_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(bot_chla_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(bot_chla_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(bot_chla_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = bot_chla_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	bot_chla_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(bot_chla_in);
numberOfNans = sum(nanLocations(:));


%% NH4

sur_NH4= griddata(lon_obs, lat_obs, squeeze(NH4_obs(1,:)),lon_me,lat_me);
sur_NH4_in = griddata(lon_obs, lat_obs, sur_NH4,lon,lat);
% Now find the nan's
nanLocations = isnan(sur_NH4_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(sur_NH4_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(sur_NH4_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(sur_NH4_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = sur_NH4_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	sur_NH4_in(x,y,z) = goodValue;
end
% u should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(sur_NH4_in);
numberOfNans = sum(nanLocations(:));

mid_NH4 = griddata(lon_obs, lat_obs, squeeze(NH4_obs(2,:)),lon_me,lat_me);
mid_NH4_in = griddata(lon_obs, lat_obs, mid_NH4,lon,lat);
% Now find the nan's
nanLocations = isnan(mid_NH4_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(mid_NH4_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(mid_NH4_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(mid_NH4_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = mid_NH4_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	mid_NH4_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(mid_NH4_in);
numberOfNans = sum(nanLocations(:));

bot_NH4 = griddata(lon_obs, lat_obs, squeeze(NH4_obs(3,:)),lon_me,lat_me);
bot_NH4_in = griddata(lon_obs, lat_obs, bot_NH4,lon,lat);
% Now find the nan's
nanLocations = isnan(bot_NH4_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(bot_NH4_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(bot_NH4_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(bot_NH4_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = bot_NH4_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	bot_NH4_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(bot_NH4_in);
numberOfNans = sum(nanLocations(:));

%% tPO4

sur_tPO4= griddata(lon_obs, lat_obs, squeeze(tPO4_obs(1,:)),lon_me,lat_me);
sur_tPO4_in = griddata(lon_obs, lat_obs, sur_tPO4,lon,lat);
% Now find the nan's
nanLocations = isnan(sur_tPO4_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(sur_tPO4_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(sur_tPO4_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(sur_tPO4_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = sur_tPO4_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	sur_tPO4_in(x,y,z) = goodValue;
end
% u should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(sur_tPO4_in);
numberOfNans = sum(nanLocations(:));

mid_tPO4 = griddata(lon_obs, lat_obs, squeeze(tPO4_obs(2,:)),lon_me,lat_me);
mid_tPO4_in = griddata(lon_obs, lat_obs, mid_tPO4,lon,lat);
% Now find the nan's
nanLocations = isnan(mid_tPO4_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(mid_tPO4_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(mid_tPO4_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(mid_tPO4_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = mid_tPO4_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	mid_tPO4_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(mid_tPO4_in);
numberOfNans = sum(nanLocations(:));

bot_tPO4 = griddata(lon_obs, lat_obs, squeeze(tPO4_obs(3,:)),lon_me,lat_me);
bot_tPO4_in = griddata(lon_obs, lat_obs, bot_tPO4,lon,lat);
% Now find the nan's
nanLocations = isnan(bot_tPO4_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(bot_tPO4_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(bot_tPO4_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(bot_tPO4_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = bot_tPO4_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	bot_tPO4_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(bot_tPO4_in);
numberOfNans = sum(nanLocations(:));

%% DO

sur_DO= griddata(lon_obs, lat_obs, squeeze(DO_obs(1,:)),lon_me,lat_me);
sur_DO_in = griddata(lon_obs, lat_obs, sur_DO,lon,lat);
% Now find the nan's
nanLocations = isnan(sur_DO_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(sur_DO_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(sur_DO_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(sur_DO_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = sur_DO_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	sur_DO_in(x,y,z) = goodValue;
end
% u should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(sur_DO_in);
numberOfNans = sum(nanLocations(:));

mid_DO = griddata(lon_obs, lat_obs, squeeze(DO_obs(2,:)),lon_me,lat_me);
mid_DO_in = griddata(lon_obs, lat_obs, mid_DO,lon,lat);
% Now find the nan's
nanLocations = isnan(mid_DO_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(mid_DO_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(mid_DO_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(mid_DO_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = mid_DO_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	mid_DO_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(mid_DO_in);
numberOfNans = sum(nanLocations(:));

bot_DO = griddata(lon_obs, lat_obs, squeeze(DO_obs(3,:)),lon_me,lat_me);
bot_DO_in = griddata(lon_obs, lat_obs, bot_DO,lon,lat);
% Now find the nan's
nanLocations = isnan(bot_DO_in);
nanLinearIndexes = find(nanLocations)
nonNanLinearIndexes = setdiff(1:numel(bot_DO_in), nanLinearIndexes);
% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(bot_DO_in), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(bot_DO_in), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the var value there.
	goodValue = bot_DO_in(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in var with the good value.
	bot_DO_in(x,y,z) = goodValue;
end
% mid should be fixed now - no nans in it.
% Double check.  Sum of nans should be zero now.
nanLocations = isnan(bot_DO_in);
numberOfNans = sum(nanLocations(:));

%% make interp. data to 3D
temp_ob = NaN(length(lon(:,1)),length(lon(1,:)),3);
temp_ob(:,:,1) = sur_temp_in;
temp_ob(:,:,2) = mid_temp_in;
temp_ob(:,:,3) = bot_temp_in;

salt_ob = NaN(length(lon(:,1)),length(lon(1,:)),3);
salt_ob(:,:,1) = sur_salt_in;
salt_ob(:,:,2) =mid_salt_in;
salt_ob(:,:,3) =bot_salt_in;

NH4_ob = NaN(length(lon(:,1)),length(lon(1,:)),3);
NH4_ob(:,:,1)= sur_NH4_in;
NH4_ob(:,:,2)= mid_NH4_in;
NH4_ob(:,:,3) = bot_NH4_in;

NO3_ob = NaN(length(lon(:,1)),length(lon(1,:)),3);
NO3_ob(:,:,1) = sur_no3_in;
NO3_ob(:,:,2) =mid_no3_in;
NO3_ob(:,:,3) =bot_no3_in;

tPO4_ob = NaN(length(lon(:,1)),length(lon(1,:)),3);
tPO4_ob(:,:,1) = sur_tPO4_in;
tPO4_ob(:,:,2) =mid_tPO4_in;
tPO4_ob(:,:,3) =bot_tPO4_in;

DO_ob = NaN(length(lon(:,1)),length(lon(1,:)),3);
DO_ob(:,:,1) = sur_DO_in;
DO_ob(:,:,2) = mid_DO_in;
DO_ob(:,:,3) = bot_DO_in;

chla_ob = NaN(length(lon(:,1)),length(lon(1,:)),3);
chla_ob(:,:,1) = sur_chla_in;
chla_ob(:,:,2) = mid_chla_in;
chla_ob(:,:,3) = bot_chla_in;

%raw grid (before interped)
for i =1:3
    lon_bf(:,:,i) = lon; 
    lat_bf(:,:,i) = lat; 
end

dep_raw = [0,5,10];
dep_bf = NaN(length(lon(:,1)),length(lon(1,:)),length(dep_raw));
for i = 1:length(dep_raw)
    dep_bf(:,:,i) = dep_raw(i);
end

%% interp 3D (depth)
dep = (0:15); %depth

for i =1:16
    lon_ex(:,:,i) = lon; 
    lat_ex(:,:,i) = lat; 
end

dep_ex = NaN(length(lon(:,1)),length(lon(1,:)),length(dep));
for i = 1:length(dep)
    dep_ex(:,:,i) = dep(i);
end
%%%%%%%%%%%%%% depth expend 10m to 15m (same value filling) %%%%%%%%%%%%%%% 
% TEMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_v_in = griddata(lon_bf,lat_bf, dep_bf,temp_ob,lon_ex,lat_ex,dep_ex); for i=1:16; nan_in(i)=size(find(isnan(temp_v_in(:,:,i))==1),1);end;nan_in   %confirm NaN size
for i = 11:16; temp_v_in(:,:,i) = squeeze(temp_v_in(:,:,10)); end;  
for i=1:16; nan_in(i)=size(find(isnan(temp_v_in(:,:,i))==1),1);end; nan_in   %confirm NaN size

% SALT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
salt_v_in = griddata(lon_bf,lat_bf, dep_bf,salt_ob,lon_ex,lat_ex,dep_ex); for i=1:16; nan_in(i)=size(find(isnan(salt_v_in(:,:,i))==1),1);end; nan_in   %confirm NaN size
for i = 11:16; salt_v_in(:,:,i) = squeeze(salt_v_in(:,:,10)); end; 
for i=1:16; nan_in(i)=size(find(isnan(salt_v_in(:,:,i))==1),1);end; nan_in   %confirm NaN size

% NO3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NO3_v_in = griddata(lon_bf,lat_bf, dep_bf,NO3_ob,lon_ex,lat_ex,dep_ex); for i=1:16; nan_in(i)=size(find(isnan(NO3_v_in(:,:,i))==1),1);end; nan_in   %confirm NaN size
for i = 11:16; NO3_v_in(:,:,i) = squeeze(NO3_v_in(:,:,10)); end; 
for i=1:16; nan_in(i)=size(find(isnan(NO3_v_in(:,:,i))==1),1);end; nan_in   %confirm NaN size

% NH4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NH4_v_in = griddata(lon_bf,lat_bf, dep_bf,NH4_ob,lon_ex,lat_ex,dep_ex); for i=1:16; nan_in(i)=size(find(isnan(NH4_v_in(:,:,i))==1),1);end;nan_in   %confirm NaN size
for i = 11:16; NH4_v_in(:,:,i) = squeeze(NH4_v_in(:,:,10)); end; 
for i=1:16; nan_in(i)=size(find(isnan(NH4_v_in(:,:,i))==1),1);end;nan_in   %confirm NaN size

% tPO4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tPO4_v_in = griddata(lon_bf,lat_bf, dep_bf,tPO4_ob,lon_ex,lat_ex,dep_ex); for i=1:16; nan_in(i)=size(find(isnan(tPO4_v_in(:,:,i))==1),1);end;nan_in   %confirm NaN size
for i = 11:16; tPO4_v_in(:,:,i) = squeeze(tPO4_v_in(:,:,10)); end; 
for i=1:16; nan_in(i)=size(find(isnan(tPO4_v_in(:,:,i))==1),1);end;nan_in   %confirm NaN size

% DO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DO_v_in = griddata(lon_bf,lat_bf, dep_bf,DO_ob,lon_ex,lat_ex,dep_ex); for i=1:16; nan_in(i)=size(find(isnan(DO_v_in(:,:,i))==1),1);end;nan_in   %confirm NaN size
for i = 11:16; DO_v_in(:,:,i) = squeeze(DO_v_in(:,:,10)); end; 
for i=1:16; nan_in(i)=size(find(isnan(DO_v_in(:,:,i))==1),1);end;nan_in   %confirm NaN size

% chla %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chla_v_in = griddata(lon_bf,lat_bf, dep_bf,chla_ob,lon_ex,lat_ex,dep_ex); for i=1:16; nan_in(i)=size(find(isnan(chla_v_in(:,:,i))==1),1);end;nan_in   %confirm NaN size
for i = 11:16; chla_v_in(:,:,i) = squeeze(chla_v_in(:,:,10)); end; 
for i=1:16; nan_in(i)=size(find(isnan(chla_v_in(:,:,i))==1),1);end;nan_in   %confirm NaN size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('interped_variable_fix_po4_fill_const_mean.mat','-v7.3');
%%%%%%%%%%%%%%%%    Make Netcdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
close all; clear all; clc;
load interped_variable_fix_po4_fill_const_mean.mat
lonlon=lon_ex(:,1,1);
latlat=lat_ex(1,:,1);
depdep = dep_ex(1,1,:);
% input format 
temp_v_in = permute(temp_v_in, [3 2 1]);
salt_v_in = permute(salt_v_in, [3 2 1]);
NO3_v_in = permute(NO3_v_in, [3 2 1]);
NH4_v_in = permute(NH4_v_in, [3 2 1]);
DO_v_in = permute(DO_v_in, [3 2 1]);
chla_v_in = permute(chla_v_in, [3 2 1]);
tPO4_v_in = permute(tPO4_v_in, [3 2 1]);
%% temp
nw=netcdf('temp_pohang_ch.nc','clobber');

nw('X') = length(lonlon);
nw('Y') = length(latlat);
% nw('T') = 12;
nw('Z') = length(depdep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=lonlon;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=latlat;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=depdep;

nw{'temperature'} = ncfloat('Z','Y','X');
nw{'temperature'}.long_name = ncchar('temperature');
nw{'temperature'}.units = ncchar('deg.C');
nw{'temperature'}.missing_value=ncfloat(-99.9999008178711);
nw{'temperature'}(:,:,:)=temp_v_in;

close(nw);

%% salt
nw=netcdf('salt_pohang_ch.nc','clobber');

nw('X') = length(lonlon);
nw('Y') = length(latlat);
% nw('T') = 12;
nw('Z') = length(depdep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=lonlon;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=latlat;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=depdep;

nw{'salinity'} = ncfloat('Z','Y','X');
nw{'salinity'}.long_name = ncchar('salinity');
nw{'salinity'}.units = ncchar('salinity');
nw{'salinity'}.missing_value=ncfloat(-99.9999008178711);
nw{'salinity'}(:,:,:)=salt_v_in;

close(nw);

%% NO3
nw=netcdf('NO3_pohang_ch.nc','clobber');

nw('X') = length(lonlon);
nw('Y') = length(latlat);
% nw('T') = 12;
nw('Z') = length(depdep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=lonlon;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=latlat;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=depdep;

nw{'NO3'} = ncfloat('Z','Y','X');
nw{'NO3'}.long_name = ncchar('NO3');
nw{'NO3'}.units = ncchar('micro mole N/L');
nw{'NO3'}.missing_value=ncfloat(-99.9999008178711);
nw{'NO3'}(:,:,:)=NO3_v_in;

close(nw);

%% NH4
nw=netcdf('NH4_pohang_ch.nc','clobber');

nw('X') = length(lonlon);
nw('Y') = length(latlat);
% nw('T') = 12;
nw('Z') = length(depdep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=lonlon;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=latlat;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=depdep;

nw{'NH4'} = ncfloat('Z','Y','X');
nw{'NH4'}.long_name = ncchar('NH4');
nw{'NH4'}.units = ncchar('micro mole N/L');
nw{'NH4'}.missing_value=ncfloat(-99.9999008178711);
nw{'NH4'}(:,:,:)=NH4_v_in;

close(nw);

%% tPO4
nw=netcdf('PO4_pohang_ch_fix_fill_const_mean.nc','clobber');

nw('X') = length(lonlon);
nw('Y') = length(latlat);
% nw('T') = 12;
nw('Z') = length(depdep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=lonlon;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=latlat;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=depdep;

nw{'tPO4'} = ncfloat('Z','Y','X');
nw{'tPO4'}.long_name = ncchar('tPO4');
nw{'tPO4'}.units = ncchar('micro mole P/L');
nw{'tPO4'}.missing_value=ncfloat(-99.9999008178711);
nw{'tPO4'}(:,:,:)=tPO4_v_in;

close(nw);

%% DO
nw=netcdf('DO_pohang_ch.nc','clobber');

nw('X') = length(lonlon);
nw('Y') = length(latlat);
% nw('T') = 12;
nw('Z') = length(depdep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=lonlon;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=latlat;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=depdep;

nw{'oxygen'} = ncfloat('Z','Y','X');
nw{'oxygen'}.long_name = ncchar('dissolved oxygen concentration');
nw{'oxygen'}.units = ncchar('millimole_oxygen meter-3');
nw{'oxygen'}.missing_value=ncfloat(-99.9999008178711);
nw{'oxygen'}(:,:,:)=DO_v_in;

close(nw);

%% chla
nw=netcdf('chla_pohang_ch.nc','clobber');

nw('X') = length(lonlon);
nw('Y') = length(latlat);
% nw('T') = 12;
nw('Z') = length(depdep);

nw{'X'} = ncfloat('X');
nw{'X'}.long_name = ncchar('Longitude');
nw{'X'}.units = ncchar('degree_east');
nw{'X'}(:)=lonlon;

nw{'Y'} = ncfloat('Y');
nw{'Y'}.long_name = ncchar('Latitude');
nw{'Y'}.units = ncchar('degree_north');
nw{'Y'}(:)=latlat;

nw{'Z'} = ncfloat('Z');
nw{'Z'}.long_name = ncchar('Depth');
nw{'Z'}.units = ncchar('meter');
nw{'Z'}(:)=depdep;

nw{'chlorophyll'} = ncfloat('Z','Y','X');
nw{'chlorophyll'}.long_name = ncchar('chlorophyll');
nw{'chlorophyll'}.units = ncchar('milligrams chlorophyll meter-3');
nw{'chlorophyll'}.missing_value=ncfloat(-99.9999008178711);
nw{'chlorophyll'}(:,:,:)=chla_v_in;

close(nw);
