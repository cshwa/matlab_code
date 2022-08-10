
clc; clear all; close all;

load coastline_80s.mat

% depth=load('cshwa_fix_merge_silver.txt');
s02 = load('step02_coastline_utm_with_nan.txt');
lon = s02(:,1); lat = s02(:,2);
% save coastline_80s.mat lon lat 
%--- utm 2 degree -----
data = s02;
x = data(:,1);
y = data(:,2);

l = length(x);
A = ['52 S'];
utmzone = repmat(A, l, 1);
[lat,lon] = utm2deg(x,y,utmzone);

% save coastline_80s.mat lon lat

temp(:,2) = lat;
temp(:,1) = lon;
temp(:,3) = data(:,3);


%--- degree to bln --------
% k = 0;
% for i =1:length(s02)
%     if ~isnan(s02(i,1))
%         k = k+1;
%     else 
%         s02(i,1) = k;
%     end
% end
    