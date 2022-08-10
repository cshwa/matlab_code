
clc; clear all; close all;

s02 = load('step02_coastline_utm_with_nan.txt');
lon = s02(:,1); lat = s02(:,2);
%--- utm 2 degree -----
data = s02;
x = data(:,1);
y = data(:,2);

l = length(x);
A = ['52 S'];
utmzone = repmat(A, l, 1);
[lat,lon] = utm2deg(x,y,utmzone);

temp(:,2) = lat;
temp(:,1) = lon;

%--- degree to bln --------
k = 1;
j = 0;
for i =1:length(temp)
    j = 1+1;
    if isnan(temp(i,1))
        i
        no = i
%         if k == 1
            ss(k,1) = j-1;
            ss(k,2) = 1;
            ss(k,3) = 99999;
            ss(k+1:i+k,1) = temp(k:i,1);
            ss(k+1:i+k,2) = temp(k:i,2);
            ss(k+1:i+k,3) = 0;
            k = k+i+1;        
%         elseif k > 2
%             ss(k,1) = i;
%             ss(k,2) = 1;
%             ss(k,3) = 99999;
%             ss(k+1:i+k,1) = temp(k:i,1);
%             ss(k+1:i+k,2) = temp(k:i,2);
%             ss(k+1:i+k,3) = 0;
%             k = k+i  
%         end
    end
end
    