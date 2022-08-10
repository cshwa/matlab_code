clc; clear all; close all

%% 과거 자료에서 새로 매립된 곳 추가하여 수심 변경
data=load('..\depth_from_yjg\depth.xyz');

lon=data(:,1);
lat=data(:,2);
zz=data(:,3);

data=load('..\depth_from_yjg\2012depth.xyz');
lon1=data(:,1);
lat1=data(:,2);
z1=data(:,3);

x=[lon; lon1];
y=[lat; lat1];
z=[zz; z1];


l = length(x);
A = ['52 S'];
utmzone = repmat(A, l, 1);
[lat,lon] = utm2deg(x,y,utmzone);

temp(:,2) = lat;
temp(:,1) = lon;
temp(:,3) = z;


save gjb_dep_deg.dat temp -ascii

save('gjb_dep_deg.mat','lon','lat','z')

