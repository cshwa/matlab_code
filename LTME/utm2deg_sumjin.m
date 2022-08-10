clc; clear all; clear all;

% 섬진강의 수심 자료의 좌표가 utm 이므로 이를 degree로 변환하기
% 변화하는 함수는 utm2deg.m

data = load('cshwa_fix_merge_silver.txt');
% [Lat,Lon] = utm2deg(x,y,utmzone)
% 대한민국은 UTM 좌표계에서 51S, 51T, 52S, 52T 
% 섬진강은 52S

x = data(:,1);
y = data(:,2);
% x = 827.3;
% y = 1560.3;

l = length(x);
A = ['52 S'];
utmzone = repmat(A, l, 1);
[lat,lon] = utm2deg(x,y,utmzone);
% 1번 파일 수정값;
% lat = lat+34.9518;
% lon = lon+3.25;
% 2번 파일 수정값;
% lat = lat+34.9518+0.0303;
% lon = lon+3.25+0.0086;
% 3번 파일 수정값;
% lat = lat+34.9518+0.0303+0.033;
% lon = lon+3.25+0.0086;
% 4번 파일 수정값;
% lat = lat+34.9518+0.0303+0.033+0.0336;
% lon = lon+3.228+0.0095;
% 5번 파일 수정값;
% lat = lat+34.9518+0.0303+0.033+0.0336+0.0165;
% lon = lon+3.2+0.0125;

% lat = lat+0.0092;
% lon = lon-0.0028;

temp(:,2) = lat;
temp(:,1) = lon;
temp(:,3) = -data(:,3);
save depth_80s_degree.dat temp -ascii

% ss = temp;
% [a b] = find(ss(:,1)<127.56);ss(a,: ) = [];
% [a b] = find(ss(:,1)>128.2);ss(a,: ) = [];
% [a b] = find(ss(:,2)<34.58);ss(a,: ) = [];
% [a b] = find(ss(:,2)>35.12);ss(a,: ) = [];
% save depth_degree_silver.dat ss -ascii




