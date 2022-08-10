clc; clear all; clear all;

% �������� ���� �ڷ��� ��ǥ�� utm �̹Ƿ� �̸� degree�� ��ȯ�ϱ�
% ��ȭ�ϴ� �Լ��� utm2deg.m

data = load('cshwa_fix_merge_silver.txt');
% [Lat,Lon] = utm2deg(x,y,utmzone)
% ���ѹα��� UTM ��ǥ�迡�� 51S, 51T, 52S, 52T 
% �������� 52S

x = data(:,1);
y = data(:,2);
% x = 827.3;
% y = 1560.3;

l = length(x);
A = ['52 S'];
utmzone = repmat(A, l, 1);
[lat,lon] = utm2deg(x,y,utmzone);
% 1�� ���� ������;
% lat = lat+34.9518;
% lon = lon+3.25;
% 2�� ���� ������;
% lat = lat+34.9518+0.0303;
% lon = lon+3.25+0.0086;
% 3�� ���� ������;
% lat = lat+34.9518+0.0303+0.033;
% lon = lon+3.25+0.0086;
% 4�� ���� ������;
% lat = lat+34.9518+0.0303+0.033+0.0336;
% lon = lon+3.228+0.0095;
% 5�� ���� ������;
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




