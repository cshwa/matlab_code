clc; clear all; clear all;

% �������� ���� �ڷ��� ��ǥ�� utm �̹Ƿ� �̸� degree�� ��ȯ�ϱ�
% ��ȭ�ϴ� �Լ��� utm2deg.m

data = load('cshwa_fix_land_nolabel_silver.txt');
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

% 
% lat = lat+0.0092;;
% lon = lon-0.0028;



temp(:,2) = lat;
temp(:,1) = lon;
% temp(:,3) = data(:,3);

save gjb_coast_deg_1970.dat temp -ascii




