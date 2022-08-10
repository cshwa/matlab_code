% ������ adcp mooring
% N34 59 13.3 
% E127 46 31.8
% �������� 1m
% First cell 1.5m
% ���� �ð� ����: 10 min
% ����� �ڷ�: 2015-6-17-06:00 - 2015-3-06-16:19

clc; clear all; close all;

% mat ���� ������ load
load all.mat

st = 1; en = 2929;
pg1 = SerPG4([st:en],:);

% pg1���� ���� ����/�� ���� ã��
% 100�� 2400�� ���۰� �� �������� ������� ������ ������ ���Ƿ� ����
st = min(find(pg1(1:100,1) >= 99));
en = min(find(pg1(2400:end,1) < 99))+2400-2;

depth = AnDepthmm([st:en],:); % unit: mm
mag = SerMagmmpersec([st:en],:);
dir = SerDir10thDeg([st:en],:)/10;
v = SerNmmpersec([st:en],:);
u = SerEmmpersec([st:en],:);
w = SerVmmpersec([st:en],:);
pg1 = SerPG4([st:en],:);
pg4 = SerPG4([st:en],:);
error = SerErmmpersec([st:en],:);
coabg = SerCAcnt([st:en],:);
echoavg = SerEAAcnt([st:en],:);
dd = depth/1000;
depth = floor(depth/1000);  % unit: m

save dd.txt dd -ascii