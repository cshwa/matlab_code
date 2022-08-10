% 섬진강 adcp mooring
% N34 59 13.3 
% E127 46 31.8
% 관측간격 1m
% First cell 1.5m
% 관측 시간 간격: 10 min
% 사용한 자료: 2015-6-17-06:00 - 2015-3-06-16:19

clc; clear all; close all;

% mat 형식 데이터 load
load all.mat

st = 1; en = 2929;
pg1 = SerPG4([st:en],:);

% pg1으로 관측 시작/끝 시점 찾기
% 100과 2400은 시작과 끝 지점에서 어느정도 떨어진 값으로 임의로 잡음
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