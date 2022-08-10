clc; clear all; close all

temp = load('location_gy.dat');
lon = temp(:,3);
lat = temp(:,2);
clear temp;

filename = 'ctd_gy.xls';
sheet = 1;
[num,txt,raw] = xlsread(filename,sheet);

% 원하는 년도: syear
syear = 2012;
[a b] = find(num(:,2) == syear);

data = num(a,:);