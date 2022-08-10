close all; clear all; clc; 
[raw1 txt1]=xlsread('KOEM_관측위치(deg_min_sec).xlsx','항만','');
[raw2 txt2]=xlsread('KOEM_관측위치(deg_min_sec).xlsx','하천','');
[raw3 txt3]=xlsread('KOEM_관측위치(deg_min_sec).xlsx','연안','');

lat_deg_1=char(txt1(:,2)); lat_deg_1(:,end) = [];
lat_min_1=char(txt1(:,3)); lat_min_1(:,end) = [];
lat_sec_1=char(txt1(:,4)); lat_sec_1(:,end) = [];
lon_deg_1=char(txt1(:,6)); lon_deg_1(:,end) = [];
lon_min_1=char(txt1(:,7)); lon_min_1(:,end) = [];
lon_sec_1=char(txt1(:,8)); lon_sec_1(:,end) = [];

lat_deg_2=char(txt2(:,3)); lat_deg_2(:,end) = [];
lat_min_2=char(txt2(:,4)); lat_min_2(:,end) = [];
lat_sec_2=char(txt2(:,5)); lat_sec_2(:,end) = [];
lon_deg_2=char(txt2(:,7)); lon_deg_2(:,end) = [];
lon_min_2=char(txt2(:,8)); lon_min_2(:,end) = [];
lon_sec_2=char(txt2(:,9)); lon_sec_2(:,end) = [];

lat_deg_3=char(txt3(:,3)); lat_deg_3(:,end) = [];
lat_min_3=char(txt3(:,4)); lat_min_3(:,end) = [];
lat_sec_3=char(txt3(:,5)); lat_sec_3(:,end) = [];
lon_deg_3=char(txt3(:,7)); lon_deg_3(:,end) = [];
lon_min_3=char(txt3(:,8)); lon_min_3(:,end) = [];
lon_sec_3=char(txt3(:,9)); lon_sec_3(:,end) = [];

% Decimal Degrees = degrees + (minutes/60) + (seconds/3600)

lat_1 = str2num(lat_deg_1) + (str2num(lat_min_1)/60) + (str2num(lat_sec_1)/3600);
lon_1 = str2num(lon_deg_1) + (str2num(lon_min_1)/60) + (str2num(lon_sec_1)/3600);
lat_2 = str2num(lat_deg_2) + (str2num(lat_min_2)/60) + (str2num(lat_sec_2)/3600);
lon_2 = str2num(lon_deg_2) + (str2num(lon_min_2)/60) + (str2num(lon_sec_2)/3600);
lat_3 = str2num(lat_deg_3) + (str2num(lat_min_3)/60) + (str2num(lat_sec_3)/3600);
lon_3 = str2num(lon_deg_3) + (str2num(lon_min_3)/60) + (str2num(lon_sec_3)/3600);

clearvars *_deg_* *_min_* *_sec_*
save('KOEM_st_info_(decimal_deg).mat','lon_*','lat_*');