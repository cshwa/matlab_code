clc; clear all; close all;
directory = pwd;
filename = 'surface_salinity.xlsx';
filename = fullfile(directory,filename);
S_S =  xlsread(filename);

% figure()
% plot(6:27, data(6:27,1));
% hold on;
% plot(6:30, data(6:30,2));
% plot(6:30, data(6:30,3));
% plot(6:30, data(6:30,4));
% plot(1:24, data(1:24,5));
% plot(1:24, data(1:24,6));
% plot(1:30, data(1:30,7));
% plot(1:30, data(1:30,8));

%%
load 'D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\pm_ebb\neap_ebb_sal.dat';
load 'D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\am_flood\neap_flood_sal.dat';
