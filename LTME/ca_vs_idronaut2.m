
clc; clear all; close all;


% temperature
figure() 
subplot(1,14,1)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0047.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f1.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,2)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0048.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f2.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,3)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0049.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f3.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,4)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0050.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f4.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,5)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0051.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f5.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,6)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0052.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f6.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,7)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0053.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f7.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,8)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0054.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f8.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,9)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0055.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f9.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,10)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0056.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f10.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,11)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0057.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f11.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,12)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0058.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f12.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,13)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0059.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f13.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,14)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0060.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f14.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);


% salinity
figure() 
subplot(1,14,1)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0047.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f1.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,2)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0048.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f2.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,3)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0049.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f3.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,4)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0050.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f4.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,5)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0051.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f5.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,6)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0052.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f6.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,7)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0053.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f7.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,8)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0054.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f8.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,9)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0055.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f9.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,10)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0056.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f10.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,11)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0057.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f11.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,12)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0058.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f12.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,13)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0059.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f13.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,14)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\4flood\raw/CAST0060.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\flood/f14.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);


%% ebb ½Ã±â

% temperature
figure() 
subplot(1,14,1)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0016.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e1.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,2)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0017.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e2.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,3)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0018.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e3.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,4)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0019.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e4.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,5)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0020.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e5.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,6)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0021.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e6.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,7)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0023.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e7.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,8)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0031.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e8.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);
subplot(1,14,9)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0024.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e9.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,10)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0025.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e10.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,11)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0026.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e11.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,12)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0027.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e12.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,13)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0028.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e13.mat');
plot(Temperature,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([6 10 0 15]);

subplot(1,14,14)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0029.txt');
plot(data(:,4),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e14.mat');
plot(Temperature,Depth,'r.');
az = 180;
el = 90;
view(az, el);
axis([6 10 0 15]);


% salinity
figure() 
subplot(1,14,1)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0016.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e1.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,2)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0017.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e2.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,3)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0018.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e3.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,4)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0019.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e4.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,5)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0020.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e5.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,6)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0021.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e6.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,7)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0023.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e7.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,8)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0031.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e8.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,9)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0024.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e9.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,10)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0025.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e10.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,11)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0026.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e11.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,12)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0027.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e12.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,13)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0028.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e13.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);

subplot(1,14,14)
data = load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\2ebb\raw/CAST0029.txt');
plot(data(:,6),data(:,3),'b.');
hold on;
load('D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\raw_mat\ebb/e14.mat');
plot(Salinity,Depth,'r.');
az = 0;
el = 270;
view(az, el);
axis([10 35 0 15]);


