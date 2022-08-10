close all; clear; clc; 

[raw txt]=xlsread('D:\장기생태\Dynamic\06_river\광양_진월_수질측정망_최하류.xlsx','Sheet1','');
[rawg txtg]=xlsread('D:\장기생태\Dynamic\06_river\광양_진월_수질측정망_최하류.xlsx','Gure_cut','');
[rawh txth]=xlsread('D:\장기생태\Dynamic\06_river\광양_진월_수질측정망_최하류.xlsx','hadong_cut','');

% Nitrate Nitrogen (NO3-N)
% MW NO3 = 62.005010
% MW N = 14.0067201 
% 1μg NO3/l = 1/ MW NO3 μg/l = 0.016128 μmol NO3/l
% 1 μg NO3/l = MW N/MW NO3 = 0.225897 μg N/l
% 1 μg N/l = 1/MW N = 0.071394 μmol N/l

% Ammonium Nitrogen (NH4-N)
% MW NH4 = 18.038508
% MW N = 14.0067201 
% 1μg NH4/l = 1/ MW NH4 = 0.055437 μmol NH4/l
% 1 μg NH4/l = MW N/MW NH4 = 0.776490 μg N/l
% 1 μg N/l = 1/MW N = 0.071394 μmol N/l

% Oxygen (O2)
% Molar volume at STP = 22.391 l
% Molar weight of oxygen = 31.998 g
% Atomic Mass of oxygen = 15.994 g/mol
% 1 μmol O2 = .022391 ml
% 1 ml/l = 103/22.391 = 44.661 μmol/l
% 1 mg/l = 22.391 ml/31.998 = 0.700 ml/l
% 1 mg-at/l = 15.994x22.391/31.998 = 11.192 ml

% 9 : NO3-N mg/L
% 1 mg/L to 1000 ug/L 

% no3_plot=raw(:,9).*1000 ./14; % 1~12: 2007
% 21~32: 2016
% 33~44: 2017
no3_plot=raw(:,9).*1000 .* 0.071394; %  0.071394 = 1/14
nh4_plot=raw(:,8).*1000 .* 0.071394; %  0.071394 = 1/14
         
% gure_no3_plot= rawg(:,9).*1000 ./14;          
gure_no3_plot= rawg(:,9).*1000 .* 0.071394;
gure_nh4_plot= rawg(:,8).*1000 .* 0.071394;

hadong_no3_plot= rawh(:,9).*1000 .* 0.071394;
hadong_nh4_plot= rawh(:,8).*1000 .* 0.071394;

% 1 mg/l (ppm) NH3-N = 1 mg/l (ppm) NH4-N
% 1 mg/l (ppm) NH4 = 0.94413 mg/l (ppm) NH3


figure;
plot(no3_plot(1:12),'r'); hold on;
plot(no3_plot(21:32),'g')
plot(no3_plot(33:44),'b')
plot(gure_no3_plot(1:12),'r-*'); hold on;
plot(gure_no3_plot(13:24),'g-*')
plot(gure_no3_plot(25:36),'b-*')
plot(hadong_no3_plot(1:12),'r-d'); hold on;
plot(hadong_no3_plot(13:24),'g-d')
plot(hadong_no3_plot(25:36),'b-d')
xlim([1 12]);
xlabel('시간 (월)','fontsize',13)
ylabel('NO3-N (mmol N / m^3)','fontsize',13)
set(gca,'xtick',[1:12]);
% set(gca,'xlim',[1 22]);
% set(gca,'xticklabel',1997:2:2018);
title('진월(하류) vs. 하동(중) vs. 구례(상류) NO3-N 농도 비교','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('진월-07','진월-16','진월-17','구례-07','구례-16','구례-17')
le = legend('진월-07','진월-16','진월-17','구례-07','구례-16',...
    '구례-17','하동-07','하동-16','하동-17');
set(le,'fontsize',8)


figure;
plot(nh4_plot(1:12),'r'); hold on;
plot(nh4_plot(21:32),'g')
plot(nh4_plot(33:44),'b')
plot(gure_nh4_plot(1:12),'r-*'); hold on;
plot(gure_nh4_plot(13:24),'g-*')
plot(gure_nh4_plot(25:36),'b-*')
plot(hadong_nh4_plot(1:12),'r-d'); hold on;
plot(hadong_nh4_plot(13:24),'g-d')
plot(hadong_nh4_plot(25:36),'b-d')
xlim([1 12]);
xlabel('시간 (월)','fontsize',13)
ylabel('NH4-N (mmol N / m^3)','fontsize',13)
set(gca,'xtick',[1:12]);
% set(gca,'xlim',[1 22]);
% set(gca,'xticklabel',1997:2:2018);
title('진월(하) vs. 하동(중) vs. 구례(상) NH4-N 농도 비교','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('진월-07','진월-16','진월-17','구례-07','구례-16','구례-17')
le = legend('진월-07','진월-16','진월-17','구례-07','구례-16',...
    '구례-17','하동-07','하동-16','하동-17');
set(le,'fontsize',8)

