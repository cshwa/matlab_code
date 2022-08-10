%% merge discharge
% close all; clear all;
% load('D:\장기생태\Dynamic\06_river\data\sj_1980to1996\songjung_discharge_1980to2018.mat');  % pre_merg_dis is discharge
% 
% % load excel 2019 trans
% [raw txt]=xlsread('D:\장기생태\Dynamic\06_river\환경과학원\수리·수문·기상_유량_구례_20190101-20191231.xls','유량');
% [raw1 txt1]=xlsread('D:\장기생태\Dynamic\06_river\환경과학원\수리·수문·기상_유량_구례_20200101-20201231.xls','유량');
% 
% t_year = 2019:2020
% 
% order_i = length(1980:2018);
% 
% dis_pre_total{order_i+1} = raw(:,4);
% dis_pre_total{order_i+2} = [raw1(1:59,4); raw1(59,4); raw1(60:end,4)]; % 2020 songjumg discharge has 365 but 2020 is leap year 
% 
% save('songjung_discharge_1980to2020.mat')

close all; clear; clc;
% -------------------------------------------------------------------------
t_year = 1982:2018
% _________________________________________________________________________
% cd D:\장기생태\Dynamic\06_river
% load jinwall_raw_data_daily.mat %there has to be song vs. song_transp
%jinwall vs. jinwall_transp
% load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','raw_*')
cd D:\장기생태\Dynamic\06_river\환경과학원
% load songjung_yymm_koem_data_89to19_3sig.mat  % polynomial_fitting_on_ecological_variable_to2004_new_v6_3sigma.m
load songjung_yymm_monthly_data_to06to15_3sig_1993to2004.mat  % polynomial_fitting_on_ecological_variable_to06to15_new_v7_3sigma_1993to2004.m
% addP=load('songjung_bio_v3_day_koem_yr.mat'); % 2004_regime_v3.m

cd D:\장기생태\Dynamic\06_river
yj=load('koem_timeseires_monthly_gy_only_1to3sig_v2_each_st_16p.mat'); % 16 points 

cd D:\장기생태\Dynamic\06_river
clearvars merg_recon_w_c
% sjtemp=load('D:\장기생태\Dynamic\06_river_physics\sumjin_recons_water_temp_present_2019.mat');
% sjtemp=load('D:\장기생태\Dynamic\06_river_physics\환경과학원\sumjin_recons_water_temp(agyang)_present_2019_high_airT.mat');
sjtemp=load('D:\장기생태\Dynamic\06_river_physics\환경과학원\sumjin_recons_water_temp(agyang)_present_1980to2020_high_airT_mixed_daily.mat');
sumjin_re_w_c = sjtemp.merg_recon_w_c;

cd D:\장기생태\Dynamic\06_river\data\sj_1980to1996
%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
load songjung_discharge_1980to2018.mat  % pre_merg_dis is discharge

in_date = yj.ref_date;

% load excel 2019 trans
[raw txt]=xlsread('D:\장기생태\Dynamic\06_river\환경과학원\수리·수문·기상_유량_구례_20190101-20191231.xls','유량');
[raw1 txt1]=xlsread('D:\장기생태\Dynamic\06_river\환경과학원\수리·수문·기상_유량_구례_20200101-20201231.xls','유량');

clearvars tempo_trans sj_trans_out_365
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
% tempo_temp=sumjin_re_w_c{t_year(order_i)-1989};
tempo_trans = dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
%     temperature=tempo_temp;
    sj_trans_out_365=tempo_trans;
else
%     temperature=cat(1,temperature,tempo_temp);
    if length(tempo_trans) == 366
        tempo_trans(60)=[];
    end
    sj_trans_out_365=cat(1,sj_trans_out_365,tempo_trans);
end
end
sj_trans_out_365(end+1:end+365) = raw(:,4);  % 2019 songjung discharge
% sj_trans_out_365(end+1:end+365 + leapyear(2020)) = [raw1(1:59,4); raw1(59,4); raw1(60:end,4)]; % 2020 songjumg discharge has 365 but 2020 is leap year 
sj_trans_out_365(end+1:end+365) = raw1(:,4);

clearvars tempo_trans sj_trans_out
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
% tempo_temp=sumjin_re_w_c{t_year(order_i)-1989};
tempo_trans = dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
%     temperature=tempo_temp;
    sj_trans_out=tempo_trans;
else
%     temperature=cat(1,temperature,tempo_temp);
    sj_trans_out=cat(1,sj_trans_out,tempo_trans);
end
end
sj_trans_out(end+1:end+365) = raw(:,4); % 2019 songjung discharge
sj_trans_out(end+1:end+365 + leapyear(2020)) = [raw1(1:59,4); raw1(59,4); raw1(60:end,4)]; % 2020 songjumg discharge has 365 but 2020 is leap year 
% sj_trans_out(end+1:end+365) = raw(:,4); % 2020 songjung discharge

clearvars tempo_temp riv_temp riv_temp_pre
order_i = 0
for i = t_year(1):t_year(end) % 1980~2018
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_temp=sumjin_re_w_c{t_year(order_i)-1979};
if order_i == 1
    riv_temp=tempo_temp;
else
    riv_temp=cat(1,riv_temp,tempo_temp);
end
end
riv_temp(end+1:end+365) = sumjin_re_w_c{40}; % 2019
temp_temp=sumjin_re_w_c{41};
% temp_temp(60)=[];
riv_temp(end+1:end+365 + leapyear(2020)) = temp_temp; % 2020
% riv_temp form it to 1989 ~ 2020
% riv_temp = NaN(1,length(sj_trans_out));
% riv_temp(1,365+leapyear(1989)+1:length(sj_trans_out)) = riv_temp_pre;

% riv_temp(end+1:end+365) = temp_temp; % 2020
 
 clearvars tempo_temp riv_temp_365 riv_temp_365_pre
order_i = 0
for i = t_year(1):t_year(end) % 1990~2018
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_temp=sumjin_re_w_c{t_year(order_i)-1979};
if order_i == 1
    riv_temp_365=tempo_temp;
else
     if length(tempo_temp) == 366
        tempo_temp(60)=[];
    end
    riv_temp_365=cat(1,riv_temp_365,tempo_temp);
end
end
riv_temp_365(end+1:end+365) = sumjin_re_w_c{40}; % 2019
temp_temp=sumjin_re_w_c{41};
temp_temp(60)=[];
riv_temp_365(end+1:end+365) = temp_temp; % 2020
 % riv_temp form it to 1989 ~ 2020
% riv_temp_365 = NaN(1,length(sj_trans_out_365));
% riv_temp_365(1,365+leapyear(1989)+1:length(sj_trans_out_365)) = riv_temp_365_pre;
 

%% monthly mean (yy.mm) form after over 3sigma value extraction.
for i=1:length(t_indx) % time for num. of months (until 2019.12.31)
if i ==1
   monthly_trans(i) = nanmean(sj_trans_out(1:t_indx(i)));
else
   monthly_trans(i) = nanmean(sj_trans_out(t_indx(i-1)+1:t_indx(i)));
end
end

riv_temp_t_ind=t_indx; % 1982~2020
k=0;
clearvars monthly_riv_temp
% for i=13:length(t_indx) % time for num. of months (until 2019.12.31)
for i=1:length(t_indx) % time for num. of months (1990.01,01 ~ 2020.12.31)
    k=k+1
    if k == 1
%         monthly_riv_temp(k) = nanmean(riv_temp(1:riv_temp_t_ind(i)));
          monthly_riv_temp(k) = nanmean(riv_temp(1:riv_temp_t_ind(i)));
    else
%         monthly_riv_temp(k) = nanmean(riv_temp(riv_temp_t_ind(i-1)+1:riv_temp_t_ind(i)));
           monthly_riv_temp(k) = nanmean(riv_temp(riv_temp_t_ind(i-1)+1:riv_temp_t_ind(i)));
    end
end

%% 1982 ~ 1986 : 1st
%% 1987 ~ 1996 : 2nd
%% 1997 ~ 2006 : 3rd
%% 2007 ~ 2015 : 4th
%% 2016 ~ 2020 : 5th

% t_day_cut1=365*length(1982:1986);
% t_day_cut2=365*length(1982:1996);
% t_day_cut3=365*length(1982:2006);
% t_day_cut4=365*length(1982:2015);

% sj_day_trans_1=sj_trans_out_365(1:t_day_cut1);
% sj_day_trans_2=sj_trans_out_365(t_day_cut1+1:t_day_cut2);
% sj_day_trans_3=sj_trans_out_365(t_day_cut2+1:t_day_cut3);
% sj_day_trans_4=sj_trans_out_365(t_day_cut3+1:t_day_cut4);
% sj_day_trans_5=sj_trans_out_365(t_day_cut4+1:end);
% 
% sj_day_riv_temp_1=riv_temp_365(1:t_day_cut1);
% sj_day_riv_temp_2=riv_temp_365(t_day_cut1+1:t_day_cut2);
% sj_day_riv_temp_3=riv_temp_365(t_day_cut2+1:t_day_cut3);
% sj_day_riv_temp_4=riv_temp_365(t_day_cut3+1:t_day_cut4);
% sj_day_riv_temp_5=riv_temp_365(t_day_cut4+1:end);


%% 1982 ~ 1986 : 1st
%% 1987 ~ 1992 : 2nd
%% 1993 ~ 2006 : 3rd
%% 2007 ~ 2015 : 4th
%% 2016 ~ 2020 : 5th

t_day_cut1=365*length(1982:1986);
t_day_cut2=365*length(1982:1992);
t_day_cut3=365*length(1982:2004);
t_day_cut4=365*length(1982:2006);
t_day_cut5=365*length(1982:2015);

sj_day_trans_1=sj_trans_out_365(1:t_day_cut1);
sj_day_trans_2=sj_trans_out_365(t_day_cut1+1:t_day_cut2);
sj_day_trans_3=sj_trans_out_365(t_day_cut2+1:t_day_cut3);
sj_day_trans_4=sj_trans_out_365(t_day_cut4+1:t_day_cut5);
sj_day_trans_5=sj_trans_out_365(t_day_cut5+1:end);

sj_day_riv_temp_1=riv_temp_365(1:t_day_cut1);
sj_day_riv_temp_2=riv_temp_365(t_day_cut1+1:t_day_cut2);
sj_day_riv_temp_3=riv_temp_365(t_day_cut2+1:t_day_cut3);
sj_day_riv_temp_4=riv_temp_365(t_day_cut4+1:t_day_cut5);
sj_day_riv_temp_5=riv_temp_365(t_day_cut5+1:end);

% sj_day_riv_temp_1=riv_temp_365(t_day_cut2-365+1:t_day_cut1-365);
% sj_day_riv_temp_2=riv_temp_365(t_day_cut1-365+1:t_day_cut2-365);
% sj_day_riv_temp_3=riv_temp_365(t_day_cut2-365+1:end);

for i = 1:365 %days  
clim_day_trans_1(i) = nanmean(sj_day_trans_1(i:365:end));
clim_day_trans_2(i) = nanmean(sj_day_trans_2(i:365:end));
clim_day_trans_3(i) = nanmean(sj_day_trans_3(i:365:end));
clim_day_trans_4(i) = nanmean(sj_day_trans_4(i:365:end));
clim_day_trans_5(i) = nanmean(sj_day_trans_5(i:365:end));
clim_day_riv_temp_1(i) = nanmean(sj_day_riv_temp_1(i:365:end));
clim_day_riv_temp_2(i) = nanmean(sj_day_riv_temp_2(i:365:end));
clim_day_riv_temp_3(i) = nanmean(sj_day_riv_temp_3(i:365:end));
clim_day_riv_temp_4(i) = nanmean(sj_day_riv_temp_4(i:365:end));
clim_day_riv_temp_5(i) = nanmean(sj_day_riv_temp_5(i:365:end)); 
end

%% additional check for extream discharge
% for i = 1:365 %days  
% clim_day_cum_trans_1(i) = nanmean(sj_day_trans_1(i:365:end)).*86400;
% clim_day_cum_trans_2(i) = nanmean(sj_day_trans_2(i:365:end)).*86400;
% clim_day_cum_trans_3(i) = nanmean(sj_day_trans_3(i:365:end)).*86400;
% % clim_day_riv_temp_1(i) = nanmean(sj_day_riv_temp_1(i:365:end)).*86400;
% % clim_day_riv_temp_2(i) = nanmean(sj_day_riv_temp_2(i:365:end)).*86400;
% % clim_day_riv_temp_3(i) = nanmean(sj_day_riv_temp_3(i:365:end)).*86400; 
% end
% 
% for i = 1:365
% clim_day_cumsum_trans_1(i) = sum(clim_day_cum_trans_1(1:i));
% clim_day_cumsum_trans_2(i) = sum(clim_day_cum_trans_2(1:i));
% clim_day_cumsum_trans_3(i) = sum(clim_day_cum_trans_3(1:i));
% end
%%

%% already cut it out on 'polynomial_fitting_on_ecological_variable_to06to15_new_v6_3sigma.m'
t_s = 8*12; %cut to 1996
% monthly_do(1:t_s) = [];
% monthly_chl(1:t_s) = [];
% monthly_nh4(1:t_s) = [];
% monthly_no3(1:t_s) = [];
% monthly_po4(1:t_s) = [];
% monthly_ss(1:t_s) = [];
monthly_trans(1:t_s) = [];

%% cut 3 regime
% cut_1 = 120 % 2006-12
cut_0 = 120 % 2006-12
cut_1 = 120 % 2006-12
cut_2 = 228 % 2015-12

mon_trans_1 = monthly_trans(1:cut_1);
mon_trans_2 = monthly_trans(cut_1+1:cut_2);
mon_trans_3 = monthly_trans(cut_2+1:end);

mon_do_1 = monthly_do;
% mon_do_2 = monthly_do(cut_1+1:cut_2);
% mon_do_3 = monthly_do(cut_2+1:end);

mon_chl_1 = monthly_chl;
% mon_chl_2 = monthly_chl(cut_1+1:cut_2);
% mon_chl_3 = monthly_chl(cut_2+1:end);

mon_nh4_1 = monthly_nh4;
% mon_nh4_2 = monthly_nh4(cut_1+1:cut_2);
% mon_nh4_3 = monthly_nh4(cut_2+1:end);

mon_no3_1 = monthly_no3;
% mon_no3_2 = monthly_no3(cut_1+1:cut_2);
% mon_no3_3 = monthly_no3(cut_2+1:end);

mon_po4_1 = monthly_po4;
% mon_po4_2 = monthly_po4(cut_1+1:cut_2);
% mon_po4_3 = monthly_po4(cut_2+1:end);

mon_ss_1 = monthly_ss;
% mon_ss_2 = monthly_ss(cut_1+1:cut_2);
% mon_ss_3 = monthly_ss(cut_2+1:end);

mon_tn_1 = monthly_tn;
% mon_tn_2 = monthly_tn(cut_1+1:cut_2);
% mon_tn_3 = monthly_tn(cut_2+1:end);

mon_tp_1 = monthly_tp;
% mon_tp_2 = monthly_tp(cut_1+1:cut_2);
% mon_tp_3 = monthly_tp(cut_2+1:end);


for i = 1:12 %month   
mon_clim_trans_1(i) = nanmean(mon_trans_1(i:12:end));
mon_clim_trans_2(i) = nanmean(mon_trans_2(i:12:end));
mon_clim_trans_3(i) = nanmean(mon_trans_3(i:12:end)); 

mon_clim_do_1(i) = nanmean(mon_do_1(i:12:end));
% mon_clim_do_2(i) = nanmean(mon_do_2(i:12:end));
% mon_clim_do_3(i) = nanmean(mon_do_3(i:12:end));
    
mon_clim_chl_1(i) = nanmean(mon_chl_1(i:12:end));
% mon_clim_chl_2(i) = nanmean(mon_chl_2(i:12:end));
% mon_clim_chl_3(i) = nanmean(mon_chl_3(i:12:end));

mon_clim_nh4_1(i) = nanmean(mon_nh4_1(i:12:end));
% mon_clim_nh4_2(i) = nanmean(mon_nh4_2(i:12:end));
% mon_clim_nh4_3(i) = nanmean(mon_nh4_3(i:12:end));

mon_clim_no3_1(i) = nanmean(mon_no3_1(i:12:end));
% mon_clim_no3_2(i) = nanmean(mon_no3_2(i:12:end));
% mon_clim_no3_3(i) = nanmean(mon_no3_3(i:12:end));

mon_clim_po4_1(i) = nanmean(mon_po4_1(i:12:end));
% mon_clim_po4_2(i) = nanmean(mon_po4_2(i:12:end));
% mon_clim_po4_3(i) = nanmean(mon_po4_3(i:12:end));

mon_clim_ss_1(i) = nanmean(mon_ss_1(i:12:end));
% mon_clim_ss_2(i) = nanmean(mon_ss_2(i:12:end));
% mon_clim_ss_3(i) = nanmean(mon_ss_3(i:12:end));

mon_clim_tn_1(i) = nanmean(mon_tn_1(i:12:end));
% mon_clim_tn_2(i) = nanmean(mon_tn_2(i:12:end));
% mon_clim_tn_3(i) = nanmean(mon_tn_3(i:12:end));

mon_clim_tp_1(i) = nanmean(mon_tp_1(i:12:end));
% mon_clim_tp_2(i) = nanmean(mon_tp_2(i:12:end));
% mon_clim_tp_3(i) = nanmean(mon_tp_3(i:12:end));
end

figure; hold on;
plot(clim_day_trans_1,'r'); plot(clim_day_trans_2,'g'); plot(clim_day_trans_3,'b');
plot(clim_day_trans_4,'m'); plot(clim_day_trans_5,'c');
grid on;
legend('1982~1986','1987~1996','1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_do_1,'r');  grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_chl_1,'r');  grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_nh4_1,'r');  grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_no3_1,'r');  grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_po4_1,'r');  grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_ss_1,'r');  grid on;
legend('1997~2006','2007~2015','2016~2020');


t_indx(1:12) -15
day_indx=[1; 45; 75; 105; 136; 166; 197; 228; 258; 289; 319; 365;];

clearvars day_clim_*
day_clim_do_1 = NaN(365,1);

day_clim_chl_1 = NaN(365,1);

day_clim_nh4_1 = NaN(365,1);

day_clim_no3_1 = NaN(365,1);

day_clim_po4_1 = NaN(365,1);

day_clim_ss_1 = NaN(365,1);

day_clim_tn_1 = NaN(365,1);

day_clim_tp_1 = NaN(365,1);


%%
day_clim_do_1(day_indx) = mon_clim_do_1;

day_clim_chl_1(day_indx) = mon_clim_chl_1;

day_clim_nh4_1(day_indx) = mon_clim_nh4_1;

day_clim_no3_1(day_indx) = mon_clim_no3_1;

day_clim_po4_1(day_indx) = mon_clim_po4_1;

day_clim_ss_1(day_indx) = mon_clim_ss_1;

day_clim_tn_1(day_indx) = mon_clim_tn_1;

day_clim_tp_1(day_indx) = mon_clim_tp_1;


%%
t=1:length(day_clim_do_1);
day_clim_do_1(isnan(day_clim_do_1)) = interp1(t(~isnan(day_clim_do_1)),day_clim_do_1(~isnan(day_clim_do_1)),t(isnan(day_clim_do_1)));

t=1:length(day_clim_chl_1);
day_clim_chl_1(isnan(day_clim_chl_1)) = interp1(t(~isnan(day_clim_chl_1)),day_clim_chl_1(~isnan(day_clim_chl_1)),t(isnan(day_clim_chl_1)));

t=1:length(day_clim_nh4_1);
day_clim_nh4_1(isnan(day_clim_nh4_1)) = interp1(t(~isnan(day_clim_nh4_1)),day_clim_nh4_1(~isnan(day_clim_nh4_1)),t(isnan(day_clim_nh4_1)));

t=1:length(day_clim_no3_1);
day_clim_no3_1(isnan(day_clim_no3_1))  = interp1(t(~isnan(day_clim_no3_1)),day_clim_no3_1(~isnan(day_clim_no3_1)),t(isnan(day_clim_no3_1)));

t=1:length(day_clim_po4_1);
day_clim_po4_1(isnan(day_clim_po4_1)) = interp1(t(~isnan(day_clim_po4_1)),day_clim_po4_1(~isnan(day_clim_po4_1)),t(isnan(day_clim_po4_1)));

t=1:length(day_clim_ss_1);
day_clim_ss_1(isnan(day_clim_ss_1)) = interp1(t(~isnan(day_clim_ss_1)),day_clim_ss_1(~isnan(day_clim_ss_1)),t(isnan(day_clim_ss_1)));

t=1:length(day_clim_tn_1);
day_clim_tn_1(isnan(day_clim_tn_1)) = interp1(t(~isnan(day_clim_tn_1)),day_clim_tn_1(~isnan(day_clim_tn_1)),t(isnan(day_clim_tn_1)));

t=1:length(day_clim_tp_1);
day_clim_tp_1(isnan(day_clim_tp_1)) = interp1(t(~isnan(day_clim_tp_1)),day_clim_tp_1(~isnan(day_clim_tp_1)),t(isnan(day_clim_tp_1)));

plt_filter_DJF = [1:31, 32:59, 335:365];

save('songjung_climate_days_bio_1982to2020_high_airT_mixed_daily_v2.mat','clim_day_*','day_clim_*');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NAMGANG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
% -------------------------------------------------------------------------
t_year = 1982:2018
% _________________________________________________________________________
% cd D:\장기생태\Dynamic\06_river
% load jinwall_raw_data_daily.mat %there has to be song vs. song_transp
%jinwall vs. jinwall_transp
% load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','raw_*')
cd D:\장기생태\Dynamic\06_river\환경과학원
% load songjung_yymm_koem_data_89to19_3sig.mat  % polynomial_fitting_on_ecological_variable_to2004_new_v6_3sigma.m
load namgang_yymm_monthly_data_to06to15_3sig_1993to2004.mat  % polynomial_fitting_on_ecological_variable_to06to15_new_v7_3sigma_1993to2004.m
% addP=load('songjung_bio_v3_day_koem_yr.mat'); % 2004_regime_v3.m

cd D:\장기생태\Dynamic\06_river
yj=load('koem_timeseires_monthly_gy_only_1to3sig_v2_each_st_16p.mat'); % 16 points 

cd D:\장기생태\Dynamic\06_river
clearvars merg_recon_w_c
ghtemp=load('D:\장기생태\Dynamic\06_river_physics\gawha_recons_water_temp_present_1980to2020_high_airT_mixed_daily.mat');
gawha_re_w_c = ghtemp.merg_recon_w_c;

cd D:\장기생태\Dynamic\06_river\data\sj_1980to1996
%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
load songjung_discharge_1980to2018.mat  % pre_merg_dis is discharge

in_date = yj.ref_date;

% load excel 2019 trans
[raw txt]=xlsread('D:\장기생태\Dynamic\06_river\환경과학원\수리·수문·기상_유량_구례_20190101-20191231.xls','유량');
[raw1 txt1]=xlsread('D:\장기생태\Dynamic\06_river\환경과학원\수리·수문·기상_유량_구례_20200101-20201231.xls','유량');

clearvars tempo_trans sj_trans_out_365
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
% tempo_temp=sumjin_re_w_c{t_year(order_i)-1989};
tempo_trans = dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
%     temperature=tempo_temp;
    sj_trans_out_365=tempo_trans;
else
%     temperature=cat(1,temperature,tempo_temp);
    if length(tempo_trans) == 366
        tempo_trans(60)=[];
    end
    sj_trans_out_365=cat(1,sj_trans_out_365,tempo_trans);
end
end
sj_trans_out_365(end+1:end+365) = raw(:,4);  % 2019 songjung discharge
% sj_trans_out_365(end+1:end+365 + leapyear(2020)) = [raw1(1:59,4); raw1(59,4); raw1(60:end,4)]; % 2020 songjumg discharge has 365 but 2020 is leap year 
sj_trans_out_365(end+1:end+365) = raw1(:,4);

clearvars tempo_trans sj_trans_out
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
% tempo_temp=sumjin_re_w_c{t_year(order_i)-1989};
tempo_trans = dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
%     temperature=tempo_temp;
    sj_trans_out=tempo_trans;
else
%     temperature=cat(1,temperature,tempo_temp);
    sj_trans_out=cat(1,sj_trans_out,tempo_trans);
end
end
sj_trans_out(end+1:end+365) = raw(:,4); % 2019 songjung discharge
sj_trans_out(end+1:end+365 + leapyear(2020)) = [raw1(1:59,4); raw1(59,4); raw1(60:end,4)]; % 2020 songjumg discharge has 365 but 2020 is leap year 
% sj_trans_out(end+1:end+365) = raw(:,4); % 2020 songjung discharge

clearvars tempo_temp riv_temp riv_temp_pre
order_i = 0
for i = t_year(1):t_year(end) % 1980~2018
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_temp=gawha_re_w_c{t_year(order_i)-1979};
if order_i == 1
    riv_temp=tempo_temp;
else
    riv_temp=cat(1,riv_temp,tempo_temp);
end
end
riv_temp(end+1:end+365) = gawha_re_w_c{40}; % 2019
temp_temp=gawha_re_w_c{41};
% temp_temp(60)=[];
riv_temp(end+1:end+365 + leapyear(2020)) = temp_temp; % 2020
% riv_temp form it to 1989 ~ 2020
% riv_temp = NaN(1,length(sj_trans_out));
% riv_temp(1,365+leapyear(1989)+1:length(sj_trans_out)) = riv_temp_pre;

% riv_temp(end+1:end+365) = temp_temp; % 2020
 
 clearvars tempo_temp riv_temp_365 riv_temp_365_pre
order_i = 0
for i = t_year(1):t_year(end) % 1990~2018
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_temp=gawha_re_w_c{t_year(order_i)-1979};
if order_i == 1
    riv_temp_365=tempo_temp;
else
     if length(tempo_temp) == 366
        tempo_temp(60)=[];
    end
    riv_temp_365=cat(1,riv_temp_365,tempo_temp);
end
end
riv_temp_365(end+1:end+365) = gawha_re_w_c{40}; % 2019
temp_temp=gawha_re_w_c{41};
temp_temp(60)=[];
riv_temp_365(end+1:end+365) = temp_temp; % 2020
 % riv_temp form it to 1989 ~ 2020
% riv_temp_365 = NaN(1,length(sj_trans_out_365));
% riv_temp_365(1,365+leapyear(1989)+1:length(sj_trans_out_365)) = riv_temp_365_pre;

%% monthly mean (yy.mm) form after over 3sigma value extraction.
for i=1:length(t_indx) % time for num. of months (until 2019.12.31)
if i ==1
   monthly_trans(i) = nanmean(sj_trans_out(1:t_indx(i)));
else
   monthly_trans(i) = nanmean(sj_trans_out(t_indx(i-1)+1:t_indx(i)));
end
end

riv_temp_t_ind=t_indx; % 1990~2020
k=0;
clearvars monthly_riv_temp
for i=13:length(t_indx) % time for num. of months (until 2019.12.31)
    k=k+1
   monthly_riv_temp(k) = nanmean(riv_temp(riv_temp_t_ind(i-1)+1:riv_temp_t_ind(i)));
end


%% 1982 ~ 1986 : 1st
%% 1987 ~ 1996 : 2nd
%% 1997 ~ 2006 : 3rd
%% 2007 ~ 2015 : 4th
%% 2016 ~ 2020 : 5th

% t_day_cut1=365*length(1982:1986);
% t_day_cut2=365*length(1982:1996);
% t_day_cut3=365*length(1982:2006);
% t_day_cut4=365*length(1982:2015);
% 
% sj_day_trans_1=sj_trans_out_365(1:t_day_cut1);
% sj_day_trans_2=sj_trans_out_365(t_day_cut1+1:t_day_cut2);
% sj_day_trans_3=sj_trans_out_365(t_day_cut2+1:t_day_cut3);
% sj_day_trans_4=sj_trans_out_365(t_day_cut3+1:t_day_cut4);
% sj_day_trans_5=sj_trans_out_365(t_day_cut4+1:end);
% 
% sj_day_riv_temp_1=riv_temp_365(1:t_day_cut1);
% sj_day_riv_temp_2=riv_temp_365(t_day_cut1+1:t_day_cut2);
% sj_day_riv_temp_3=riv_temp_365(t_day_cut2+1:t_day_cut3);
% sj_day_riv_temp_4=riv_temp_365(t_day_cut3+1:t_day_cut4);
% sj_day_riv_temp_5=riv_temp_365(t_day_cut4+1:end);


%% 1982 ~ 1986 : 1st
%% 1987 ~ 1992 : 2nd
%% 1993 ~ 2006 : 3rd
%% 2007 ~ 2015 : 4th
%% 2016 ~ 2020 : 5th

t_day_cut1=365*length(1982:1986);
t_day_cut2=365*length(1982:1992);
t_day_cut3=365*length(1982:2004);
t_day_cut4=365*length(1982:2006);
t_day_cut5=365*length(1982:2015);

sj_day_trans_1=sj_trans_out_365(1:t_day_cut1);
sj_day_trans_2=sj_trans_out_365(t_day_cut1+1:t_day_cut2);
sj_day_trans_3=sj_trans_out_365(t_day_cut2+1:t_day_cut3);
sj_day_trans_4=sj_trans_out_365(t_day_cut4+1:t_day_cut5);
sj_day_trans_5=sj_trans_out_365(t_day_cut5+1:end);

sj_day_riv_temp_1=riv_temp_365(1:t_day_cut1);
sj_day_riv_temp_2=riv_temp_365(t_day_cut1+1:t_day_cut2);
sj_day_riv_temp_3=riv_temp_365(t_day_cut2+1:t_day_cut3);
sj_day_riv_temp_4=riv_temp_365(t_day_cut4+1:t_day_cut5);
sj_day_riv_temp_5=riv_temp_365(t_day_cut5+1:end);

% sj_day_riv_temp_1=riv_temp_365(t_day_cut2-365+1:t_day_cut1-365);
% sj_day_riv_temp_2=riv_temp_365(t_day_cut1-365+1:t_day_cut2-365);
% sj_day_riv_temp_3=riv_temp_365(t_day_cut2-365+1:end);

for i = 1:365 %days  
clim_day_trans_1(i) = nanmean(sj_day_trans_1(i:365:end));
clim_day_trans_2(i) = nanmean(sj_day_trans_2(i:365:end));
clim_day_trans_3(i) = nanmean(sj_day_trans_3(i:365:end));
clim_day_trans_4(i) = nanmean(sj_day_trans_4(i:365:end));
clim_day_trans_5(i) = nanmean(sj_day_trans_5(i:365:end));
clim_day_riv_temp_1(i) = nanmean(sj_day_riv_temp_1(i:365:end));
clim_day_riv_temp_2(i) = nanmean(sj_day_riv_temp_2(i:365:end));
clim_day_riv_temp_3(i) = nanmean(sj_day_riv_temp_3(i:365:end));
clim_day_riv_temp_4(i) = nanmean(sj_day_riv_temp_4(i:365:end));
clim_day_riv_temp_5(i) = nanmean(sj_day_riv_temp_5(i:365:end)); 
end

t_s = 8*12; %cut to 1996
monthly_trans(1:t_s) = [];
% monthly_do(1:t_s) = [];
% monthly_chl(1:t_s) = [];
% monthly_nh4(1:t_s) = [];
% monthly_no3(1:t_s) = [];
% monthly_po4(1:t_s) = [];
% monthly_ss(1:t_s) = [];

%% cut 3 regime
% cut_1_dis = 120 % 2006-12
% cut_2_dis = 228 % 2015-12

cut_1 = 120 % 2006-12
cut_2 = 228 % 2015-12

mon_trans_1 = monthly_trans(1:cut_1);
mon_trans_2 = monthly_trans(cut_1+1:cut_2);
mon_trans_3 = monthly_trans(cut_2+1:end);

mon_do_1 = monthly_do(1:cut_1);

mon_chl_1 = monthly_chl(1:cut_1);

mon_nh4_1 = monthly_nh4(1:cut_1);

mon_no3_1 = monthly_no3(1:cut_1);

mon_po4_1 = monthly_po4(1:cut_1);

mon_ss_1 = monthly_ss(1:cut_1);


for i = 1:12 %month   
mon_clim_trans_1(i) = nanmean(mon_trans_1(i:12:end));
mon_clim_trans_2(i) = nanmean(mon_trans_2(i:12:end));
mon_clim_trans_3(i) = nanmean(mon_trans_3(i:12:end)); 

mon_clim_do_1(i) = nanmean(mon_do_1(i:12:end));

mon_clim_chl_1(i) = nanmean(mon_chl_1(i:12:end));

mon_clim_nh4_1(i) = nanmean(mon_nh4_1(i:12:end));

mon_clim_no3_1(i) = nanmean(mon_no3_1(i:12:end));

mon_clim_po4_1(i) = nanmean(mon_po4_1(i:12:end));

mon_clim_ss_1(i) = nanmean(mon_ss_1(i:12:end));
end

figure; hold on;
plot(clim_day_trans_1,'r');  grid on;

figure; hold on;
plot(mon_clim_do_1,'r');  grid on;

figure; hold on;
plot(mon_clim_chl_1,'r');  grid on;

figure; hold on;
plot(mon_clim_nh4_1,'r');  grid on;

figure; hold on;
plot(mon_clim_no3_1,'r');  grid on;

figure; hold on;
plot(mon_clim_po4_1,'r'); grid on;

figure; hold on;
plot(mon_clim_ss_1,'r'); grid on;


t_indx(1:12) -15
day_indx=[1; 45; 75; 105; 136; 166; 197; 228; 258; 289; 319; 365;];

clearvars day_clim_*
day_clim_do_1 = NaN(365,1);

day_clim_chl_1 = NaN(365,1);

day_clim_nh4_1 = NaN(365,1);

day_clim_no3_1 = NaN(365,1);

day_clim_po4_1 = NaN(365,1);

day_clim_ss_1 = NaN(365,1);


%%
day_clim_do_1(day_indx) = mon_clim_do_1;

day_clim_chl_1(day_indx) = mon_clim_chl_1;

day_clim_nh4_1(day_indx) = mon_clim_nh4_1;

day_clim_no3_1(day_indx) = mon_clim_no3_1;

day_clim_po4_1(day_indx) = mon_clim_po4_1;

day_clim_ss_1(day_indx) = mon_clim_ss_1;

%%
t=1:length(day_clim_do_1);
day_clim_do_1(isnan(day_clim_do_1)) = interp1(t(~isnan(day_clim_do_1)),day_clim_do_1(~isnan(day_clim_do_1)),t(isnan(day_clim_do_1)));

t=1:length(day_clim_chl_1);
day_clim_chl_1(isnan(day_clim_chl_1)) = interp1(t(~isnan(day_clim_chl_1)),day_clim_chl_1(~isnan(day_clim_chl_1)),t(isnan(day_clim_chl_1)));

t=1:length(day_clim_nh4_1);
day_clim_nh4_1(isnan(day_clim_nh4_1)) = interp1(t(~isnan(day_clim_nh4_1)),day_clim_nh4_1(~isnan(day_clim_nh4_1)),t(isnan(day_clim_nh4_1)));

t=1:length(day_clim_no3_1);
day_clim_no3_1(isnan(day_clim_no3_1))  = interp1(t(~isnan(day_clim_no3_1)),day_clim_no3_1(~isnan(day_clim_no3_1)),t(isnan(day_clim_no3_1)));

t=1:length(day_clim_po4_1);
day_clim_po4_1(isnan(day_clim_po4_1)) = interp1(t(~isnan(day_clim_po4_1)),day_clim_po4_1(~isnan(day_clim_po4_1)),t(isnan(day_clim_po4_1)));

t=1:length(day_clim_ss_1);
day_clim_ss_1(isnan(day_clim_ss_1)) = interp1(t(~isnan(day_clim_ss_1)),day_clim_ss_1(~isnan(day_clim_ss_1)),t(isnan(day_clim_ss_1)));

save('namgang_climate_days_bio_1982to2020_high_airT_mixed_daily_v2.mat','clim_day_*','day_clim_*');


figure; hold on;
plot(clim_day_trans_1,'r'); plot(clim_day_trans_2,'g'); plot(clim_day_trans_3,'b'); 
plot(clim_day_trans_4,'m'); plot(clim_day_trans_5,'c');
grid on;
text(2,650,num2str(mean(clim_day_trans_1),'%0.2f'),'color','r');
text(2,600,num2str(mean(clim_day_trans_2),'%0.2f'),'color','g');
text(2,550,num2str(mean(clim_day_trans_3),'%0.2f'),'color','b');
text(2,500,num2str(mean(clim_day_trans_4),'%0.2f'),'color','b');
text(2,450,num2str(mean(clim_day_trans_5),'%0.2f'),'color','b');
ylabel('discharge(m^3/s)'); xlabel('days'); xlim([1 365])
ylim([0 1000])

total_d=eomday(2021,[1:12])
%
winter_t=[1:sum(total_d(1:2)),sum(total_d(1:11))+1:sum(total_d(1:12))];
jan_t=[1:31];
feb_t=[total_d(1)+1:sum(total_d(1:2))];
dec_t=[sum(total_d(1:11))+1:sum(total_d(1:12))];
%
summer_t=[sum(total_d(1:5))+1:sum(total_d(1:8))];
june_t=[sum(total_d(1:5))+1:sum(total_d(1:6))];
july_t=[sum(total_d(1:6))+1:sum(total_d(1:7))];
aug_t=[sum(total_d(1:7))+1:sum(total_d(1:8))];
aug_t_half=[sum(total_d(1:7))+15:sum(total_d(1:8))];

feb_cum_dis=[sum(clim_day_trans_1(feb_t) .* 86400), sum(clim_day_trans_2(feb_t) .* 86400), ...
    sum(clim_day_trans_3(feb_t) .* 86400),sum(clim_day_trans_4(feb_t) .* 86400),sum(clim_day_trans_5(feb_t) .* 86400)]

aug_cum_dis=[sum(clim_day_trans_1(aug_t) .* 86400), sum(clim_day_trans_2(aug_t) .* 86400), ...
    sum(clim_day_trans_3(aug_t) .* 86400),sum(clim_day_trans_4(aug_t) .* 86400),sum(clim_day_trans_5(aug_t) .* 86400)]

aug_cum_dis_half=[sum(clim_day_trans_1(aug_t_half) .* 86400), sum(clim_day_trans_2(aug_t_half) .* 86400), ...
    sum(clim_day_trans_3(aug_t_half) .* 86400),sum(clim_day_trans_4(aug_t_half) .* 86400),sum(clim_day_trans_5(aug_t_half) .* 86400)]


winter_cum_dis=[sum(clim_day_trans_1(winter_t) .* 86400), sum(clim_day_trans_2(winter_t) .* 86400), ...
    sum(clim_day_trans_3(winter_t) .* 86400),sum(clim_day_trans_4(winter_t) .* 86400),sum(clim_day_trans_5(winter_t) .* 86400)]

summer_cum_dis=[sum(clim_day_trans_1(summer_t) .* 86400), sum(clim_day_trans_2(summer_t) .* 86400), ...
    sum(clim_day_trans_3(summer_t) .* 86400),sum(clim_day_trans_4(summer_t) .* 86400),sum(clim_day_trans_5(summer_t) .* 86400)]

jan_cum_dis=[sum(clim_day_trans_1(jan_t) .* 86400), sum(clim_day_trans_2(jan_t) .* 86400), ...
    sum(clim_day_trans_3(jan_t) .* 86400),sum(clim_day_trans_4(jan_t) .* 86400),sum(clim_day_trans_5(jan_t) .* 86400)]

july_cum_dis=[sum(clim_day_trans_1(july_t) .* 86400), sum(clim_day_trans_2(july_t) .* 86400), ...
    sum(clim_day_trans_3(july_t) .* 86400),sum(clim_day_trans_4(july_t) .* 86400),sum(clim_day_trans_5(july_t) .* 86400)]

dec_cum_dis=[sum(clim_day_trans_1(dec_t) .* 86400), sum(clim_day_trans_2(dec_t) .* 86400), ...
    sum(clim_day_trans_3(dec_t) .* 86400),sum(clim_day_trans_4(dec_t) .* 86400),sum(clim_day_trans_5(dec_t) .* 86400)]

june_cum_dis=[sum(clim_day_trans_1(june_t) .* 86400), sum(clim_day_trans_2(june_t) .* 86400), ...
    sum(clim_day_trans_3(june_t) .* 86400),sum(clim_day_trans_4(june_t) .* 86400),sum(clim_day_trans_5(june_t) .* 86400)]


plot(sum([aug_cum_dis_half; july_cum_dis;]))


figure; hold on;
plot(feb_cum_dis./1000000,'linew',2);
ylabel('2월 누적유량(Mt/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');


figure; hold on;
plot(aug_cum_dis./1000000,'linew',2);
ylabel('8월 누적유량(Mt/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');


figure; hold on;
plot(winter_cum_dis./1000000,'linew',2);
ylabel('겨울(12~2월) 누적유량(Mt/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');


figure; hold on;
plot(summer_cum_dis./1000000,'linew',2);
ylabel('여름(6~8월) 누적유량(Mt/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');

%%
figure; hold on;
% plot(dec_cum_dis./1000000,'r','linew',2); % 1,000,000
% plot(jan_cum_dis./1000000,'b','linew',2);
% plot(feb_cum_dis./1000000,'k','linew',2);
plot(dec_cum_dis,'r','linew',2); % 1,000,000
plot(jan_cum_dis,'b','linew',2);
plot(feb_cum_dis,'k','linew',2);
% ylabel('겨울철 누적유량(Mt/월)');
ylabel('월별 평균 방류량(톤/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');
legend('12월','1월','2월');

figure; hold on;
% plot(june_cum_dis./1000000,'r','linew',2);
% plot(july_cum_dis./1000000,'b','linew',2);
% plot(aug_cum_dis./1000000,'k','linew',2);
plot(june_cum_dis,'r','linew',2);
plot(july_cum_dis,'b','linew',2);
plot(aug_cum_dis,'k','linew',2);
% ylabel('여름철 누적유량(Mt/월)');
ylabel('월별 평균 방류량(톤/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');
legend('6월','7월','8월');

close all; clear all; clc;

pathff = 'D:\장기생태\Dynamic\06_river\data\GY_input\';

r_list1=dir([pathff,'river_*to*_reg_bio_*']);
r_list2=dir([pathff,'river_*th_reg_bio_*']);

r_list = [r_list1; r_list2]

for i = 1:length(r_list)
    r_tran(:,:,i)=ncread([pathff,r_list(i).name],'river_transport');
    r_nh4(:,:,:,i)=ncread([pathff,r_list(i).name],'river_NH4');
    r_no3(:,:,:,i)=ncread([pathff,r_list(i).name],'river_NO3');
end

r_s_tran = squeeze(r_tran(1,:,:));
r_s_nh4 = squeeze(r_nh4(1,20,:,:));
r_s_no3 = squeeze(r_no3(1,20,:,:));


% load nut. 
clearvars plt_nh4 plt_no3 plt_po4
plt_nh4=r_s_nh4 ./1000 .*14.006720.* r_s_tran .* 86400 ./ 10^6 
plt_no3=r_s_no3 ./1000 .*14.006720.* r_s_tran .* 86400 ./ 10^6 
plt_din = plt_nh4 + plt_no3;
% plt_po4=sj.yr_po4 ./1000 .*30.973762.* r_s_tran .* 86400 ./ 10^6
fig = figure; hold on;
for i = 1:5
%         plot((i-1)*365+1:365*i,plt_no3(:,i),'ro','linew',2);
%         plot((i-1)*365+1:365*i,plt_nh4(:,i),'ko','linew',2);
         plot((i-1)*365+1:365*i,plt_din(:,i),'b','linew',2);
         plot((i-1)*365+1:365*i,plt_no3(:,i),'r','linew',2);
        plot((i-1)*365+1:365*i,plt_nh4(:,i),'k','linew',2);
end


fig = figure; hold on;
for i = 1:5
%         plot((i-1)*365+1:365*i,plt_no3(:,i),'ro','linew',2);
%         plot((i-1)*365+1:365*i,plt_nh4(:,i),'ko','linew',2);
         plot((i-1)*365+1:365*i,plt_din(:,i),'b','linew',2);
end

clim_day_trans_1 = plt_din(:,1);
clim_day_trans_2 = plt_din(:,2);
clim_day_trans_3 = plt_din(:,3);
clim_day_trans_4 = plt_din(:,4);
clim_day_trans_5 = plt_din(:,5);


total_d=eomday(2021,[1:12])
%
winter_t=[1:sum(total_d(1:2)),sum(total_d(1:11))+1:sum(total_d(1:12))];
winter_t_1=[1:sum(total_d(1)),sum(total_d(1:11))+1:sum(total_d(1:12))];
winter_t_2=[sum(total_d(1:10))+1:sum(total_d(1:12))];

jan_t=[1:31];
feb_t=[total_d(1)+1:sum(total_d(1:2))];
dec_t=[sum(total_d(1:11))+1:sum(total_d(1:12))];
%
summer_t=[sum(total_d(1:5))+1:sum(total_d(1:8))];
june_t=[sum(total_d(1:5))+1:sum(total_d(1:6))];
july_t=[sum(total_d(1:6))+1:sum(total_d(1:7))];
aug_t=[sum(total_d(1:7))+1:sum(total_d(1:8))];
aug_t_half=[sum(total_d(1:7))+15:sum(total_d(1:8))];
aug_t_5=[sum(total_d(1:7))+26:sum(total_d(1:8))];


feb_cum_dis=[sum(clim_day_trans_1(feb_t) ), sum(clim_day_trans_2(feb_t) ), ...
    sum(clim_day_trans_3(feb_t) ),sum(clim_day_trans_4(feb_t) ),sum(clim_day_trans_5(feb_t) )]

aug_cum_dis=[sum(clim_day_trans_1(aug_t) ), sum(clim_day_trans_2(aug_t) ), ...
    sum(clim_day_trans_3(aug_t) ),sum(clim_day_trans_4(aug_t) ),sum(clim_day_trans_5(aug_t) )]

aug_cum_dis_half=[sum(clim_day_trans_1(aug_t_half) ), sum(clim_day_trans_2(aug_t_half) ), ...
    sum(clim_day_trans_3(aug_t_half) ),sum(clim_day_trans_4(aug_t_half) ),sum(clim_day_trans_5(aug_t_half) )]

aug_cum_dis_5=[sum(clim_day_trans_1(aug_t_5) ), sum(clim_day_trans_2(aug_t_5) ), ...
    sum(clim_day_trans_3(aug_t_5) ),sum(clim_day_trans_4(aug_t_5) ),sum(clim_day_trans_5(aug_t_5) )]


winter_cum_dis=[sum(clim_day_trans_1(winter_t) ), sum(clim_day_trans_2(winter_t) ), ...
    sum(clim_day_trans_3(winter_t) ),sum(clim_day_trans_4(winter_t) ),sum(clim_day_trans_5(winter_t) )]

winter_cum_dis_1=[sum(clim_day_trans_1(winter_t_1) ), sum(clim_day_trans_2(winter_t_1) ), ...
    sum(clim_day_trans_3(winter_t_1) ),sum(clim_day_trans_4(winter_t_1) ),sum(clim_day_trans_5(winter_t_1) )]

winter_cum_dis_2=[sum(clim_day_trans_1(winter_t_2) ), sum(clim_day_trans_2(winter_t_2) ), ...
    sum(clim_day_trans_3(winter_t_2) ),sum(clim_day_trans_4(winter_t_2) ),sum(clim_day_trans_5(winter_t_2) )]


summer_cum_dis=[sum(clim_day_trans_1(summer_t) ), sum(clim_day_trans_2(summer_t) ), ...
    sum(clim_day_trans_3(summer_t) ),sum(clim_day_trans_4(summer_t) ),sum(clim_day_trans_5(summer_t) )]

jan_cum_dis=[sum(clim_day_trans_1(jan_t) ), sum(clim_day_trans_2(jan_t) ), ...
    sum(clim_day_trans_3(jan_t) ),sum(clim_day_trans_4(jan_t) ),sum(clim_day_trans_5(jan_t) )]

july_cum_dis=[sum(clim_day_trans_1(july_t) ), sum(clim_day_trans_2(july_t) ), ...
    sum(clim_day_trans_3(july_t) ),sum(clim_day_trans_4(july_t) ),sum(clim_day_trans_5(july_t) )]

dec_cum_dis=[sum(clim_day_trans_1(dec_t) ), sum(clim_day_trans_2(dec_t) ), ...
    sum(clim_day_trans_3(dec_t) ),sum(clim_day_trans_4(dec_t) ),sum(clim_day_trans_5(dec_t) )]

june_cum_dis=[sum(clim_day_trans_1(june_t) ), sum(clim_day_trans_2(june_t) ), ...
    sum(clim_day_trans_3(june_t) ),sum(clim_day_trans_4(june_t) ),sum(clim_day_trans_5(june_t) )]


figure; hold on;
plot(sum([july_cum_dis;aug_cum_dis;],1),'linew',2);
ylabel('7 ~ 8월 누적 유입량(ton/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');


figure; hold on;
plot(feb_cum_dis,'linew',2);
ylabel('2월 누적 DIN유입량(ton/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');


figure; hold on;
plot(aug_cum_dis,'linew',2);
ylabel('8월 누적 DIN유입량(ton/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');

figure; hold on;
plot(([aug_cum_dis_5;]),'linew',2);
ylabel('7~8월 누적 DIN유입량(ton/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');



figure; hold on;
plot(winter_cum_dis,'linew',2);
ylabel('겨울(12~2월) DIN유입량(ton/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');


figure; hold on;
plot(winter_cum_dis_1,'linew',2);
ylabel('겨울(12~1월) DIN유입량(ton/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');

figure; hold on;
plot(winter_cum_dis_2,'linew',2);
ylabel('겨울(11~12월) DIN유입량(ton/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');


figure; hold on;
plot(summer_cum_dis,'linew',2);
ylabel('여름(6~8월) DIN유입량(ton/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');

%%
figure; hold on;
plot(dec_cum_dis,'r','linew',2);
plot(jan_cum_dis,'b','linew',2);
plot(feb_cum_dis,'k','linew',2);
% ylabel('겨울철 누적 DIN유입량(톤/월)');
ylabel('월별 평균 용존무기질소 유입량(톤/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');
legend('12월','1월','2월');

figure; hold on;
plot(june_cum_dis,'r','linew',2);
plot(july_cum_dis,'b','linew',2);
plot(aug_cum_dis,'k','linew',2);
% ylabel('여름철 누적 DIN유입량(톤/월)');
ylabel('월별 평균 용존무기질소 유입량(톤/월)');
xticks(1:5); xticklabels({'82~86','87~92','93~04','07~15','16~20'});
grid on; box on; set(gca,'fontsize',15,'fontweight','bold');
legend('6월','7월','8월');


fig = figure; hold on;
         plot(sum(plt_din,1),'b','linew',2);      




fig = figure; hold on;
        plot(plt_no3,'ro','linew',2);
        plot(plt_nh4,'ko','linew',2);
         plot(plt_no3,'r','linew',2);
        plot(plt_nh4,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_no3(1:18)),18,1),':','color','r','linew',2);
        plot(19:27,repmat(nanmean(plt_no3(19:27)),length(19:27),1),':','color','r','linew',2);
        plot(28:31,repmat(nanmean(plt_no3(28:31)),length(28:31),1),':','color','r','linew',2);
        plot(1:18,repmat(nanmean(plt_nh4(1:18)),18,1),':','color','k','linew',2);
        plot(19:27,repmat(nanmean(plt_nh4(19:27)),length(19:27),1),':','color','k','linew',2);
        plot(28:31,repmat(nanmean(plt_nh4(28:31)),length(28:31),1),':','color','k','linew',2);
        title(['NO3-N, NH4-N, PO4-P load songjung']);
        xlabel('time','fontsize',13)
        ylabel('load (ton/yr)','fontsize',13)
        legend('NO3-N','NH4-N');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([7 31])
        ylim([0 inf])
%         ylim([0 5000])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_po4,'bo','linew',2);
plot(plt_po4,'b','linew',2);
plot(13:18,repmat(nanmean(plt_po4(13:18)),length(13:18),1),':','color','b','linew',2);
plot(19:27,repmat(nanmean(plt_po4(19:27)),length(19:27),1),':','color','b','linew',2);
plot(28:31,repmat(nanmean(plt_po4(28:31)),length(28:31),1),':','color','b','linew',2);
ylabel('PO4-P load (ton/yr)','fontsize',13)
ylim([0 inf])
xtickangle(45)
print(fig,strcat(['load_nut_river_yearly']),'-dpng')
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)


