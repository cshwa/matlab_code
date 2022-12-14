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
t_year = 1989:2018
% _________________________________________________________________________
% cd D:\장기생태\Dynamic\06_river
% load jinwall_raw_data_daily.mat %there has to be song vs. song_transp
%jinwall vs. jinwall_transp
% load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','raw_*')
cd D:\장기생태\Dynamic\06_river\환경과학원
% load songjung_yymm_koem_data_89to19_3sig.mat  % polynomial_fitting_on_ecological_variable_to2004_new_v6_3sigma.m
load songjung_yymm_monthly_data_to06to15_3sig_2020.mat  % polynomial_fitting_on_ecological_variable_to06to15_new_v7_3sigma_2020.m
% addP=load('songjung_bio_v3_day_koem_yr.mat'); % 2004_regime_v3.m

cd D:\장기생태\Dynamic\06_river
yj=load('koem_timeseires_monthly_gy_only_1to3sig_v2_each_st_16p.mat'); % 16 points 

cd D:\장기생태\Dynamic\06_river
clearvars merg_recon_w_c
% sjtemp=load('D:\장기생태\Dynamic\06_river_physics\sumjin_recons_water_temp_present_2019.mat');
% sjtemp=load('D:\장기생태\Dynamic\06_river_physics\환경과학원\sumjin_recons_water_temp(agyang)_present_2019_high_airT.mat');
sjtemp=load('D:\장기생태\Dynamic\06_river_physics\환경과학원\sumjin_recons_water_temp(agyang)_present_2020_high_airT_mixed_daily.mat');
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
order_i = 1
for i = t_year(2):t_year(end) % 1990~2018
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_temp=sumjin_re_w_c{t_year(order_i)-1989};
if order_i == 2
    riv_temp_pre=tempo_temp;
else
    riv_temp_pre=cat(1,riv_temp_pre,tempo_temp);
end
end
riv_temp_pre(end+1:end+365) = sumjin_re_w_c{30}; % 2019
temp_temp=sumjin_re_w_c{31};
% temp_temp(60)=[];
riv_temp_pre(end+1:end+365 + leapyear(2020)) = temp_temp; % 2020
% riv_temp form it to 1989 ~ 2020
riv_temp = NaN(1,length(sj_trans_out));
riv_temp(1,365+leapyear(1989)+1:length(sj_trans_out)) = riv_temp_pre;

% riv_temp(end+1:end+365) = temp_temp; % 2020
 
 clearvars tempo_temp riv_temp_365 riv_temp_365_pre
order_i = 1
for i = t_year(2):t_year(end) % 1990~2018
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_temp=sumjin_re_w_c{t_year(order_i)-1989};
if order_i == 2
    riv_temp_365_pre=tempo_temp;
else
     if length(tempo_temp) == 366
        tempo_temp(60)=[];
    end
    riv_temp_365_pre=cat(1,riv_temp_365_pre,tempo_temp);
end
end
riv_temp_365_pre(end+1:end+365) = sumjin_re_w_c{30}; % 2019
temp_temp=sumjin_re_w_c{31};
temp_temp(60)=[];
riv_temp_365_pre(end+1:end+365) = temp_temp; % 2020
 % riv_temp form it to 1989 ~ 2020
riv_temp_365 = NaN(1,length(sj_trans_out_365));
riv_temp_365(1,365+leapyear(1989)+1:length(sj_trans_out_365)) = riv_temp_365_pre;
 

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

t_day_cut0=365*length(1989:1996);
t_day_cut1=365*length(1989:2006);
t_day_cut2=365*length(1989:2015);
sj_day_trans_1=sj_trans_out_365(t_day_cut0+1:t_day_cut1);
sj_day_trans_2=sj_trans_out_365(t_day_cut1+1:t_day_cut2);
sj_day_trans_3=sj_trans_out_365(t_day_cut2+1:end);
sj_day_riv_temp_1=riv_temp_365(t_day_cut0+1:t_day_cut1);
sj_day_riv_temp_2=riv_temp_365(t_day_cut1+1:t_day_cut2);
sj_day_riv_temp_3=riv_temp_365(t_day_cut2+1:end);

% sj_day_riv_temp_1=riv_temp_365(t_day_cut0-365+1:t_day_cut1-365);
% sj_day_riv_temp_2=riv_temp_365(t_day_cut1-365+1:t_day_cut2-365);
% sj_day_riv_temp_3=riv_temp_365(t_day_cut2-365+1:end);

for i = 1:365 %days  
clim_day_trans_1(i) = nanmean(sj_day_trans_1(i:365:end));
clim_day_trans_2(i) = nanmean(sj_day_trans_2(i:365:end));
clim_day_trans_3(i) = nanmean(sj_day_trans_3(i:365:end));
clim_day_riv_temp_1(i) = nanmean(sj_day_riv_temp_1(i:365:end));
clim_day_riv_temp_2(i) = nanmean(sj_day_riv_temp_2(i:365:end));
clim_day_riv_temp_3(i) = nanmean(sj_day_riv_temp_3(i:365:end)); 
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
cut_1 = 120 % 2006-12
cut_2 = 228 % 2015-12

mon_trans_1 = monthly_trans(1:cut_1);
mon_trans_2 = monthly_trans(cut_1+1:cut_2);
mon_trans_3 = monthly_trans(cut_2+1:end);

mon_do_1 = monthly_do(1:cut_1);
mon_do_2 = monthly_do(cut_1+1:cut_2);
mon_do_3 = monthly_do(cut_2+1:end);

mon_chl_1 = monthly_chl(1:cut_1);
mon_chl_2 = monthly_chl(cut_1+1:cut_2);
mon_chl_3 = monthly_chl(cut_2+1:end);

mon_nh4_1 = monthly_nh4(1:cut_1);
mon_nh4_2 = monthly_nh4(cut_1+1:cut_2);
mon_nh4_3 = monthly_nh4(cut_2+1:end);

mon_no3_1 = monthly_no3(1:cut_1);
mon_no3_2 = monthly_no3(cut_1+1:cut_2);
mon_no3_3 = monthly_no3(cut_2+1:end);

mon_po4_1 = monthly_po4(1:cut_1);
mon_po4_2 = monthly_po4(cut_1+1:cut_2);
mon_po4_3 = monthly_po4(cut_2+1:end);

mon_ss_1 = monthly_ss(1:cut_1);
mon_ss_2 = monthly_ss(cut_1+1:cut_2);
mon_ss_3 = monthly_ss(cut_2+1:end);

mon_tn_1 = monthly_tn(1:cut_1);
mon_tn_2 = monthly_tn(cut_1+1:cut_2);
mon_tn_3 = monthly_tn(cut_2+1:end);

mon_tp_1 = monthly_tp(1:cut_1);
mon_tp_2 = monthly_tp(cut_1+1:cut_2);
mon_tp_3 = monthly_tp(cut_2+1:end);


for i = 1:12 %month   
mon_clim_trans_1(i) = nanmean(mon_trans_1(i:12:end));
mon_clim_trans_2(i) = nanmean(mon_trans_2(i:12:end));
mon_clim_trans_3(i) = nanmean(mon_trans_3(i:12:end)); 

mon_clim_do_1(i) = nanmean(mon_do_1(i:12:end));
mon_clim_do_2(i) = nanmean(mon_do_2(i:12:end));
mon_clim_do_3(i) = nanmean(mon_do_3(i:12:end));
    
mon_clim_chl_1(i) = nanmean(mon_chl_1(i:12:end));
mon_clim_chl_2(i) = nanmean(mon_chl_2(i:12:end));
mon_clim_chl_3(i) = nanmean(mon_chl_3(i:12:end));

mon_clim_nh4_1(i) = nanmean(mon_nh4_1(i:12:end));
mon_clim_nh4_2(i) = nanmean(mon_nh4_2(i:12:end));
mon_clim_nh4_3(i) = nanmean(mon_nh4_3(i:12:end));

mon_clim_no3_1(i) = nanmean(mon_no3_1(i:12:end));
mon_clim_no3_2(i) = nanmean(mon_no3_2(i:12:end));
mon_clim_no3_3(i) = nanmean(mon_no3_3(i:12:end));

mon_clim_po4_1(i) = nanmean(mon_po4_1(i:12:end));
mon_clim_po4_2(i) = nanmean(mon_po4_2(i:12:end));
mon_clim_po4_3(i) = nanmean(mon_po4_3(i:12:end));

mon_clim_ss_1(i) = nanmean(mon_ss_1(i:12:end));
mon_clim_ss_2(i) = nanmean(mon_ss_2(i:12:end));
mon_clim_ss_3(i) = nanmean(mon_ss_3(i:12:end));

mon_clim_tn_1(i) = nanmean(mon_tn_1(i:12:end));
mon_clim_tn_2(i) = nanmean(mon_tn_2(i:12:end));
mon_clim_tn_3(i) = nanmean(mon_tn_3(i:12:end));

mon_clim_tp_1(i) = nanmean(mon_tp_1(i:12:end));
mon_clim_tp_2(i) = nanmean(mon_tp_2(i:12:end));
mon_clim_tp_3(i) = nanmean(mon_tp_3(i:12:end));
end

figure; hold on;
plot(clim_day_trans_1,'r'); plot(clim_day_trans_2,'g'); plot(clim_day_trans_3,'b'); grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_do_1,'r'); plot(mon_clim_do_2,'g'); plot(mon_clim_do_3,'b'); grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_chl_1,'r'); plot(mon_clim_chl_2,'g'); plot(mon_clim_chl_3,'b'); grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_nh4_1,'r'); plot(mon_clim_nh4_2,'g'); plot(mon_clim_nh4_3,'b'); grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_no3_1,'r'); plot(mon_clim_no3_2,'g'); plot(mon_clim_no3_3,'b'); grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_po4_1,'r'); plot(mon_clim_po4_2,'g'); plot(mon_clim_po4_3,'b'); grid on;
legend('1997~2006','2007~2015','2016~2020');

figure; hold on;
plot(mon_clim_ss_1,'r'); plot(mon_clim_ss_2,'g'); plot(mon_clim_ss_3,'b'); grid on;
legend('1997~2006','2007~2015','2016~2020');


t_indx(1:12) -15
day_indx=[1; 45; 75; 105; 136; 166; 197; 228; 258; 289; 319; 365;];

clearvars day_clim_*
day_clim_do_1 = NaN(365,1);
day_clim_do_2 = NaN(365,1);
day_clim_do_3 = NaN(365,1);

day_clim_chl_1 = NaN(365,1);
day_clim_chl_2 = NaN(365,1);
day_clim_chl_3 = NaN(365,1);

day_clim_nh4_1 = NaN(365,1);
day_clim_nh4_2 = NaN(365,1);
day_clim_nh4_3 = NaN(365,1);

day_clim_no3_1 = NaN(365,1);
day_clim_no3_2 = NaN(365,1);
day_clim_no3_3 = NaN(365,1);

day_clim_po4_1 = NaN(365,1);
day_clim_po4_2 = NaN(365,1);
day_clim_po4_3 = NaN(365,1);

day_clim_ss_1 = NaN(365,1);
day_clim_ss_2 = NaN(365,1);
day_clim_ss_3 = NaN(365,1);

day_clim_tn_1 = NaN(365,1);
day_clim_tn_2 = NaN(365,1);
day_clim_tn_3 = NaN(365,1);

day_clim_tp_1 = NaN(365,1);
day_clim_tp_2 = NaN(365,1);
day_clim_tp_3 = NaN(365,1);

%%
day_clim_do_1(day_indx) = mon_clim_do_1;
day_clim_do_2(day_indx) = mon_clim_do_2;
day_clim_do_3(day_indx) = mon_clim_do_3;

day_clim_chl_1(day_indx) = mon_clim_chl_1;
day_clim_chl_2(day_indx) = mon_clim_chl_2;
day_clim_chl_3(day_indx) = mon_clim_chl_3;

day_clim_nh4_1(day_indx) = mon_clim_nh4_1;
day_clim_nh4_2(day_indx) = mon_clim_nh4_2;
day_clim_nh4_3(day_indx) = mon_clim_nh4_3;

day_clim_no3_1(day_indx) = mon_clim_no3_1;
day_clim_no3_2(day_indx) = mon_clim_no3_2;
day_clim_no3_3(day_indx) = mon_clim_no3_3;

day_clim_po4_1(day_indx) = mon_clim_po4_1;
day_clim_po4_2(day_indx) = mon_clim_po4_2;
day_clim_po4_3(day_indx) = mon_clim_po4_3;

day_clim_ss_1(day_indx) = mon_clim_ss_1;
day_clim_ss_2(day_indx) = mon_clim_ss_2;
day_clim_ss_3(day_indx) = mon_clim_ss_3;

day_clim_tn_1(day_indx) = mon_clim_tn_1;
day_clim_tn_2(day_indx) = mon_clim_tn_2;
day_clim_tn_3(day_indx) = mon_clim_tn_3;

day_clim_tp_1(day_indx) = mon_clim_tp_1;
day_clim_tp_2(day_indx) = mon_clim_tp_2;
day_clim_tp_3(day_indx) = mon_clim_tp_3;

%%
t=1:length(day_clim_do_1);
day_clim_do_1(isnan(day_clim_do_1)) = interp1(t(~isnan(day_clim_do_1)),day_clim_do_1(~isnan(day_clim_do_1)),t(isnan(day_clim_do_1)));
t=1:length(day_clim_do_2);
day_clim_do_2(isnan(day_clim_do_2)) = interp1(t(~isnan(day_clim_do_2)),day_clim_do_2(~isnan(day_clim_do_2)),t(isnan(day_clim_do_2)));
t=1:length(day_clim_do_3);
day_clim_do_3(isnan(day_clim_do_3)) = interp1(t(~isnan(day_clim_do_3)),day_clim_do_3(~isnan(day_clim_do_3)),t(isnan(day_clim_do_3)));

t=1:length(day_clim_chl_1);
day_clim_chl_1(isnan(day_clim_chl_1)) = interp1(t(~isnan(day_clim_chl_1)),day_clim_chl_1(~isnan(day_clim_chl_1)),t(isnan(day_clim_chl_1)));
t=1:length(day_clim_chl_2);
day_clim_chl_2(isnan(day_clim_chl_2)) = interp1(t(~isnan(day_clim_chl_2)),day_clim_chl_2(~isnan(day_clim_chl_2)),t(isnan(day_clim_chl_2)));
t=1:length(day_clim_chl_3);
day_clim_chl_3(isnan(day_clim_chl_3)) = interp1(t(~isnan(day_clim_chl_3)),day_clim_chl_3(~isnan(day_clim_chl_3)),t(isnan(day_clim_chl_3)));

t=1:length(day_clim_nh4_1);
day_clim_nh4_1(isnan(day_clim_nh4_1)) = interp1(t(~isnan(day_clim_nh4_1)),day_clim_nh4_1(~isnan(day_clim_nh4_1)),t(isnan(day_clim_nh4_1)));
t=1:length(day_clim_nh4_2);
day_clim_nh4_2(isnan(day_clim_nh4_2)) = interp1(t(~isnan(day_clim_nh4_2)),day_clim_nh4_2(~isnan(day_clim_nh4_2)),t(isnan(day_clim_nh4_2)));
t=1:length(day_clim_nh4_3);
day_clim_nh4_3(isnan(day_clim_nh4_3)) = interp1(t(~isnan(day_clim_nh4_3)),day_clim_nh4_3(~isnan(day_clim_nh4_3)),t(isnan(day_clim_nh4_3)));

t=1:length(day_clim_no3_1);
day_clim_no3_1(isnan(day_clim_no3_1))  = interp1(t(~isnan(day_clim_no3_1)),day_clim_no3_1(~isnan(day_clim_no3_1)),t(isnan(day_clim_no3_1)));
t=1:length(day_clim_no3_2);
day_clim_no3_2(isnan(day_clim_no3_2))  = interp1(t(~isnan(day_clim_no3_2)),day_clim_no3_2(~isnan(day_clim_no3_2)),t(isnan(day_clim_no3_2)));
t=1:length(day_clim_no3_3);
day_clim_no3_3(isnan(day_clim_no3_3))  = interp1(t(~isnan(day_clim_no3_3)),day_clim_no3_3(~isnan(day_clim_no3_3)),t(isnan(day_clim_no3_3)));

t=1:length(day_clim_po4_1);
day_clim_po4_1(isnan(day_clim_po4_1)) = interp1(t(~isnan(day_clim_po4_1)),day_clim_po4_1(~isnan(day_clim_po4_1)),t(isnan(day_clim_po4_1)));
t=1:length(day_clim_po4_2);
day_clim_po4_2(isnan(day_clim_po4_2)) = interp1(t(~isnan(day_clim_po4_2)),day_clim_po4_2(~isnan(day_clim_po4_2)),t(isnan(day_clim_po4_2)));
t=1:length(day_clim_po4_3);
day_clim_po4_3(isnan(day_clim_po4_3)) = interp1(t(~isnan(day_clim_po4_3)),day_clim_po4_3(~isnan(day_clim_po4_3)),t(isnan(day_clim_po4_3)));

t=1:length(day_clim_ss_1);
day_clim_ss_1(isnan(day_clim_ss_1)) = interp1(t(~isnan(day_clim_ss_1)),day_clim_ss_1(~isnan(day_clim_ss_1)),t(isnan(day_clim_ss_1)));
t=1:length(day_clim_ss_2);
day_clim_ss_2(isnan(day_clim_ss_2)) = interp1(t(~isnan(day_clim_ss_2)),day_clim_ss_2(~isnan(day_clim_ss_2)),t(isnan(day_clim_ss_2)));
t=1:length(day_clim_ss_3);
day_clim_ss_3(isnan(day_clim_ss_3)) = interp1(t(~isnan(day_clim_ss_3)),day_clim_ss_3(~isnan(day_clim_ss_3)),t(isnan(day_clim_ss_3)));

t=1:length(day_clim_tn_1);
day_clim_tn_1(isnan(day_clim_tn_1)) = interp1(t(~isnan(day_clim_tn_1)),day_clim_tn_1(~isnan(day_clim_tn_1)),t(isnan(day_clim_tn_1)));
t=1:length(day_clim_tn_2);
day_clim_tn_2(isnan(day_clim_tn_2)) = interp1(t(~isnan(day_clim_tn_2)),day_clim_tn_2(~isnan(day_clim_tn_2)),t(isnan(day_clim_tn_2)));
t=1:length(day_clim_tn_3);
day_clim_tn_3(isnan(day_clim_tn_3)) = interp1(t(~isnan(day_clim_tn_3)),day_clim_tn_3(~isnan(day_clim_tn_3)),t(isnan(day_clim_tn_3)));

t=1:length(day_clim_tp_1);
day_clim_tp_1(isnan(day_clim_tp_1)) = interp1(t(~isnan(day_clim_tp_1)),day_clim_tp_1(~isnan(day_clim_tp_1)),t(isnan(day_clim_tp_1)));
t=1:length(day_clim_tp_2);
day_clim_tp_2(isnan(day_clim_tp_2)) = interp1(t(~isnan(day_clim_tp_2)),day_clim_tp_2(~isnan(day_clim_tp_2)),t(isnan(day_clim_tp_2)));
t=1:length(day_clim_tp_3);
day_clim_tp_3(isnan(day_clim_tp_3)) = interp1(t(~isnan(day_clim_tp_3)),day_clim_tp_3(~isnan(day_clim_tp_3)),t(isnan(day_clim_tp_3)));

plt_filter_DJF = [1:31, 32:59, 335:365];

figure; hold on;
plot(clim_day_trans_2(plt_filter_DJF),'g'); plot(clim_day_trans_3(plt_filter_DJF),'b'); grid on;
yline(mean(clim_day_trans_2(plt_filter_DJF)),'g--','linew',2)
yline(mean(clim_day_trans_3(plt_filter_DJF)),'b--','linew',2)
text(5,45,num2str(mean(clim_day_trans_2(plt_filter_DJF)),'%0.2f'),'color','g');
text(5,40,num2str(mean(clim_day_trans_3(plt_filter_DJF)),'%0.2f'),'color','b');
ylabel('discharge(m^3/s)'); xlabel('time(days)'); xlim([1 length(plt_filter_DJF)]); ylim([0 50])
xticks([32, 60]); xticklabels({})
hA = gca;
%Tick properties can be set in the X/Y/ZRuler.MinorTicks
hA.XAxis.MinorTick='on';
%Tick locations can be set in the X/Y/ZRuler.MinorTick
hA.XAxis.MinorTickValues = [32, 60]
text(13,2,'Jan','fontsize',15);
text(43,2,'Feb','fontsize',15);
text(71,2,'Dec','fontsize',15);

figure; hold on;
plot(clim_day_trans_2,'g'); plot(clim_day_trans_3,'b'); grid on;
text(2,600,num2str(mean(clim_day_trans_2),'%0.2f'),'color','g');
text(2,550,num2str(mean(clim_day_trans_3),'%0.2f'),'color','b');


figure; hold on;
plot(clim_day_trans_1,'r'); plot(clim_day_trans_2,'g'); plot(clim_day_trans_3,'b'); grid on;
text(2,650,num2str(mean(clim_day_trans_1),'%0.2f'),'color','r');
text(2,600,num2str(mean(clim_day_trans_2),'%0.2f'),'color','g');
text(2,550,num2str(mean(clim_day_trans_3),'%0.2f'),'color','b');
ylabel('discharge(m^3/s)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(clim_day_riv_temp_1,'r'); plot(clim_day_riv_temp_2,'g'); plot(clim_day_riv_temp_3,'b'); grid on;
text(2,28,num2str(mean(clim_day_riv_temp_1),'%0.2f'),'color','r');
text(2,25,num2str(mean(clim_day_riv_temp_2),'%0.2f'),'color','g');
text(2,22,num2str(mean(clim_day_riv_temp_3),'%0.2f'),'color','b');
ylabel('temp(^oC)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_chl_1,'r'); plot(day_clim_chl_2,'g'); plot(day_clim_chl_3,'b'); grid on;
text(2,13.5,num2str(mean(day_clim_chl_1),'%0.2f'),'color','r');
text(2,12.5,num2str(mean(day_clim_chl_2),'%0.2f'),'color','g');
text(2,11.5,num2str(mean(day_clim_chl_3),'%0.2f'),'color','b');
ylabel('Chl.a(ug/L)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_do_1,'r'); plot(day_clim_do_2,'g'); plot(day_clim_do_3,'b'); grid on;
text(2,300,num2str(mean(day_clim_do_1),'%0.2f'),'color','r');
text(2,280,num2str(mean(day_clim_do_2),'%0.2f'),'color','g');
text(2,260,num2str(mean(day_clim_do_3),'%0.2f'),'color','b');
ylabel('DO(mmol/m^3)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_nh4_1,'r'); plot(day_clim_nh4_2,'g'); plot(day_clim_nh4_3,'b'); grid on;
text(2,13.5,num2str(mean(day_clim_nh4_1),'%0.2f'),'color','r');
text(2,12.5,num2str(mean(day_clim_nh4_2),'%0.2f'),'color','g');
text(2,11.5,num2str(mean(day_clim_nh4_3),'%0.2f'),'color','b');
ylabel('NH4-N (mmol N /m^3)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_no3_1,'r'); plot(day_clim_no3_2,'g'); plot(day_clim_no3_3,'b'); grid on;
text(2,80,num2str(mean(day_clim_no3_1),'%0.2f'),'color','r');
text(2,75,num2str(mean(day_clim_no3_2),'%0.2f'),'color','g');
text(2,70,num2str(mean(day_clim_no3_3),'%0.2f'),'color','b');
ylabel('NO3-N (mmol N /m^3)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_po4_1,'r'); plot(day_clim_po4_2,'g'); plot(day_clim_po4_3,'b'); grid on;
text(2,1.7,num2str(mean(day_clim_po4_1),'%0.2f'),'color','r');
text(2,1.5,num2str(mean(day_clim_po4_2),'%0.2f'),'color','g');
text(2,1.3,num2str(mean(day_clim_po4_3),'%0.2f'),'color','b');
% plot((day_clim_po4_2 + day_clim_po4_3)./2 ,'m');
% text(2,1.1,num2str(mean((day_clim_po4_2 + day_clim_po4_3)./2),'%0.2f'),'color','m');
ylabel('PO4-P (mmol P /m^3)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_ss_1,'r'); plot(day_clim_ss_2,'g'); plot(day_clim_ss_3,'b'); grid on;
text(2,11,num2str(mean(day_clim_ss_1),'%0.2f'),'color','r');
text(2,10,num2str(mean(day_clim_ss_2),'%0.2f'),'color','g');
text(2,9,num2str(mean(day_clim_ss_3),'%0.2f'),'color','b');
ylabel('SS (mg/L)'); xlabel('days'); xlim([1 365])

% load climate
figure; hold on;
plot(clim_day_trans_1,'r'); plot(clim_day_trans_2,'g'); 
plot(clim_day_trans_3,'b'); grid on;
text(2,200,num2str(mean(clim_day_trans_1),'%0.2f'),'color','r');
text(2,150,num2str(mean(clim_day_trans_2),'%0.2f'),'color','g');
text(2,100,num2str(mean(clim_day_trans_3),'%0.2f'),'color','b');
yline(mean(clim_day_trans_1),'r--'); yline(mean(clim_day_trans_2),'g--'); 
yline(mean(clim_day_trans_3),'b--'); 
ylabel('discharge(m^3/s)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(clim_day_trans_1.*(day_clim_po4_1'./1000 .*30.973762),'r'); plot(clim_day_trans_2.*(day_clim_po4_2'./1000 .*30.973762),'g'); 
plot(clim_day_trans_3.*(day_clim_po4_3'./1000 .*30.973762),'b'); grid on;
text(2,18,num2str(mean(clim_day_trans_1.*(day_clim_po4_1'./1000 .*30.973762)),'%0.2f'),'color','r');
text(2,16,num2str(mean(clim_day_trans_2.*(day_clim_po4_2'./1000 .*30.973762)),'%0.2f'),'color','g');
text(2,14,num2str(mean(clim_day_trans_3.*(day_clim_po4_3'./1000 .*30.973762)),'%0.2f'),'color','b');
ylabel('load P (g/s)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(clim_day_trans_1.*(day_clim_nh4_1'./1000 .*14.006720),'r'); plot(clim_day_trans_2.*(day_clim_nh4_2'./1000 .*14.006720),'g'); 
plot(clim_day_trans_3.*(day_clim_nh4_3'./1000 .*14.006720),'b'); grid on;
text(2,80,num2str(mean(clim_day_trans_1.*(day_clim_nh4_1'./1000 .*14.006720)),'%0.2f'),'color','r');
text(2,70,num2str(mean(clim_day_trans_2.*(day_clim_nh4_2'./1000 .*14.006720)),'%0.2f'),'color','g');
text(2,60,num2str(mean(clim_day_trans_3.*(day_clim_nh4_3'./1000 .*14.006720)),'%0.2f'),'color','b');
ylabel('load NH4-N (g/s)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(clim_day_trans_1.*(day_clim_no3_1'./1000 .*14.006720),'r'); plot(clim_day_trans_2.*(day_clim_no3_2'./1000 .*14.006720),'g'); 
plot(clim_day_trans_3.*(day_clim_no3_3'./1000 .*14.006720),'b'); grid on;
text(2,800,num2str(mean(clim_day_trans_1.*(day_clim_no3_1'./1000 .*14.006720)),'%0.2f'),'color','r');
text(2,700,num2str(mean(clim_day_trans_2.*(day_clim_no3_2'./1000 .*14.006720)),'%0.2f'),'color','g');
text(2,600,num2str(mean(clim_day_trans_3.*(day_clim_no3_3'./1000 .*14.006720)),'%0.2f'),'color','b');
ylabel('load NO3-N (g/s)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(clim_day_trans_1.*(day_clim_ss_1'),'r'); plot(clim_day_trans_2.*(day_clim_ss_2'),'g'); 
plot(clim_day_trans_3.*(day_clim_ss_3'),'b'); grid on;
text(2,5000,num2str(mean(clim_day_trans_1.*(day_clim_ss_1')),'%0.2f'),'color','r');
text(2,4000,num2str(mean(clim_day_trans_2.*(day_clim_ss_2')),'%0.2f'),'color','g');
text(2,3000,num2str(mean(clim_day_trans_3.*(day_clim_ss_3')),'%0.2f'),'color','b');
ylabel('load SS (g/s)'); xlabel('days'); xlim([1 365])

figure; hold on; 
plot((monthly_po4./1000 .*14.006720).* monthly_trans(1:end))


% sj_tn_load_1 = clim_day_trans_1.*(day_clim_tn_1'); sj_tn_load_2 = clim_day_trans_2.*(day_clim_tn_2');
% sj_tn_load_3 = clim_day_trans_3.*(day_clim_tn_3');
% sj_tp_load_1 = clim_day_trans_1.*(day_clim_tp_1'); sj_tp_load_2 = clim_day_trans_2.*(day_clim_tp_2');
% sj_tp_load_3 = clim_day_trans_3.*(day_clim_tp_3');
% save('songjung_climate_days_load_to07to16.mat','sj_*_load_*');

figure; hold on;
plot(clim_day_trans_1.*(day_clim_tn_1'),'r'); plot(clim_day_trans_2.*(day_clim_tn_2'),'g'); 
plot(clim_day_trans_3.*(day_clim_tn_3'),'b'); grid on;
text(2,5000,num2str(mean(clim_day_trans_1.*(day_clim_tn_1')),'%0.2f'),'color','r');
text(2,4000,num2str(mean(clim_day_trans_2.*(day_clim_tn_2')),'%0.2f'),'color','g');
text(2,3000,num2str(mean(clim_day_trans_3.*(day_clim_tn_3')),'%0.2f'),'color','b');
ylabel('load TN (g/s)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(clim_day_trans_1.*(day_clim_tp_1'),'r'); plot(clim_day_trans_2.*(day_clim_tp_2'),'g'); 
plot(clim_day_trans_3.*(day_clim_tp_3'),'b'); grid on;
text(2,5000,num2str(mean(clim_day_trans_1.*(day_clim_tp_1')),'%0.2f'),'color','r');
text(2,4000,num2str(mean(clim_day_trans_2.*(day_clim_tp_2')),'%0.2f'),'color','g');
text(2,3000,num2str(mean(clim_day_trans_3.*(day_clim_tp_3')),'%0.2f'),'color','b');
ylabel('load TP (g/s)'); xlabel('days'); xlim([1 365])


clearvars yr_trans yr_po4
monthly_trans_2 = monthly_trans;
% monthly_trans_2(1:t_s) = [];
monthly_po4_2 = monthly_po4;
% monthly_po4_2(1:t_s) = [];
for i = 1:23 %year
yr_trans(i) = nanmean(monthly_trans_2((i-1)*12+1:12*i));
yr_po4(i) = nanmean(monthly_po4_2((i-1)*12+1:12*i));
end

figure; hold on; 
plot((yr_po4./1000 .*14.006720).* yr_trans)

figure; hold on; 
mon_load_po4=(monthly_po4_2./1000 .*14.006720).* monthly_trans_2;
plot(mon_load_po4,'bo')
plot(find(isnan(mon_load_po4)==0),mon_load_po4(~isnan(mon_load_po4)))
for i = 1:23
line((i-1)*12+1:12*i,repmat(((yr_po4(i)./1000 .*14.006720).* yr_trans(i)),12,1),'color','r')
end
line(1:12*10,repmat(nanmean((yr_po4(1:10)./1000 .*14.006720).* yr_trans(1:10)),length(1:12*10),1),'color','g')
line(12*10+1:19*12,repmat(nanmean((yr_po4(11:19)./1000 .*14.006720).* yr_trans(11:19)),length(12*10+1:19*12),1),'color','g')
line(12*19+1:23*12,repmat(nanmean((yr_po4(20:23)./1000 .*14.006720).* yr_trans(20:23)),length(12*19+1:23*12),1),'color','g')
ylabel('PO4-P (mmol P /m^3)'); xlabel('time'); xticks(1:12:length(monthly_po4_2)); xticklabels(1997:2020);
xtickangle(45); grid on; xlim([1 length(monthly_po4_2)])

figure; hold on;
for i = 1:10
plot(monthly_po4_2((i-1)*12+1:12*i),'ro')
end
for i = 11:19
plot(monthly_po4_2((i-1)*12+1:12*i),'go')
end
for i = 20:23
plot(monthly_po4_2((i-1)*12+1:12*i),'bo')
end

figure; hold on;
for i = 1:10
    plt_po4=monthly_po4_2((i-1)*12+1:12*i);
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'r');
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'ro');
end
for i = 11:19
    plt_po4=monthly_po4_2((i-1)*12+1:12*i);
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'g')
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'go')
end
for i = 20:23
    plt_po4=monthly_po4_2((i-1)*12+1:12*i);
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'b')
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'bo')
end
grid on; ylabel('PO4-P (mmol P /m^3)'); xlabel('time(month)'); xticks(1:12); xticklabels(1:12);
 xlim([1 12])

save('songjung_climate_days_bio_to07to16_high_airT_mixed_daily_2020.mat','clim_day_*','day_clim_*');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NAMGANG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
% -------------------------------------------------------------------------
t_year = 1989:2018
% _________________________________________________________________________
% cd D:\장기생태\Dynamic\06_river
% load jinwall_raw_data_daily.mat %there has to be song vs. song_transp
%jinwall vs. jinwall_transp
% load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','raw_*')
cd D:\장기생태\Dynamic\06_river\환경과학원
% load songjung_yymm_koem_data_89to19_3sig.mat  % polynomial_fitting_on_ecological_variable_to2004_new_v6_3sigma.m
load namgang_yymm_monthly_data_to06to15_3sig_2020.mat  % polynomial_fitting_on_ecological_variable_to06to15_new_v7_3sigma_2020.m
% addP=load('songjung_bio_v3_day_koem_yr.mat'); % 2004_regime_v3.m

cd D:\장기생태\Dynamic\06_river
yj=load('koem_timeseires_monthly_gy_only_1to3sig_v2_each_st_16p.mat'); % 16 points 

cd D:\장기생태\Dynamic\06_river
clearvars merg_recon_w_c
ghtemp=load('D:\장기생태\Dynamic\06_river_physics\gawha_recons_water_temp_present_2020_high_airT_mixed_daily.mat');
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
sj_trans_out_365(end+1:end+365) = raw(:,4);

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

clearvars tempo_temp riv_temp riv_temp_365
order_i = 1
for i = t_year(2):t_year(end) % 1990~2018
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_temp=gawha_re_w_c{t_year(order_i)-1989};
if order_i == 2
    riv_temp_pre=tempo_temp;
else
    riv_temp_pre=cat(1,riv_temp_pre,tempo_temp);
end
end
 riv_temp_pre(end+1:end+365) = gawha_re_w_c{30}; % 2019
 temp_temp=gawha_re_w_c{31};
% temp_temp(60)=[];
riv_temp_pre(end+1:end+365 + leapyear(2020)) = temp_temp; % 2020
% riv_temp form it to 1989 ~ 2020
riv_temp = NaN(1,11688); % 11688 means 1989:2020
riv_temp(1,365+leapyear(1989)+1:end) = riv_temp_pre;
% riv_temp(1,length(1989:1999)*365+sum(leapyear(1989:1999))+1:end) = riv_temp_pre;

 clearvars tempo_temp riv_temp_365 riv_temp_365_pre
order_i = 1
for i = t_year(2):t_year(end) % 1990~2018
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_temp=gawha_re_w_c{t_year(order_i)-1989};
if order_i == 2
    riv_temp_365_pre=tempo_temp;
else
     if length(tempo_temp) == 366
        tempo_temp(60)=[];
    end
    riv_temp_365_pre=cat(1,riv_temp_365_pre,tempo_temp);
end
end
 riv_temp_365_pre(end+1:end+365) = gawha_re_w_c{30}; % 2019
temp_temp=gawha_re_w_c{31};
temp_temp(60)=[];
riv_temp_365_pre(end+1:end+365) = temp_temp; % 2020
 % riv_temp form it to 1989 ~ 2020
riv_temp_365 = NaN(1,11688 - sum(leapyear(1989:2020)));
riv_temp_365(1,365+leapyear(1989)+1:end) = riv_temp_365_pre;

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


t_day_cut0=365*length(1989:1996);
t_day_cut1=365*length(1989:2006);
t_day_cut2=365*length(1989:2015);

sj_day_trans_1=sj_trans_out_365(t_day_cut0+1:t_day_cut1);
sj_day_trans_2=sj_trans_out_365(t_day_cut1+1:t_day_cut2);
sj_day_trans_3=sj_trans_out_365(t_day_cut2+1:end);
gw_day_riv_temp_1=riv_temp_365(t_day_cut0+1:t_day_cut1);
gw_day_riv_temp_2=riv_temp_365(t_day_cut1+1:t_day_cut2);
gw_day_riv_temp_3=riv_temp_365(t_day_cut2+1:end);

for i = 1:365 %days  
clim_day_trans_1(i) = nanmean(sj_day_trans_1(i:365:end));
clim_day_trans_2(i) = nanmean(sj_day_trans_2(i:365:end));
clim_day_trans_3(i) = nanmean(sj_day_trans_3(i:365:end));
clim_day_riv_temp_1(i) = nanmean(gw_day_riv_temp_1(i:365:end));
clim_day_riv_temp_2(i) = nanmean(gw_day_riv_temp_2(i:365:end));
clim_day_riv_temp_3(i) = nanmean(gw_day_riv_temp_3(i:365:end)); 
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
mon_do_2 = monthly_do(cut_1+1:cut_2);
mon_do_3 = monthly_do(cut_2+1:end);

mon_chl_1 = monthly_chl(1:cut_1);
mon_chl_2 = monthly_chl(cut_1+1:cut_2);
mon_chl_3 = monthly_chl(cut_2+1:end);

mon_nh4_1 = monthly_nh4(1:cut_1);
mon_nh4_2 = monthly_nh4(cut_1+1:cut_2);
mon_nh4_3 = monthly_nh4(cut_2+1:end);

mon_no3_1 = monthly_no3(1:cut_1);
mon_no3_2 = monthly_no3(cut_1+1:cut_2);
mon_no3_3 = monthly_no3(cut_2+1:end);

mon_po4_1 = monthly_po4(1:cut_1);
mon_po4_2 = monthly_po4(cut_1+1:cut_2);
mon_po4_3 = monthly_po4(cut_2+1:end);

mon_ss_1 = monthly_ss(1:cut_1);
mon_ss_2 = monthly_ss(cut_1+1:cut_2);
mon_ss_3 = monthly_ss(cut_2+1:end);


for i = 1:12 %month   
mon_clim_trans_1(i) = nanmean(mon_trans_1(i:12:end));
mon_clim_trans_2(i) = nanmean(mon_trans_2(i:12:end));
mon_clim_trans_3(i) = nanmean(mon_trans_3(i:12:end)); 

mon_clim_do_1(i) = nanmean(mon_do_1(i:12:end));
mon_clim_do_2(i) = nanmean(mon_do_2(i:12:end));
mon_clim_do_3(i) = nanmean(mon_do_3(i:12:end));
    
mon_clim_chl_1(i) = nanmean(mon_chl_1(i:12:end));
mon_clim_chl_2(i) = nanmean(mon_chl_2(i:12:end));
mon_clim_chl_3(i) = nanmean(mon_chl_3(i:12:end));

mon_clim_nh4_1(i) = nanmean(mon_nh4_1(i:12:end));
mon_clim_nh4_2(i) = nanmean(mon_nh4_2(i:12:end));
mon_clim_nh4_3(i) = nanmean(mon_nh4_3(i:12:end));

mon_clim_no3_1(i) = nanmean(mon_no3_1(i:12:end));
mon_clim_no3_2(i) = nanmean(mon_no3_2(i:12:end));
mon_clim_no3_3(i) = nanmean(mon_no3_3(i:12:end));

mon_clim_po4_1(i) = nanmean(mon_po4_1(i:12:end));
mon_clim_po4_2(i) = nanmean(mon_po4_2(i:12:end));
mon_clim_po4_3(i) = nanmean(mon_po4_3(i:12:end));

mon_clim_ss_1(i) = nanmean(mon_ss_1(i:12:end));
mon_clim_ss_2(i) = nanmean(mon_ss_2(i:12:end));
mon_clim_ss_3(i) = nanmean(mon_ss_3(i:12:end));
end

figure; hold on;
plot(clim_day_trans_1,'r'); plot(clim_day_trans_2,'g'); plot(clim_day_trans_3,'b'); grid on;

figure; hold on;
plot(mon_clim_do_1,'r'); plot(mon_clim_do_2,'g'); plot(mon_clim_do_3,'b'); grid on;

figure; hold on;
plot(mon_clim_chl_1,'r'); plot(mon_clim_chl_2,'g'); plot(mon_clim_chl_3,'b'); grid on;

figure; hold on;
plot(mon_clim_nh4_1,'r'); plot(mon_clim_nh4_2,'g'); plot(mon_clim_nh4_3,'b'); grid on;

figure; hold on;
plot(mon_clim_no3_1,'r'); plot(mon_clim_no3_2,'g'); plot(mon_clim_no3_3,'b'); grid on;

figure; hold on;
plot(mon_clim_po4_1,'r'); plot(mon_clim_po4_2,'g'); plot(mon_clim_po4_3,'b'); grid on;

figure; hold on;
plot(mon_clim_ss_1,'r'); plot(mon_clim_ss_2,'g'); plot(mon_clim_ss_3,'b'); grid on;


t_indx(1:12) -15
day_indx=[1; 45; 75; 105; 136; 166; 197; 228; 258; 289; 319; 365;];

clearvars day_clim_*
day_clim_do_1 = NaN(365,1);
day_clim_do_2 = NaN(365,1);
day_clim_do_3 = NaN(365,1);

day_clim_chl_1 = NaN(365,1);
day_clim_chl_2 = NaN(365,1);
day_clim_chl_3 = NaN(365,1);

day_clim_nh4_1 = NaN(365,1);
day_clim_nh4_2 = NaN(365,1);
day_clim_nh4_3 = NaN(365,1);

day_clim_no3_1 = NaN(365,1);
day_clim_no3_2 = NaN(365,1);
day_clim_no3_3 = NaN(365,1);

day_clim_po4_1 = NaN(365,1);
day_clim_po4_2 = NaN(365,1);
day_clim_po4_3 = NaN(365,1);

day_clim_ss_1 = NaN(365,1);
day_clim_ss_2 = NaN(365,1);
day_clim_ss_3 = NaN(365,1);

%%
day_clim_do_1(day_indx) = mon_clim_do_1;
day_clim_do_2(day_indx) = mon_clim_do_2;
day_clim_do_3(day_indx) = mon_clim_do_3;

day_clim_chl_1(day_indx) = mon_clim_chl_1;
day_clim_chl_2(day_indx) = mon_clim_chl_2;
day_clim_chl_3(day_indx) = mon_clim_chl_3;

day_clim_nh4_1(day_indx) = mon_clim_nh4_1;
day_clim_nh4_2(day_indx) = mon_clim_nh4_2;
day_clim_nh4_3(day_indx) = mon_clim_nh4_3;

day_clim_no3_1(day_indx) = mon_clim_no3_1;
day_clim_no3_2(day_indx) = mon_clim_no3_2;
day_clim_no3_3(day_indx) = mon_clim_no3_3;

day_clim_po4_1(day_indx) = mon_clim_po4_1;
day_clim_po4_2(day_indx) = mon_clim_po4_2;
day_clim_po4_3(day_indx) = mon_clim_po4_3;

day_clim_ss_1(day_indx) = mon_clim_ss_1;
day_clim_ss_2(day_indx) = mon_clim_ss_2;
day_clim_ss_3(day_indx) = mon_clim_ss_3;

%%
t=1:length(day_clim_do_1);
day_clim_do_1(isnan(day_clim_do_1)) = interp1(t(~isnan(day_clim_do_1)),day_clim_do_1(~isnan(day_clim_do_1)),t(isnan(day_clim_do_1)));
t=1:length(day_clim_do_2);
day_clim_do_2(isnan(day_clim_do_2)) = interp1(t(~isnan(day_clim_do_2)),day_clim_do_2(~isnan(day_clim_do_2)),t(isnan(day_clim_do_2)));
t=1:length(day_clim_do_3);
day_clim_do_3(isnan(day_clim_do_3)) = interp1(t(~isnan(day_clim_do_3)),day_clim_do_3(~isnan(day_clim_do_3)),t(isnan(day_clim_do_3)));

t=1:length(day_clim_chl_1);
day_clim_chl_1(isnan(day_clim_chl_1)) = interp1(t(~isnan(day_clim_chl_1)),day_clim_chl_1(~isnan(day_clim_chl_1)),t(isnan(day_clim_chl_1)));
t=1:length(day_clim_chl_2);
day_clim_chl_2(isnan(day_clim_chl_2)) = interp1(t(~isnan(day_clim_chl_2)),day_clim_chl_2(~isnan(day_clim_chl_2)),t(isnan(day_clim_chl_2)));
t=1:length(day_clim_chl_3);
day_clim_chl_3(isnan(day_clim_chl_3)) = interp1(t(~isnan(day_clim_chl_3)),day_clim_chl_3(~isnan(day_clim_chl_3)),t(isnan(day_clim_chl_3)));

t=1:length(day_clim_nh4_1);
day_clim_nh4_1(isnan(day_clim_nh4_1)) = interp1(t(~isnan(day_clim_nh4_1)),day_clim_nh4_1(~isnan(day_clim_nh4_1)),t(isnan(day_clim_nh4_1)));
t=1:length(day_clim_nh4_2);
day_clim_nh4_2(isnan(day_clim_nh4_2)) = interp1(t(~isnan(day_clim_nh4_2)),day_clim_nh4_2(~isnan(day_clim_nh4_2)),t(isnan(day_clim_nh4_2)));
t=1:length(day_clim_nh4_3);
day_clim_nh4_3(isnan(day_clim_nh4_3)) = interp1(t(~isnan(day_clim_nh4_3)),day_clim_nh4_3(~isnan(day_clim_nh4_3)),t(isnan(day_clim_nh4_3)));

t=1:length(day_clim_no3_1);
day_clim_no3_1(isnan(day_clim_no3_1))  = interp1(t(~isnan(day_clim_no3_1)),day_clim_no3_1(~isnan(day_clim_no3_1)),t(isnan(day_clim_no3_1)));
t=1:length(day_clim_no3_2);
day_clim_no3_2(isnan(day_clim_no3_2))  = interp1(t(~isnan(day_clim_no3_2)),day_clim_no3_2(~isnan(day_clim_no3_2)),t(isnan(day_clim_no3_2)));
t=1:length(day_clim_no3_3);
day_clim_no3_3(isnan(day_clim_no3_3))  = interp1(t(~isnan(day_clim_no3_3)),day_clim_no3_3(~isnan(day_clim_no3_3)),t(isnan(day_clim_no3_3)));

t=1:length(day_clim_po4_1);
day_clim_po4_1(isnan(day_clim_po4_1)) = interp1(t(~isnan(day_clim_po4_1)),day_clim_po4_1(~isnan(day_clim_po4_1)),t(isnan(day_clim_po4_1)));
t=1:length(day_clim_po4_2);
day_clim_po4_2(isnan(day_clim_po4_2)) = interp1(t(~isnan(day_clim_po4_2)),day_clim_po4_2(~isnan(day_clim_po4_2)),t(isnan(day_clim_po4_2)));
t=1:length(day_clim_po4_3);
day_clim_po4_3(isnan(day_clim_po4_3)) = interp1(t(~isnan(day_clim_po4_3)),day_clim_po4_3(~isnan(day_clim_po4_3)),t(isnan(day_clim_po4_3)));

t=1:length(day_clim_ss_1);
day_clim_ss_1(isnan(day_clim_ss_1)) = interp1(t(~isnan(day_clim_ss_1)),day_clim_ss_1(~isnan(day_clim_ss_1)),t(isnan(day_clim_ss_1)));
t=1:length(day_clim_ss_2);
day_clim_ss_2(isnan(day_clim_ss_2)) = interp1(t(~isnan(day_clim_ss_2)),day_clim_ss_2(~isnan(day_clim_ss_2)),t(isnan(day_clim_ss_2)));
t=1:length(day_clim_ss_3);
day_clim_ss_3(isnan(day_clim_ss_3)) = interp1(t(~isnan(day_clim_ss_3)),day_clim_ss_3(~isnan(day_clim_ss_3)),t(isnan(day_clim_ss_3)));


figure; hold on;
plot(clim_day_trans_1,'r'); plot(clim_day_trans_2,'g'); plot(clim_day_trans_3,'b'); grid on;
text(2,650,num2str(mean(clim_day_trans_1),'%0.2f'),'color','r');
text(2,600,num2str(mean(clim_day_trans_2),'%0.2f'),'color','g');
text(2,550,num2str(mean(clim_day_trans_3),'%0.2f'),'color','b');
ylabel('discharge(m^3/s)'); xlabel('days'); xlim([1 365])


figure; hold on;
plot(clim_day_riv_temp_1,'r'); plot(clim_day_riv_temp_2,'g'); plot(clim_day_riv_temp_3,'b'); grid on;
text(2,28,num2str(mean(clim_day_riv_temp_1),'%0.2f'),'color','r');
text(2,25,num2str(mean(clim_day_riv_temp_2),'%0.2f'),'color','g');
text(2,22,num2str(mean(clim_day_riv_temp_3),'%0.2f'),'color','b');
ylabel('temp(^oC)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_chl_1,'r'); plot(day_clim_chl_2,'g'); plot(day_clim_chl_3,'b'); grid on;
text(2,12.5,num2str(mean(day_clim_chl_1),'%0.2f'),'color','r');
text(2,11.5,num2str(mean(day_clim_chl_2),'%0.2f'),'color','g');
text(2,10.5,num2str(mean(day_clim_chl_3),'%0.2f'),'color','b');
ylabel('Chl.a(ug/L)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_do_1,'r'); plot(day_clim_do_2,'g'); plot(day_clim_do_3,'b'); grid on;
text(2,300,num2str(mean(day_clim_do_1),'%0.2f'),'color','r');
text(2,280,num2str(mean(day_clim_do_2),'%0.2f'),'color','g');
text(2,260,num2str(mean(day_clim_do_3),'%0.2f'),'color','b');
ylabel('DO(mmol/m^3)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_nh4_1,'r'); plot(day_clim_nh4_2,'g'); plot(day_clim_nh4_3,'b'); grid on;
text(2,13.5,num2str(mean(day_clim_nh4_1),'%0.2f'),'color','r');
text(2,12.5,num2str(mean(day_clim_nh4_2),'%0.2f'),'color','g');
text(2,11.5,num2str(mean(day_clim_nh4_3),'%0.2f'),'color','b');
ylabel('NH4-N (mmol N /m^3)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_no3_1,'r'); plot(day_clim_no3_2,'g'); plot(day_clim_no3_3,'b'); grid on;
text(2,80,num2str(mean(day_clim_no3_1),'%0.2f'),'color','r');
text(2,75,num2str(mean(day_clim_no3_2),'%0.2f'),'color','g');
text(2,70,num2str(mean(day_clim_no3_3),'%0.2f'),'color','b');
ylabel('NO3-N (mmol N /m^3)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_po4_1,'r'); plot(day_clim_po4_2,'g'); plot(day_clim_po4_3,'b'); grid on;
text(2,0.27,num2str(mean(day_clim_po4_1),'%0.2f'),'color','r');
text(2,0.25,num2str(mean(day_clim_po4_2),'%0.2f'),'color','g');
text(2,0.23,num2str(mean(day_clim_po4_3),'%0.2f'),'color','b');
% plot((day_clim_po4_2 + day_clim_po4_3)./2 ,'m');
% text(2,1.1,num2str(mean((day_clim_po4_2 + day_clim_po4_3)./2),'%0.2f'),'color','m');
ylabel('PO4-P (mmol p /m^3)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(day_clim_ss_1,'r'); plot(day_clim_ss_2,'g'); plot(day_clim_ss_3,'b'); grid on;
text(2,8,num2str(mean(day_clim_ss_1),'%0.2f'),'color','r');
text(2,7.5,num2str(mean(day_clim_ss_2),'%0.2f'),'color','g');
text(2,7,num2str(mean(day_clim_ss_3),'%0.2f'),'color','b');
ylabel('SS (mg/L)'); xlabel('days'); xlim([1 365])

% load climate
figure; hold on;
plot(clim_day_trans_1,'r'); plot(clim_day_trans_2,'g'); 
plot(clim_day_trans_3,'b'); grid on;
text(2,200,num2str(mean(clim_day_trans_1),'%0.2f'),'color','r');
text(2,150,num2str(mean(clim_day_trans_2),'%0.2f'),'color','g');
text(2,100,num2str(mean(clim_day_trans_3),'%0.2f'),'color','b');
yline(mean(clim_day_trans_1),'r--'); yline(mean(clim_day_trans_2),'g--'); 
yline(mean(clim_day_trans_3),'b--'); 
ylabel('discharge(m^3/s)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(clim_day_trans_1.*(day_clim_po4_1'./1000 .*30.973762),'r'); plot(clim_day_trans_2.*(day_clim_po4_2'./1000 .*30.973762),'g'); 
plot(clim_day_trans_3.*(day_clim_po4_3'./1000 .*30.973762),'b'); grid on;
text(2,4,num2str(mean(clim_day_trans_1.*(day_clim_po4_1'./1000 .*30.973762)),'%0.2f'),'color','r');
text(2,3.5,num2str(mean(clim_day_trans_2.*(day_clim_po4_2'./1000 .*30.973762)),'%0.2f'),'color','g');
text(2,3.0,num2str(mean(clim_day_trans_3.*(day_clim_po4_3'./1000 .*30.973762)),'%0.2f'),'color','b');
ylabel('load P (g/s)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(clim_day_trans_1.*(day_clim_nh4_1'./1000 .*14.006720),'r'); plot(clim_day_trans_2.*(day_clim_nh4_2'./1000 .*14.006720),'g'); 
plot(clim_day_trans_3.*(day_clim_nh4_3'./1000 .*14.006720),'b'); grid on;
text(2,55,num2str(mean(clim_day_trans_1.*(day_clim_nh4_1'./1000 .*14.006720)),'%0.2f'),'color','r');
text(2,50,num2str(mean(clim_day_trans_2.*(day_clim_nh4_2'./1000 .*14.006720)),'%0.2f'),'color','g');
text(2,45,num2str(mean(clim_day_trans_3.*(day_clim_nh4_3'./1000 .*14.006720)),'%0.2f'),'color','b');
ylabel('load NH4-N (g/s)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(clim_day_trans_1.*(day_clim_no3_1'./1000 .*14.006720),'r'); plot(clim_day_trans_2.*(day_clim_no3_2'./1000 .*14.006720),'g'); 
plot(clim_day_trans_3.*(day_clim_no3_3'./1000 .*14.006720),'b'); grid on;
text(2,650,num2str(mean(clim_day_trans_1.*(day_clim_no3_1'./1000 .*14.006720)),'%0.2f'),'color','r');
text(2,600,num2str(mean(clim_day_trans_2.*(day_clim_no3_2'./1000 .*14.006720)),'%0.2f'),'color','g');
text(2,550,num2str(mean(clim_day_trans_3.*(day_clim_no3_3'./1000 .*14.006720)),'%0.2f'),'color','b');
ylabel('load NO3-N (g/s)'); xlabel('days'); xlim([1 365])

figure; hold on;
plot(clim_day_trans_1.*(day_clim_ss_1'),'r'); plot(clim_day_trans_2.*(day_clim_ss_2'),'g'); 
plot(clim_day_trans_3.*(day_clim_ss_3'),'b'); grid on;
text(2,5000,num2str(mean(clim_day_trans_1.*(day_clim_ss_1')),'%0.2f'),'color','r');
text(2,4000,num2str(mean(clim_day_trans_2.*(day_clim_ss_2')),'%0.2f'),'color','g');
text(2,3000,num2str(mean(clim_day_trans_3.*(day_clim_ss_3')),'%0.2f'),'color','b');
ylabel('load SS (g/s)'); xlabel('days'); xlim([1 365])

figure; hold on; 
plot((monthly_po4./1000 .*14.006720).* monthly_trans(1:end))

clearvars yr_trans yr_po4
monthly_trans_2 = monthly_trans;
% monthly_trans_2(1:t_s) = [];
monthly_po4_2 = monthly_po4;
% monthly_po4_2(1:t_s) = [];
for i = 1:23 %year
yr_trans(i) = nanmean(monthly_trans_2((i-1)*12+1:12*i));
yr_po4(i) = nanmean(monthly_po4_2((i-1)*12+1:12*i));
end

figure; hold on; 
plot((yr_po4./1000 .*14.006720).* yr_trans)

figure; hold on; 
mon_load_po4=(monthly_po4_2./1000 .*14.006720).* monthly_trans_2;
plot(mon_load_po4,'bo')
plot(find(isnan(mon_load_po4)==0),mon_load_po4(~isnan(mon_load_po4)))
for i = 1:23
line((i-1)*12+1:12*i,repmat(((yr_po4(i)./1000 .*14.006720).* yr_trans(i)),12,1),'color','r')
end
line(1:12*10,repmat(nanmean((yr_po4(1:10)./1000 .*14.006720).* yr_trans(1:10)),length(1:12*10),1),'color','g')
line(12*10+1:19*12,repmat(nanmean((yr_po4(11:19)./1000 .*14.006720).* yr_trans(11:19)),length(12*10+1:19*12),1),'color','g')
line(12*19+1:23*12,repmat(nanmean((yr_po4(20:23)./1000 .*14.006720).* yr_trans(20:23)),length(12*19+1:23*12),1),'color','g')
ylabel('PO4-P (mmol P /m^3)'); xlabel('time'); xticks(1:12:length(monthly_po4_2)); xticklabels(1997:2020);
xtickangle(45); grid on; xlim([1 length(monthly_po4_2)])

figure; hold on;
for i = 1:10
plot(monthly_po4_2((i-1)*12+1:12*i),'ro')
end
for i = 11:19
plot(monthly_po4_2((i-1)*12+1:12*i),'go')
end
for i = 20:23
plot(monthly_po4_2((i-1)*12+1:12*i),'bo')
end

figure; hold on;
for i = 1:10
    plt_po4=monthly_po4_2((i-1)*12+1:12*i);
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'r');
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'ro');
end
for i = 11:19
    plt_po4=monthly_po4_2((i-1)*12+1:12*i);
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'g')
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'go')
end
for i = 20:23
    plt_po4=monthly_po4_2((i-1)*12+1:12*i);
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'b')
plot(find(isnan(plt_po4)==0),plt_po4(~isnan(plt_po4)),'bo')
end
grid on; ylabel('PO4-P (mmol P /m^3)'); xlabel('time(month)'); xticks(1:12); xticklabels(1:12);
 xlim([1 12])

save('namgang_climate_days_bio_to07to16_high_airT_mixed_daily_2020.mat','clim_day_*','day_clim_*');
