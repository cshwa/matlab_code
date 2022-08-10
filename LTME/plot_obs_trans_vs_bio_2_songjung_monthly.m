close all; clear; clc;
% -------------------------------------------------------------------------
t_year = 1989:2018
% _________________________________________________________________________
cd D:\장기생태\Dynamic\06_river
load jinwall_raw_data_daily.mat %there has to be song vs. song_transp
%jinwall vs. jinwall_transp
% load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','raw_*')
cd D:\장기생태\Dynamic\06_river\환경과학원
load songjung_yymm_monthly_data_89to19_3sig.mat  % polynomial_fitting_on_ecological_variable_to2004_new_v6_3sigma.m
% addP=load('songjung_bio_v3_day_mon_yr.mat'); % 2004_regime_v3.m

cd D:\장기생태\Dynamic\06_river\data\sj_1980to1996
%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
load songjung_discharge_1980to2018.mat  % pre_merg_dis is discharge
% load excel 2019 trans
[raw_sj txt_sj]=xlsread('D:\장기생태\Dynamic\06_river\환경과학원\수리·수문·기상_유량_구례_20190101-20191231.xls','유량');


cd D:\장기생태\Dynamic\06_river
% approximation : Sumjin(songjung) & Namgang have same discharge during 1980~1989.

load('sumjin_recons_water_temp_present.mat','merg_recon_w_c');
sumjin_re_w_c = merg_recon_w_c;

koem=load('plot_yr_trans_vs_koem_16points.mat'); %from plot_obs_trans_vs_bio_2_songjung_monthly_vs_koem.m


clearvars temperature sj_trans_out
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
sj_trans_out(end+1:end+365) = raw_sj(:,4);

%% monthly mean (yy.mm) form after over 3sigma value extraction.

for i=1:length(t_indx) % time for num. of months (until 2019.12.31)
if i ==1
   monthly_trans(i) = nanmean(sj_trans_out(1:t_indx(i)));
else
   monthly_trans(i) = nanmean(sj_trans_out(t_indx(i-1)+1:t_indx(i)));
end
end


figure; plot(monthly_trans)
hold on; for i = 6:12:360; xline(i,'color','r'); end

% fill missing value
mon_in_chl = monthly_chl;
mon_in_no3 = monthly_no3;
mon_in_nh4 = monthly_nh4;
mon_in_do = monthly_do;
mon_in_ss = monthly_ss;

t=1:length(t_indx);
mon_in_chl(isnan(mon_in_chl))=interp1(t(~isnan(mon_in_chl)),mon_in_chl(~isnan(mon_in_chl)),t(isnan(mon_in_chl)));
mon_in_no3(isnan(mon_in_no3))=interp1(t(~isnan(mon_in_no3)),mon_in_no3(~isnan(mon_in_no3)),t(isnan(mon_in_no3)));
mon_in_nh4(isnan(mon_in_nh4))=interp1(t(~isnan(mon_in_nh4)),mon_in_nh4(~isnan(mon_in_nh4)),t(isnan(mon_in_nh4)));
mon_in_do(isnan(mon_in_do))=interp1(t(~isnan(mon_in_do)),mon_in_do(~isnan(mon_in_do)),t(isnan(mon_in_do)));
mon_in_ss(isnan(mon_in_ss))=interp1(t(~isnan(mon_in_ss)),mon_in_ss(~isnan(mon_in_ss)),t(isnan(mon_in_ss)));


% raw_chl = raw_chl;
% raw_do = raw_do *0.7*44.661;
% raw_no3 = raw_no3 .*1000 ./14 ;
% raw_nh4 = raw_nh4 .*1000 ./14 ;
% monthly_chl(361:end) = [];
% monthly_no3(361:end) = [];
% monthly_nh4(361:end) = [];
% monthly_do(361:end) = [];
% monthly_po4(361:end) = [];
% monthly_ss(361:end) = [];
% 
% mon_in_chl(361:end) = [];
% mon_in_no3(361:end) = [];
% mon_in_nh4(361:end) = [];
% mon_in_do(361:end) = [];
% mon_in_po4(361:end) = [];
% mon_in_ss(361:end) = [];

june_chl = monthly_chl(6:12:end);
june_no3 = monthly_no3(6:12:end);
june_nh4 = monthly_nh4(6:12:end);
june_do = monthly_do(6:12:end);
june_po4 = monthly_po4(6:12:end);
june_ss = monthly_ss(6:12:end);

june_in_chl = mon_in_chl(6:12:end);
june_in_no3 = mon_in_no3(6:12:end);
june_in_nh4 = mon_in_nh4(6:12:end);
june_in_do = mon_in_do(6:12:end);
june_in_po4 = mon_in_po4(6:12:end);
june_in_ss = mon_in_ss(6:12:end);
june_trans = monthly_trans(6:12:end);

monthly_trans_x = monthly_trans - nanmean(monthly_trans);
mon_in_chl_x = mon_in_chl - nanmean(mon_in_chl);
mon_in_do_x = mon_in_do - nanmean(mon_in_do);
mon_in_no3_x = mon_in_no3 - nanmean(mon_in_no3);
mon_in_nh4_x = mon_in_nh4 - nanmean(mon_in_nh4);
mon_in_po4_x = mon_in_po4 - nanmean(mon_in_po4);
mon_in_ss_x = mon_in_ss - nanmean(mon_in_ss);

monthly_chl_x = monthly_chl - nanmean(monthly_chl);
monthly_do_x = monthly_do - nanmean(monthly_do);
monthly_no3_x = monthly_no3 - nanmean(monthly_no3);
monthly_nh4_x = monthly_nh4 - nanmean(monthly_nh4);
monthly_po4_x = monthly_po4 - nanmean(monthly_po4);
monthly_ss_x = monthly_ss - nanmean(monthly_ss);

june_trans_x = june_trans - nanmean(june_trans);
june_in_chl_x = june_in_chl - nanmean(june_in_chl);
june_in_do_x = june_in_do - nanmean(june_in_do);
june_in_no3_x = june_in_no3 - nanmean(june_in_no3);
june_in_nh4_x = june_in_nh4 - nanmean(june_in_nh4);
june_in_po4_x = june_in_po4 - nanmean(june_in_po4);
june_in_ss_x = june_in_ss - nanmean(june_in_ss);

for i = 1:31 %year
yr_trans(i) = nanmean(monthly_trans((i-1)*12+1:12*i));
yr_chl(i) = nanmean(monthly_chl((i-1)*12+1:12*i));
yr_nh4(i) = nanmean(monthly_nh4((i-1)*12+1:12*i));
yr_no3(i) = nanmean(monthly_no3((i-1)*12+1:12*i));
yr_po4(i) = nanmean(monthly_po4((i-1)*12+1:12*i));
yr_ss(i) = nanmean(monthly_ss((i-1)*12+1:12*i));
yr_in_chl(i) = nanmean(mon_in_chl((i-1)*12+1:12*i));
yr_in_nh4(i) = nanmean(mon_in_nh4((i-1)*12+1:12*i));
yr_in_no3(i) = nanmean(mon_in_no3((i-1)*12+1:12*i));
yr_in_po4(i) = nanmean(mon_in_po4((i-1)*12+1:12*i));
yr_in_ss(i) = nanmean(mon_in_ss((i-1)*12+1:12*i));
end

save('songjung_plt.mat','yr_*','monthly_*');

% cd D:\장기생태\Dynamic\06_river\환경과학원
% load 2007_2019_time_axis.mat
% t_tick(2:end) = t_tick(2:end)+1;
t_tick = 1:12:372 ; %1989.01.01~2020.01.01
%% YEARLY
%CHL vs. SS
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_ss,'b*','linew',2);
        plot(yr_ss,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS chl vs. SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
xtickangle(45)


%% SS
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_ss,'b*','linew',2);
        plot(yr_ss,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end

%% chl
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_chl,'b*','linew',2);
        plot(yr_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end

%% NH4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_nh4,'b*','linew',2);
        plot(yr_nh4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. nh4']);
        xlabel('time','fontsize',13)
        ylabel('NH4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_no3,'b*','linew',2);
        plot(yr_no3,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS transp vs. no3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_po4,'b*','linew',2);
        plot(yr_po4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS transp vs. PO4']);
        xlabel('time','fontsize',13)
        ylabel('PO4 (mmol P / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)


%% NH4 vs CHL
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_nh4 ./1000 .*14.006720.* yr_trans,'b*','linew',2);
        plot(yr_nh4 ./1000 .*14.006720.* yr_trans,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS transp vs. NH4']);
        xlabel('time','fontsize',13)
        ylabel('NH4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:30);
        xlim([1 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
xtickangle(45)

%% NO3 vs CHL
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_no3 ./1000 .*14.006720.* yr_trans,'b*','linew',2);
        plot(yr_no3 ./1000 .*14.006720.* yr_trans,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS transp vs. NO3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:30);
        xlim([1 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
xtickangle(45)

%% PO4 vs CHL
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_po4./1000 .*30.973762 .* yr_trans,'b*','linew',2);
        plot(yr_po4./1000 .*30.973762 .* yr_trans,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS transp vs. PO4']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:30);
        xlim([1 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
xtickangle(45)

%% ss vs CHL
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_ss .* yr_trans,'b*','linew',2);
        plot(yr_ss .* yr_trans,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS transp vs. SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L * m^3/s)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:30);
        xlim([1 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
xtickangle(45)

%% NH4/PO4 vs CHL
NP=((yr_no3 ./1000 .*14.006720) + (yr_nh4 ./1000 .*14.006720)) ./ (yr_po4./1000 .*30.973762)
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(NP .* yr_trans,'b*','linew',2);
        plot(NP .* yr_trans,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS transp vs. N:PO4']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:30);
        xlim([1 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(NP ,'b*','linew',2);
        plot(NP ,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS transp vs. NH4:PO4']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:30);
        xlim([1 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
xtickangle(45)


%% MONTHLY
%% SS
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(monthly_ss,'b*','linew',2);
        plot(mon_in_ss,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 monthly OBS transp vs. SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(monthly_chl)])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(monthly_trans,'g-','linew',2);
for i = 6:12:360; xline(i,'color','m'); end


%% chl
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(monthly_chl,'b*','linew',2);
        plot(mon_in_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 monthly OBS transp vs. chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(monthly_chl)])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(monthly_trans,'g-','linew',2);
for i = 6:12:360; xline(i,'color','m'); end

%% NH4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(monthly_nh4,'b*','linew',2);
        plot(mon_in_nh4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 monthly OBS transp vs. nh4']);
        xlabel('time','fontsize',13)
        ylabel('NH4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(monthly_chl)])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(monthly_trans,'g-','linew',2);
for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(monthly_no3,'b*','linew',2);
        plot(mon_in_no3,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 monthly OBS transp vs. no3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(monthly_chl)])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(monthly_trans,'g-','linew',2);
for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(monthly_po4,'b*','linew',2);
        plot(mon_in_po4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 monthly OBS transp vs. po4']);
        xlabel('time','fontsize',13)
        ylabel('po4 (mmol P / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(monthly_chl)])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(monthly_trans,'g-','linew',2);
for i = 6:12:360; xline(i,'color','m'); end

return

corrcoef(monthly_trans(~isnan(mon_in_chl)),mon_in_chl(~isnan(mon_in_chl)))
% corrcoef(sj_trans_out,interp_chl,'Rows','complete')
corrcoef(monthly_trans(~isnan(monthly_chl)),monthly_chl(~isnan(monthly_chl)))

figure;
plot(monthly_trans(~isnan(mon_in_chl)),mon_in_chl(~isnan(mon_in_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('Chla (ug/L)','fontsize',13)
grid on;
title(['GY 1989-2018 monthly interped chl vs. transport']);

figure;
plot(monthly_trans(~isnan(monthly_chl)),monthly_chl(~isnan(monthly_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('Chla (ug/L)','fontsize',13)
grid on;
title(['GY 1989-2018 monthly chl vs. transport']);


clearvars r lags
[r,lags] = xcorr(sj_trans_out_x(~isnan(interp_chl)),interp_chl_x(~isnan(interp_chl)),'normalized') 
% [r,lags] = xcorr(sj_trans_out_x(~isnan(raw_chl)),raw_chl_x(~isnan(raw_chl)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily chl vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL chl vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

sj_tran_in = sj_trans_out(~isnan(raw_chl));
raw_chl_in = raw_chl(~isnan(raw_chl));
clearvars r lags
[r,lags] = xcorr(sj_tran_in - nanmean(sj_tran_in),raw_chl_in - nanmean(raw_chl_in), 'normalized') 
% [r,lags] = xcorr(sj_trans_out(~isnan(interp_chl)),interp_chl(~isnan(interp_chl)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily chl vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');



%% nh4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(raw_nh4,'b*','linew',2);
        plot(interp_nh4,'r','linew',2);
%         plot(1:length(obm_gy_nh4),obm_gy_nh4 + obs_std_gy_nh4,'m-','linew',2);
%         plot(1:length(obm_gy_nh4),obm_gy_nh4 - obs_std_gy_nh4,'m-','linew',2);
        title(['GY 2001 daily OBS transp vs. nh4']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(raw_chl)])
        set(gca,'xticklabel',2007:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_nh4'),'-dpng') 
yyaxis right
plot(sj_trans_out,'g-','linew',2);
xlim([1 length(raw_nh4)])

return

corrcoef(sj_trans_out(~isnan(interp_nh4)),interp_nh4(~isnan(interp_nh4)))
% corrcoef(sj_trans_out,interp_nh4,'Rows','complete')
corrcoef(sj_trans_out(~isnan(raw_nh4)),raw_nh4(~isnan(raw_nh4)))

figure;
plot(sj_trans_out(~isnan(interp_nh4)),interp_nh4(~isnan(interp_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['GY 2001 daily interped nh4 vs. transport']);

figure;
plot(sj_trans_out(~isnan(raw_nh4)),raw_nh4(~isnan(raw_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['GY 2001 daily nh4 vs. transport']);



clearvars r lags
% [r,lags] = xcorr(sj_trans_out_x(~isnan(interp_nh4)),interp_nh4_x(~isnan(interp_nh4)),'normalized') 
[r,lags] = xcorr(sj_trans_out_x(~isnan(raw_nh4)),raw_nh4_x(~isnan(raw_nh4)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily nh4 vs. transport']);
n_nan_nh4 = find(isnan(raw_nh4)==0)



text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL nh4 vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');


%% no3
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(raw_no3,'b*','linew',2);
        plot(interp_no3,'r','linew',2);
%         plot(1:length(obm_gy_no3),obm_gy_no3 + obs_std_gy_no3,'m-','linew',2);
%         plot(1:length(obm_gy_no3),obm_gy_no3 - obs_std_gy_no3,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL no3 + transport']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(raw_chl)])
        set(gca,'xticklabel',2007:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3'),'-dpng') 
yyaxis right
plot(sj_trans_out,'g-','linew',2);
xlim([1 length(raw_no3)])

return

corrcoef(sj_trans_out(~isnan(interp_no3)),interp_no3(~isnan(interp_no3)))
% corrcoef(sj_trans_out,interp_no3,'Rows','complete')
corrcoef(sj_trans_out(~isnan(raw_no3)),raw_no3(~isnan(raw_no3)))

figure;
plot(sj_trans_out(~isnan(interp_no3)),interp_no3(~isnan(interp_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['GY 2001 daily interped no3 vs. transport']);

figure;
plot(sj_trans_out(~isnan(raw_no3)),raw_no3(~isnan(raw_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['GY 2001 daily no3 vs. transport']);



clearvars r lags
% [r,lags] = xcorr(sj_trans_out_x(~isnan(interp_no3)),interp_no3_x(~isnan(interp_no3)),'normalized') 
[r,lags] = xcorr(sj_trans_out_x(~isnan(raw_no3)),raw_no3_x(~isnan(raw_no3)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily no3 vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL no3 vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');



%% do
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(raw_do,'b*','linew',2);
        plot(interp_do,'r','linew',2);
%         plot(1:length(obm_gy_do),obm_gy_do + obs_std_gy_do,'m-','linew',2);
%         plot(1:length(obm_gy_do),obm_gy_do - obs_std_gy_do,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL do + transport']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('do (mmol / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(raw_chl)])
        set(gca,'xticklabel',2007:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_do'),'-dpng') 
yyaxis right
plot(sj_trans_out,'g-','linew',2);
xlim([1 length(raw_do)])

return

corrcoef(sj_trans_out(~isnan(interp_do)),interp_do(~isnan(interp_do)))
% corrcoef(sj_trans_out,interp_do,'Rows','complete')
corrcoef(sj_trans_out(~isnan(raw_do)),raw_do(~isnan(raw_do)))

figure;
plot(sj_trans_out(~isnan(interp_do)),interp_do(~isnan(interp_do)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('do (mmol/ m^3)','fontsize',13)
grid on;
title(['GY 2001 daily interped do vs. transport']);

figure;
plot(sj_trans_out(~isnan(raw_do)),raw_do(~isnan(raw_do)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('do (mmol/ m^3)','fontsize',13)
grid on;
title(['GY 2001 daily do vs. transport']);


clearvars r lags
% [r,lags] = xcorr(sj_trans_out_x(~isnan(interp_do)),interp_do_x(~isnan(interp_do)),'normalized') 
[r,lags] = xcorr(sj_trans_out_x(~isnan(raw_do)),raw_do_x(~isnan(raw_do)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily do vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL do vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');



%% po4
%% SS
corrcoef(monthly_trans(~isnan(mon_in_ss)),mon_in_ss(~isnan(mon_in_ss)))
% corrcoef(monthly_trans,mon_in_ss,'Rows','complete')
corrcoef(monthly_trans(~isnan(monthly_ss)),monthly_ss(~isnan(monthly_ss)))
corrcoef(yr_trans(~isnan(yr_ss)),yr_ss(~isnan(yr_ss)))

figure;
plot(monthly_trans(~isnan(mon_in_ss)),mon_in_ss(~isnan(mon_in_ss)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('ss (mg/L)','fontsize',13)
grid on;
title(['monthly interped ss vs. transport']);

figure;
plot(monthly_trans(~isnan(monthly_ss)),monthly_ss(~isnan(monthly_ss)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('ss (mg/L)','fontsize',13)
grid on;
title(['monthly ss vs. transport']);

figure;
plot(yr_trans(~isnan(yr_ss)),yr_ss(~isnan(yr_ss)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('ss (mg/L)','fontsize',13)
grid on;
title(['yearly ss vs. transport']);
ylim([0 11]);  xlim([0 140])

%% chl

corrcoef(monthly_trans(~isnan(mon_in_chl)),mon_in_chl(~isnan(mon_in_chl)))
% corrcoef(monthly_trans,mon_in_chl,'Rows','complete')
corrcoef(monthly_trans(~isnan(monthly_chl)),monthly_chl(~isnan(monthly_chl)))
corrcoef(yr_trans(~isnan(yr_chl)),yr_chl(~isnan(yr_chl)))

figure;
plot(monthly_trans(~isnan(mon_in_chl)),mon_in_chl(~isnan(mon_in_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
grid on;
title(['monthly interped chl vs. transport']);

figure;
plot(monthly_trans(~isnan(monthly_chl)),monthly_chl(~isnan(monthly_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
grid on;
title(['monthly chl vs. transport']);

figure;
plot(yr_trans(~isnan(yr_chl)),yr_chl(~isnan(yr_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
grid on;
title(['yearly chl vs. transport']);
ylim([0 11]);  xlim([0 140])

%% nh4

corrcoef(monthly_trans(~isnan(mon_in_nh4)),mon_in_nh4(~isnan(mon_in_nh4)))
% corrcoef(monthly_trans,mon_in_nh4,'Rows','complete')
corrcoef(monthly_trans(~isnan(monthly_nh4)),monthly_nh4(~isnan(monthly_nh4)))
corrcoef(yr_trans(~isnan(yr_nh4)),yr_nh4(~isnan(yr_nh4)))

figure;
plot(monthly_trans(~isnan(mon_in_nh4)),mon_in_nh4(~isnan(mon_in_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly interped nh4 vs. transport']);

figure;
plot(monthly_trans(~isnan(monthly_nh4)),monthly_nh4(~isnan(monthly_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly nh4 vs. transport']);


figure;
plot(yr_trans(~isnan(yr_nh4)),yr_nh4(~isnan(yr_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly nh4 vs. transport']);
xlim([0 140])


%% no3

corrcoef(monthly_trans(~isnan(mon_in_no3)),mon_in_no3(~isnan(mon_in_no3)))
% corrcoef(monthly_trans,mon_in_no3,'Rows','complete')
corrcoef(monthly_trans(~isnan(monthly_no3)),monthly_no3(~isnan(monthly_no3)))
corrcoef(yr_trans(~isnan(yr_no3)),yr_no3(~isnan(yr_no3)))

figure;
plot(monthly_trans(~isnan(mon_in_no3)),mon_in_no3(~isnan(mon_in_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly interped no3 vs. transport']);

figure;
plot(monthly_trans(~isnan(monthly_no3)),monthly_no3(~isnan(monthly_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly no3 vs. transport']);


figure;
plot(yr_trans(~isnan(yr_no3)),yr_no3(~isnan(yr_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly no3 vs. transport']);
xlim([0 140])


%% po4

corrcoef(monthly_trans(~isnan(mon_in_po4)),mon_in_po4(~isnan(mon_in_po4)))
% corrcoef(monthly_trans,mon_in_po4,'Rows','complete')
corrcoef(monthly_trans(~isnan(monthly_po4)),monthly_po4(~isnan(monthly_po4)))
corrcoef(yr_trans(~isnan(yr_po4)),yr_po4(~isnan(yr_po4)))

figure;
plot(monthly_trans(~isnan(mon_in_po4)),mon_in_po4(~isnan(mon_in_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly interped po4 vs. transport']);

figure;
plot(monthly_trans(~isnan(monthly_po4)),monthly_po4(~isnan(monthly_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly po4 vs. transport']);


figure;
plot(yr_trans(~isnan(yr_po4)),yr_po4(~isnan(yr_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly po4 vs. transport']);
xlim([0 140])


corrcoef(monthly_trans(~isnan(mon_in_po4)),mon_in_po4(~isnan(mon_in_po4)))
% corrcoef(monthly_trans,mon_in_po4,'Rows','complete')
corrcoef(monthly_trans(~isnan(monthly_po4)),monthly_po4(~isnan(monthly_po4)))
corrcoef(yr_trans(~isnan(yr_po4)),yr_po4(~isnan(yr_po4)))

figure;
plot(monthly_trans(~isnan(mon_in_po4)),mon_in_po4(~isnan(mon_in_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly interped po4 vs. transport']);

figure;
plot(monthly_trans(~isnan(monthly_po4)),monthly_po4(~isnan(monthly_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly po4 vs. transport']);


figure;
plot(yr_trans(~isnan(yr_po4)),yr_po4(~isnan(yr_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly po4 vs. transport']);



%% koem chl vs river chl

corrcoef(koem.yr_chl(~isnan(yr_chl)),yr_chl(~isnan(yr_chl)))
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(yr_chl,'b*','linew',2);
        plot(yr_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS chl vs. chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:30);
        ylim([0 12])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(koem.yr_chl,'b-','linew',2);
ylim([0 12])
xtickangle(45)
legend('river','koem');

figure;
plot(koem.yr_chl(~isnan(yr_chl)),yr_chl(~isnan(yr_chl)),'*')
xlabel('chl (ug/L)','fontsize',13)
ylabel('koem-chl (ug/L)','fontsize',13)
grid on;
title(['yearly po4 vs. transport']);
xlim([0 12])
ylim([0 12])


clearvars r lags
% [r,lags] = xcorr(monthly_trans_x(~isnan(mon_in_po4)),mon_in_po4_x(~isnan(mon_in_po4)),'normalized') 
[r,lags] = xcorr(sj_trans_out_x(~isnan(monthly_po4)),monthly_po4_x(~isnan(monthly_po4)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily po4 vs. transport']);
n_nan_po4 = find(isnan(monthly_po4)==0)



text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL po4 vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% June concentration vs. transp. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ss
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(june_ss,'b*','linew',2);
        plot(june_in_ss,'r','linew',2);
%         plot(1:length(obm_gy_ss),obm_gy_ss + obs_std_gy_ss,'m-','linew',2);
%         plot(1:length(obm_gy_ss),obm_gy_ss - obs_std_gy_ss,'m-','linew',2);
        title(['GY 1989-2018 june OBS transp vs. ss']);
        xlabel('time','fontsize',13)
        ylabel('ss (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[1:length(1989:2018)]);
        xlim([1 length(1989:2018)])
        set(gca,'xticklabel',1989:2018,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_ss'),'-dpng') 
yyaxis right
plot(june_trans,'g-','linew',2); xtickangle(45)

% CHL
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(june_chl,'b*','linew',2);
        plot(june_in_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 june OBS transp vs. chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[1:length(1989:2018)]);
        xlim([1 length(1989:2018)])
        set(gca,'xticklabel',1989:2018,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(june_trans,'g-','linew',2);
% for i = 6:12:360; xline(i,'color','m'); end

return

corrcoef(june_trans(~isnan(june_in_chl)),june_in_chl(~isnan(june_in_chl)))
% corrcoef(monthly_trans,interp_chl,'Rows','complete')
corrcoef(june_trans(~isnan(june_chl)),june_chl(~isnan(june_chl)))

figure;
plot(june_trans(~isnan(june_in_chl)),june_in_chl(~isnan(june_in_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('Chla (ug/L)','fontsize',13)
grid on;
title(['GY 1989-2018 june interped chl vs. transport']);

% NH4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(june_nh4,'b*','linew',2);
        plot(june_in_nh4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 june OBS transp vs. nh4']);
        xlabel('time','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[1:length(1989:2018)]);
        xlim([1 length(1989:2018)])
        set(gca,'xticklabel',1989:2018,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(june_trans,'g-','linew',2);

corrcoef(june_trans(~isnan(june_in_nh4)),june_in_nh4(~isnan(june_in_nh4)))
% corrcoef(sj_trans_out,interp_nh4,'Rows','complete')
corrcoef(june_trans(~isnan(june_nh4)),june_nh4(~isnan(june_nh4)))

figure;
plot(june_trans(~isnan(june_in_nh4)),june_in_nh4(~isnan(june_in_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N /m^3)','fontsize',13)
grid on;
title(['GY 1989-2018 june interped nh4 vs. transport']);

% NO3
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(june_no3,'b*','linew',2);
        plot(june_in_no3,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 june OBS transp vs. no3']);
        xlabel('time','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[1:length(1989:2018)]);
        xlim([1 length(1989:2018)])
        set(gca,'xticklabel',1989:2018,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(june_trans,'g-','linew',2);

corrcoef(june_trans(~isnan(june_in_no3)),june_in_no3(~isnan(june_in_no3)))
% corrcoef(sj_trans_out,interp_no3,'Rows','complete')
corrcoef(june_trans(~isnan(june_no3)),june_no3(~isnan(june_no3)))

figure;
plot(june_trans(~isnan(june_in_no3)),june_in_no3(~isnan(june_in_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N /m^3)','fontsize',13)
grid on;
title(['GY 1989-2018 june interped no3 vs. transport']);

% DO
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(june_do,'b*','linew',2);
        plot(june_in_do,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 june OBS transp vs. do']);
        xlabel('time','fontsize',13)
        ylabel('do (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[1:length(1989:2018)]);
        xlim([1 length(1989:2018)])
        set(gca,'xticklabel',1989:2018,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(june_trans,'g-','linew',2);

corrcoef(june_trans(~isnan(june_in_do)),june_in_do(~isnan(june_in_do)))
% corrcoef(sj_trans_out,interp_do,'Rows','complete')
corrcoef(june_trans(~isnan(june_do)),june_do(~isnan(june_do)))

figure;
plot(june_trans(~isnan(june_in_do)),june_in_do(~isnan(june_in_do)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('do (mmol/m^3)','fontsize',13)
grid on;
title(['GY 1989-2018 june interped do vs. transport']);




%% total mass through the river

figure;
plot((monthly_trans) .* mon_in_chl)  % m^3/s * mg/m^3 = mg/s
hold; 
yyaxis right
plot(monthly_trans,'r');


figure;
plot((monthly_trans) .* (mon_in_nh4 ./1000 .* 14 ))  % m^3/s * g/m^3 = g/s
hold; 
yyaxis right
plot(monthly_trans,'r');


figure;
plot((monthly_trans) .* (mon_in_no3 ./1000 .* 14 ))  % m^3/s * g/m^3 = g/s
hold; 
yyaxis right
plot(monthly_trans,'r');

% figure;
% plot((monthly_trans) .* (mon_in_no3 ./1000 .* 14 ))  % m^3/s * g/m^3 = g/s
% hold; 
% yyaxis right
% plot(monthly_trans,'r');

cd D:\장기생태\Dynamic\06_river
load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','yp_w_*_04');   % polynomial_fitting_on_ecological_variable_to2004_new_v5_3sigma.m 
% ~2003 : advanced 전체기간
% 2004~ : 04~ data's climate
cd D:\장기생태\Dynamic\06_river\환경과학원
load('sumjin(songjung)_polynomial_climate_to2004_3sig.mat','yp_w_*_af'); % ~2004 : 1st regime, 2004~2018 : 2nd regime
% sumjin_do = yp_w_do_04'*0.7*44.661;
% sumjin_chl = yp_w_chl_04';
% sumjin_no3 = yp_w_no3_04'*1000/14;
% sumjin_nh4 = yp_w_nh4_04'*1000/14;

clearvars all_recon_*
for i = 1:30
   yy = t_year(i);
   clearvars tem_*
   if i >= 15  %1989~2003
   tem_do = yp_w_do_04;
   tem_chl = yp_w_chl_04;
   tem_no3 = yp_w_no3_04;
   tem_nh4 = yp_w_nh4_04;
   else %2004~
   tem_do = yp_w_do_af;
   tem_chl = yp_w_chl_af;
   tem_no3 = yp_w_no3_af;
   tem_nh4 = yp_w_nh4_af; 
   end
if leapyear(yy) == 0
   tem_do(60) = [];
   tem_chl(60) = [];
   tem_no3(60) = [];
   tem_nh4(60) = [];
end
if i == 1
    all_recon_do = tem_do;
    all_recon_chl = tem_chl;
    all_recon_no3 = tem_no3;
    all_recon_nh4 = tem_nh4;
else
   all_recon_do=cat(2,all_recon_do,tem_do);
   all_recon_chl=cat(2,all_recon_chl,tem_chl);
   all_recon_no3=cat(2,all_recon_no3,tem_no3);
   all_recon_nh4=cat(2,all_recon_nh4,tem_nh4);
end

end

% daily transport of nutrients too many data for plot it
figure;hold; 
plot((sj_trans_out) .* all_recon_chl,'b')  % m^3/s * mg/m^3 = mg/s
hold on;
yyaxis right
plot(sj_trans_out,'r');
hold off

figure;hold; 
plot(((sj_trans_out) .* (all_recon_nh4))./1000,'b')  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on;
yyaxis right
plot(sj_trans_out,'r');
hold off


figure;hold; 
plot(((sj_trans_out) .* (all_recon_no3))./1000,'b');  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
hold on;
% g/s /1000 = kg/s
yyaxis right
plot(sj_trans_out,'r'); hold off

% figure;
% plot((sj_trans_out) .* (mon_in_no3 ./1000 .* 14 ))  % m^3/s * mg/m^3 = mg/s
% hold; 
% yyaxis right
% plot(sj_trans_out,'r');

% save('transp_daily_and_nutrients.mat','all_recon_*','sj_trans_out');

% daily on 2001 only

figure;hold; 
plot((dis_pre_total{2001-1979}) .* yp_w_chl_04,'color','b','linew',2)  % m^3/s * mg/m^3 = mg/s
ylabel('Chla x transport (mg/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
hold off

figure;hold; 
plot((dis_pre_total{2001-1979}) .* yp_w_nh4_04,'color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
ylabel('nh4 x transport (g/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
hold off


figure;hold; 
plot((dis_pre_total{2001-1979}) .* yp_w_no3_04,'color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
ylabel('no3 x transport (g/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
hold off


tt_tick = [1 t_indx(1:11)+1];
figure;hold; 
for i = 8:30
plot(double(dis_pre_total{t_year(i)-1979}),'b'); hold on;  % m^3/s * mg/m^3 = mg/s
end
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title('1996~2018 songjung transport')
% set(gca,'xtick',[tt_tick(1:end)]);
% set(gca,'xticklabel',,'fontsize',10);


%% find when 6 month's transp. is larger than summer

for i = 8:30
figure; 
plot(double(dis_pre_total{t_year(i)-1979}),'b');  % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_year(i)) ' songjung transport'])
saveas(gcf,[num2str(t_year(i)) '-songjung-transport.png']);
close all;
end

for i = 8:30
figure; 
plot(double(dis_pre_total{t_year(i)-1979}),'b');  % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_year(i)) ' songjung transport'])
ylim([0 5000])
saveas(gcf,[num2str(t_year(i)) '-songjung-transport.png']);
close all;
end

cd D:\장기생태\Dynamic\06_river\환경과학원
clearvars raw_*
load('sumjin(jinwall)_polynomial_climate_advanced(v4)_3sig.mat','raw_*'); 

r_do_jinwal = raw_do *0.7*44.661;
r_chl_jinwal = raw_chl;
r_no3_jinwal = raw_no3 .* 1000 ./14;
r_nh4_jinwal = raw_nh4 .* 1000 ./14;

t_07=load('2007_2019_time_axis.mat','t_tick');
t_07_tick = t_07.t_tick;

cd D:\장기생태\Dynamic\06_river
clearvars raw_*
load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat'); 
r_do_song = raw_do *0.7*44.661;
r_chl_song = raw_chl;
r_no3_song = raw_no3 .* 1000 ./14;
r_nh4_song = raw_nh4 .* 1000 ./14;

t_indx_year = t_indx(12:12:end); %1989~2019

% when 6 month's transp. larger
% 1996 : 8th year
% 1998 : 10th year
% 2001 : 13th year
% 2008 : 20th year 
nut_big_y_indx = [8, 10, 13, 20];
t_big_year = [1996, 1998, 2001, 2008];
for i = 1:4
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS nh4'])
ylim([0 5000])
yyaxis right
plot(r_nh4_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('NH4 (mmol N /m^3)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-nh4.png']);
close all;
end

for i = 1:4
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS no3'])
ylim([0 5000])
yyaxis right
plot(r_no3_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('NO3 (mmol N /m^3)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-no3.png']);
close all;
end

for i = 1:4
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS Chla'])
ylim([0 5000])
yyaxis right
plot(r_chl_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('Chla (ug/L)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-Chla.png']);
close all;
end

%%
figure;hold;
clearvars temp_val
temp_val = yp_w_chl_04;
temp_val(60)=[];
plot((dis_pre_total{2001-1979}) .* temp_val','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s
hold on; 
plot((dis_pre_total{2001-1979}) .* (r_chl_song(t_indx_year(nut_big_y_indx(3)-1)+1:t_indx_year(nut_big_y_indx(3))) )','r*','linew',2)
ylabel('Chla x transport (mg/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off

figure;hold; 
clearvars temp_val
temp_val = yp_w_nh4_04;
temp_val(60)=[];
plot((dis_pre_total{2001-1979}) .* temp_val','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2001-1979}) .* (r_nh4_song(t_indx_year(nut_big_y_indx(3)-1)+1:t_indx_year(nut_big_y_indx(3))) ./ 1000 .*14)','r*','linew',2)
ylabel('nh4 x transport (g/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off


figure;hold; 
clearvars temp_val
temp_val = yp_w_no3_04;
temp_val(60)=[];
plot((dis_pre_total{2001-1979}) .* temp_val','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2001-1979}) .* (r_no3_song(t_indx_year(nut_big_y_indx(3)-1)+1:t_indx_year(nut_big_y_indx(3))) ./ 1000 .*14)','r*','linew',2)
ylabel('no3 x transport (g/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off

%% 2008
figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_chl_af','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_chl_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) )','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_chl_jinwal(t_07_tick(2)+1:t_07_tick(3)) )','c*','linew',2)
ylabel('Chla x transport (mg/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off

figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_nh4_af','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_nh4_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) ./ 1000 .*14)','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_nh4_jinwal(t_07_tick(2)+1:t_07_tick(3)) ./ 1000 .*14)','c*','linew',2)
ylabel('nh4 x transport (g/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off


figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_no3_af','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_no3_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) ./ 1000 .*14)','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_no3_jinwal(t_07_tick(2)+1:t_07_tick(3)) ./ 1000 .*14)','c*','linew',2)
ylabel('no3 x transport (g/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off



% set(gca,'xtick',[tt_tick(1:end)]);
% set(gca,'xticklabel',,'fontsize',10);

t_indx_year = t_indx(12:12:end); %1989~2019
%nh4, no3 exist from 1996

nut_big_y_indx = 8:length(t_indx_year)-1;
t_big_year = 1996:2018;
for i = 1:23
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS nh4'])
ylim([0 5000])
yyaxis right
plot(r_nh4_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('NH4 (mmol N /m^3)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-nh4.png']);
close all;
end

for i = 1:23
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS no3'])
ylim([0 5000])
yyaxis right
plot(r_no3_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('NO3 (mmol N /m^3)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-no3.png']);
close all;
end

for i = 1:23
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS Chla'])
ylim([0 5000])
yyaxis right
plot(r_chl_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('Chla (ug/L)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-Chla.png']);
close all;
end

%% + jinwall all plot
%% since 2007 : 19th year 
figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_chl_af','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_chl_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) )','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_chl_jinwal(t_07_tick(2)+1:t_07_tick(3)) )','c*','linew',2)
ylabel('Chla x transport (mg/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off

figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_nh4_af','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_nh4_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) ./ 1000 .*14)','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_nh4_jinwal(t_07_tick(2)+1:t_07_tick(3)) ./ 1000 .*14)','c*','linew',2)
ylabel('nh4 x transport (g/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off


figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_no3_af','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_no3_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) ./ 1000 .*14)','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_no3_jinwal(t_07_tick(2)+1:t_07_tick(3)) ./ 1000 .*14)','c*','linew',2)
ylabel('no3 x transport (g/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off


t_indx_year = t_indx(12:12:end); %1989~2019
%nh4, no3 exist from 1996
nut_big_y_indx = 19:length(t_indx_year)-1;
t_big_year = 2007:2018;
for i = 1:12
figure; 
clearvars yp_w_*_af2
if leapyear(t_big_year(i))==1
     yp_w_chl_af2 = yp_w_chl_af;  
else
     yp_w_chl_af2 = yp_w_chl_af;
     yp_w_chl_af2(60)=[];
end
    plot((dis_pre_total{t_big_year(i)-1979}) .* yp_w_chl_af2','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s

hold on;
plot((dis_pre_total{t_big_year(i)-1979}) .* (r_chl_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))) )','r*','linew',2)
hold on;
if i == 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_chl_jinwal(t_07_tick(i):t_07_tick(i+1)))','c*','linew',2)
elseif i ~= 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_chl_jinwal(t_07_tick(i)+1:t_07_tick(i+1)))','c*','linew',2)
end
hold on; grid on;
ylabel('Chla x transport (mg/s)')
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS Chla'])
yyaxis right
plot(double(dis_pre_total{t_big_year(i)-1979}),'linew',2); % m^3/s * mg/m^3 = mg/s
ylim([0 5000])
ylabel('tansport(m^3/s)');
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-Chla-jin.png']);
close all;
end

t_indx_year = t_indx(12:12:end); %1989~2019
%nh4, no3 exist from 1996
nut_big_y_indx = 19:length(t_indx_year)-1;
t_big_year = 2007:2018;
for i = 1:12
figure; 
clearvars yp_w_*_af2
if leapyear(t_big_year(i))==1
     yp_w_nh4_af2 = yp_w_nh4_af;  
else
     yp_w_nh4_af2 = yp_w_nh4_af;
     yp_w_nh4_af2(60)=[];
end
    plot((dis_pre_total{t_big_year(i)-1979}) .* yp_w_nh4_af2','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s

hold on;
plot((dis_pre_total{t_big_year(i)-1979}) .* (r_nh4_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))) ./ 1000 .* 14 )','r*','linew',2)
hold on;
if i == 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_nh4_jinwal(t_07_tick(i):t_07_tick(i+1))./ 1000 .* 14)','c*','linew',2)
elseif i ~= 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_nh4_jinwal(t_07_tick(i)+1:t_07_tick(i+1))./ 1000 .* 14)','c*','linew',2)
end
hold on;grid on;
ylim([0 300])
ylabel('nh4 x transport (mg/s)')
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS nh4'])
yyaxis right
plot(double(dis_pre_total{t_big_year(i)-1979}),'linew',2); % m^3/s * mg/m^3 = mg/s
ylim([0 5000])
ylabel('tansport(m^3/s)');
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-nh4-jin.png']);
close all;
end


t_indx_year = t_indx(12:12:end); %1989~2019
%no3, no3 exist from 1996
nut_big_y_indx = 19:length(t_indx_year)-1;
t_big_year = 2007:2018;
for i = 1:12
figure; 
clearvars yp_w_*_af2
if leapyear(t_big_year(i))==1
     yp_w_no3_af2 = yp_w_no3_af;  
else
     yp_w_no3_af2 = yp_w_no3_af;
     yp_w_no3_af2(60)=[];
end
    plot((dis_pre_total{t_big_year(i)-1979}) .* yp_w_no3_af2','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s

hold on;
plot((dis_pre_total{t_big_year(i)-1979}) .* (r_no3_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))) ./ 1000 .* 14 )','r*','linew',2)
hold on;
if i == 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_no3_jinwal(t_07_tick(i):t_07_tick(i+1))./ 1000 .* 14)','c*','linew',2)
elseif i ~= 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_no3_jinwal(t_07_tick(i)+1:t_07_tick(i+1))./ 1000 .* 14)','c*','linew',2)
end
hold on; grid on;
ylabel('no3 x transport (mg/s)')
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS no3'])
yyaxis right
plot(double(dis_pre_total{t_big_year(i)-1979}),'linew',2); % m^3/s * mg/m^3 = mg/s
ylim([0 5000])
ylabel('tansport(m^3/s)');
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-no3-jin.png']);
close all;
end


