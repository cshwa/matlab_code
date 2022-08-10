close all; clear; clc;
% -------------------------------------------------------------------------
t_year = 2007:2018
% _________________________________________________________________________
cd D:\장기생태\Dynamic\06_river
load jinwall_raw_data_daily.mat %there has to be song vs. song_transp
%jinwall vs. jinwall_transp
load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','raw_*')

cd D:\장기생태\Dynamic\06_river\data\sj_1980to1996
%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
load songjung_discharge_1980to2018.mat  % pre_merg_dis is discharge

cd D:\장기생태\Dynamic\06_river
% approximation : Sumjin(songjung) & Namgang have same discharge during 1980~1989.

load('sumjin_recons_water_temp_present.mat','merg_recon_w_c');
sumjin_re_w_c = merg_recon_w_c;

clearvars temperature sj_trans_out
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_temp=sumjin_re_w_c{t_year(order_i)-1989};
tempo_trans = dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
    temperature=tempo_temp;
    sj_trans_out=tempo_trans;
else
    temperature=cat(1,temperature,tempo_temp);
    sj_trans_out=cat(1,sj_trans_out,tempo_trans);
end
end


raw_chl = raw_chl;
raw_do = raw_do *0.7*44.661;
raw_no3 = raw_no3 .*1000 ./14 ;
raw_nh4 = raw_nh4 .*1000 ./14 ;

interp_chl=raw_chl;
interp_do=raw_do;
interp_no3=raw_no3;
interp_nh4=raw_nh4;

t=1:length(raw_chl);
interp_chl(isnan(interp_chl))=interp1(t(~isnan(interp_chl)),interp_chl(~isnan(interp_chl)),t(isnan(interp_chl)));
interp_no3(isnan(interp_no3))=interp1(t(~isnan(interp_no3)),interp_no3(~isnan(interp_no3)),t(isnan(interp_no3)));
interp_nh4(isnan(interp_nh4))=interp1(t(~isnan(interp_nh4)),interp_nh4(~isnan(interp_nh4)),t(isnan(interp_nh4)));
interp_do(isnan(interp_do))=interp1(t(~isnan(interp_do)),interp_do(~isnan(interp_do)),t(isnan(interp_do)));

sj_trans_out_x = sj_trans_out - nanmean(sj_trans_out);
interp_chl_x = interp_chl - nanmean(interp_chl);
interp_do_x = interp_do - nanmean(interp_do);
interp_no3_x = interp_no3 - nanmean(interp_no3);
interp_nh4_x = interp_nh4 - nanmean(interp_nh4);

raw_chl_x = raw_chl - nanmean(raw_chl);
raw_do_x = raw_do - nanmean(raw_do);
raw_no3_x = raw_no3 - nanmean(raw_no3);
raw_nh4_x = raw_nh4 - nanmean(raw_nh4);

cd D:\장기생태\Dynamic\06_river\환경과학원
load 2007_2019_time_axis.mat
t_tick(2:end) = t_tick(2:end)+1;
%% chl
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(raw_chl,'b*','linew',2);
        plot(interp_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 2001 daily OBS transp vs. chl']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(raw_chl)])
        set(gca,'xticklabel',2007:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(sj_trans_out,'g-','linew',2);

return

corrcoef(sj_trans_out(~isnan(interp_chl)),interp_chl(~isnan(interp_chl)))
% corrcoef(sj_trans_out,interp_chl,'Rows','complete')
corrcoef(sj_trans_out(~isnan(raw_chl)),raw_chl(~isnan(raw_chl)))

figure;
plot(sj_trans_out(~isnan(interp_chl)),interp_chl(~isnan(interp_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('Chla (ug/L)','fontsize',13)
grid on;
title(['GY 2001 daily interped chl vs. transport']);

figure;
plot(sj_trans_out(~isnan(raw_chl)),raw_chl(~isnan(raw_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('Chla (ug/L)','fontsize',13)
grid on;
title(['GY 2001 daily chl vs. transport']);


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


