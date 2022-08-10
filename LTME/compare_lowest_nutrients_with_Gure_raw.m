close all; clear; clc;
cd D:\장기생태\Dynamic\06_river\환경과학원
load('sumjin(jinwall)_polynomial_climate_advanced(v4)_3sig.mat','raw_*'); 

r_do_jinwal = raw_do *0.7*44.661;
r_chl_jinwal = raw_chl;
r_no3_jinwal = raw_no3 .* 1000 ./14;
r_nh4_jinwal = raw_nh4 .* 1000 ./14;

%% vs.
cd D:\장기생태\Dynamic\06_river
clearvars raw_*
load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat'); 
r_do_song = raw_do *0.7*44.661;
r_chl_song = raw_chl;
r_no3_song = raw_no3 .* 1000 ./14;
r_nh4_song = raw_nh4 .* 1000 ./14;

% % % cd D:\장기생태\Dynamic\06_river\환경과학원
% % % clearvars yp_w_*_04 yp_w_*_af 
% % % load('sumjin(gure)_polynomial_climate_to2004_3sig.mat','yp_w_*'); % ~2004 : 1st regime, 2004~2018 : 2nd regime
% % % clearvars yp_w_do yp_w_chl yp_w_no3 yp_w_nh4
% % % load('sumjin(gure)_polynomial_climate_to2004_v2_3sig.mat','yp_w_*'); % ~2018 : whole regime

r_do_gure_match = r_do_song(length(r_do_song)-length(r_do_jinwal)+1:end);
r_chl_gure_match = r_chl_song(length(r_chl_song)-length(r_chl_jinwal)+1:end);
r_no3_gure_match = r_no3_song(length(r_no3_song)-length(r_no3_jinwal)+1:end);
r_nh4_gure_match = r_nh4_song(length(r_nh4_song)-length(r_nh4_jinwal)+1:end);

cd D:\장기생태\Dynamic\06_river\환경과학원
load 2007_2019_time_axis.mat
t_tick(2:end) = t_tick(2:end)+1; %time axes for years

% DO
figure; hold on;
plot(r_do_jinwal,'*','color','b','linew',1)
plot(r_do_gure_match,'o','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('DO (mmol/m^3)','fontsize',13)
title('jinwall vs. gure raw OBS DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'fontsize',13)
set(gca,'xtick',[t_tick(1:end)]);
xlim([1 length(r_do_jinwal)])
set(gca,'xticklabel',2007:2019,'fontsize',10);
legend('jinwall','gure');

clearvars yp_est_do
figure; hold on;
plot(r_do_gure_match,r_do_jinwal,'.','linew',2)
xp = r_do_gure_match;
pf_do = polyfit(r_do_gure_match,r_do_jinwal,1);
yp_est_do = polyval(pf_do,r_do_gure_match);
% plot(r_do_gure_match,yp_est_do,'r','linew',2);
xlabel('gure (mmol/m^3)','fontsize',13)
ylabel('jinwall (mmol/m^3)','fontsize',13)
title('gure vs. jinwall compare do','fontsize',13)
grid on;
set(gca,'fontsize',13)
legend('raw');


% figure; hold on;
% plot(y_do_jinwal,'r','linew',2); 
% plot(yp_est_do,'linew',2); % 04~ do
% yp_est_do_adv = polyval(pf_do,y_do_adv);
% plot(yp_est_do_adv,'g','linew',2); % ~04 do
% xlim([1 366]);
% xlabel('time(days)','fontsize',13)
% ylabel('do (mg/L)','fontsize',13)
% grid on;
% set(gca,'fontsize',13)
% legend('jinwall','jinwall-re','jinwall-re(past)');
% 
% figure; hold on;
% plot(1:366, yp_adv_do_jinwal ,'--','color','m','linew',2)
% plot(1:366, yp_adv_do ,'--','color','b','linew',2)
% plot(1:366, yp_w_do_af ,'--','color','r','linew',2)
% plot(1:366, yp_est_do_adv,'--','color','k','linew',2)
% plot(1:366, yp_est_do,'--','color','g','linew',2)
% xlabel('time(days)','fontsize',13)
% ylabel('do (mg/L)','fontsize',13)
% title('gure vs. jinwall + recon. do.','fontsize',13)
% grid on
% set(gca,'fontsize',13)
% % ylim([0 0.5])
% xlim([1 366])
% legend('jinwall','whole(adv)','2nd','jinwall-re(adv)','jinwall-re');

% CHL
figure; hold on;
plot(r_chl_jinwal,'*','color','b','linew',1)
plot(r_chl_gure_match,'o','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('Chla (ug/L)','fontsize',13)
title('jinwall vs. gure raw OBS Chla.','fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[t_tick(1:end)]);
xlim([1 length(r_chl_jinwal)])
set(gca,'xticklabel',2007:2019,'fontsize',10);
legend('jinwall','gure');

clearvars yp_est_chl
figure; hold on;
plot(r_chl_gure_match,r_chl_jinwal,'.','linew',2)
xp = r_chl_gure_match;
pf_chl = polyfit(r_chl_gure_match,r_chl_jinwal,1);
yp_est_chl = polyval(pf_chl,r_chl_gure_match);
% plot(r_chl_gure_match,yp_est_chl,'r','linew',2);
xlabel('gure (ug/L)','fontsize',13)
ylabel('jinwall (ug/L)','fontsize',13)
title('gure vs. jinwall compare chl','fontsize',13)
grid on;
set(gca,'fontsize',13)
legend('raw');

% % figure; hold on;
% % plot(y_chl_jinwal,'r','linew',2); 
% % plot(yp_est_chl,'linew',2); % 04~ chl
% % yp_est_chl_adv = polyval(pf_chl,y_chl_adv);
% % plot(yp_est_chl_adv,'g','linew',2); % ~04 chl
% % xlim([1 366]);
% % xlabel('time(days)','fontsize',13)
% % ylabel('chl (mg/m^3)','fontsize',13)
% % grid on;
% % set(gca,'fontsize',13)
% % legend('jinwall','jinwall-re','jinwall-re(past)');
% % 
% % figure; hold on;
% % plot(1:366, yp_adv_chl_jinwal ,'--','color','m','linew',2)
% % plot(1:366, yp_adv_chl ,'--','color','b','linew',2)
% % plot(1:366, yp_w_chl_af ,'--','color','r','linew',2)
% % plot(1:366, yp_est_chl_adv,'--','color','k','linew',2)
% % plot(1:366, yp_est_chl,'--','color','g','linew',2)
% % xlabel('time(days)','fontsize',13)
% % ylabel('chl (mg/m^3)','fontsize',13)
% % title('gure vs. jinwall + recon. chl.','fontsize',13)
% % grid on
% % set(gca,'fontsize',13)
% % % ylim([0 0.5])
% % xlim([1 366])
% % legend('jinwall','whole(adv)','2nd','jinwall-re(adv)','jinwall-re');

% no3
figure; hold on;
plot(r_no3_jinwal,'*','color','b','linew',1)
plot(r_no3_gure_match,'o','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (mmol N /m^3)','fontsize',13)
title('jinwall vs. gure raw OBS NO3.','fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[t_tick(1:end)]);
xlim([1 length(r_chl_jinwal)])
set(gca,'xticklabel',2007:2019,'fontsize',10);
legend('jinwall','gure');

clearvars yp_est_no3
figure; hold on;
plot(r_no3_gure_match,r_no3_jinwal,'.','linew',2)
xp = r_no3_gure_match;
pf_no3 = polyfit(r_no3_gure_match,r_no3_jinwal,1);
yp_est_no3 = polyval(pf_no3,r_no3_gure_match);
% plot(r_no3_gure_match,yp_est_no3,'r','linew',2);
xlabel('gure (mmol N /m^3)','fontsize',13)
ylabel('jinwall (mmol N /m^3)','fontsize',13)
title('gure vs. jinwall compare no3','fontsize',13)
grid on;
set(gca,'fontsize',13)
legend('raw');


% % figure; hold on;
% % plot(y_no3_jinwal,'r','linew',2); 
% % plot(yp_est_no3,'linew',2); % 04~ NO3
% % yp_est_no3_adv = polyval(pf_no3,y_no3_adv);
% % plot(yp_est_no3_adv,'g','linew',2); % ~04 NO3
% % xlim([1 366]);
% % xlabel('time(days)','fontsize',13)
% % ylabel('NO3-N (mmol N /m^3)','fontsize',13)
% % grid on;
% % set(gca,'fontsize',13)
% % legend('jinwall','jinwall-re','jinwall-re(past)');
% % 
% % figure; hold on;
% % plot(1:366, yp_adv_no3_jinwal *1000/14,'--','color','m','linew',2)
% % plot(1:366, yp_adv_no3 *1000/14,'--','color','b','linew',2)
% % plot(1:366, yp_w_no3_af *1000/14,'--','color','r','linew',2)
% % plot(1:366, yp_est_no3_adv,'--','color','k','linew',2)
% % plot(1:366, yp_est_no3,'--','color','g','linew',2)
% % xlabel('time(days)','fontsize',13)
% % ylabel('no3-N (mmol N /m^3)','fontsize',13)
% % title('gure vs. jinwall + recon. no3.','fontsize',13)
% % grid on
% % set(gca,'fontsize',13)
% % % ylim([0 0.5])
% % xlim([1 366])
% % legend('jinwall','whole(adv)','2nd','jinwall-re(adv)','jinwall-re');

% nh4
figure; hold on;
plot(r_nh4_jinwal,'*','color','b','linew',1)
plot(r_nh4_gure_match,'o','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('NH4 (mmol N /m^3)','fontsize',13)
title('jinwall vs. gure raw OBS NH4.','fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[t_tick(1:end)]);
xlim([1 length(r_nh4_jinwal)])
set(gca,'xticklabel',2007:2019,'fontsize',10);
legend('jinwall','gure');


clearvars yp_est_nh4
figure; hold on;
plot(r_nh4_gure_match,r_nh4_jinwal,'.','linew',2)
xp = r_nh4_gure_match;
pf_nh4 = polyfit(r_nh4_gure_match,r_nh4_jinwal,1);
yp_est_nh4 = polyval(pf_nh4,r_nh4_gure_match);
% plot(r_nh4_gure_match,yp_est_nh4,'r','linew',2);
xlabel('gure (mmol N /m^3)','fontsize',13)
ylabel('jinwall (mmol N /m^3)','fontsize',13)
title('gure vs. jinwall compare nh4','fontsize',13)
grid on;
set(gca,'fontsize',13)
legend('raw');


% figure; hold on;
% plot(y_nh4_jinwal,'r','linew',2); 
% plot(yp_est_nh4,'linew',2); % 04~ NO3
% yp_est_nh4_adv = polyval(pf_nh4,y_nh4_adv);
% plot(yp_est_nh4_adv,'g','linew',2); % ~04 NO3
% xlim([1 366]);
% xlabel('time(days)','fontsize',13)
% ylabel('NH4-N (mmol N /m^3)','fontsize',13)
% grid on;
% set(gca,'fontsize',13)
% legend('jinwall','jinwall-re','jinwall-re(past)');
% 
% 
% figure; hold on;
% plot(1:366, yp_adv_nh4_jinwal *1000/14,'--','color','m','linew',2)
% plot(1:366, yp_adv_nh4 *1000/14,'--','color','b','linew',2)
% plot(1:366, yp_w_nh4_af *1000/14,'--','color','r','linew',2)
% plot(1:366, yp_est_nh4_adv,'--','color','k','linew',2)
% plot(1:366, yp_est_nh4,'--','color','g','linew',2)
% xlabel('time(days)','fontsize',13)
% ylabel('nh4-N (mmol N /m^3)','fontsize',13)
% title('gure vs. jinwall + recon. nh4.','fontsize',13)
% grid on
% set(gca,'fontsize',13)
% % ylim([0 0.5])
% xlim([1 366])
% legend('jinwall','whole(adv)','2nd','jinwall-re(adv)','jinwall-re');
