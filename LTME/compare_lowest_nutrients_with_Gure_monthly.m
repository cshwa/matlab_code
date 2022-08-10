close all; clear; clc;
cd D:\장기생태\Dynamic\06_river\환경과학원
load('jinwall_yymm_monthly_data_07to19_3sig.mat'); 

%monthly mean and filled missing value
mon_in_do_jin = mon_in_do *0.7*44.661;
mon_in_chl_jin = mon_in_chl;
mon_in_no3_jin = mon_in_no3 .* 1000 ./14;
mon_in_nh4_jin = mon_in_nh4 .* 1000 ./14;

%monthly mean only
monthly_do_jin = monthly_do *0.7*44.661;
monthly_chl_jin = monthly_chl;
monthly_no3_jin = monthly_no3 .* 1000 ./14;
monthly_nh4_jin = monthly_nh4 .* 1000 ./14;

%% vs.
cd D:\장기생태\Dynamic\06_river
clearvars -except *_jin
load('songjung_yymm_monthly_data_89to19_3sig.mat'); 

%monthly mean and filled missing value
mon_in_do_song = mon_in_do *0.7*44.661;
mon_in_chl_song = mon_in_chl;
mon_in_no3_song = mon_in_no3 .* 1000 ./14;
mon_in_nh4_song = mon_in_nh4 .* 1000 ./14;

%monthly mean only
monthly_do_song = monthly_do *0.7*44.661;
monthly_chl_song = monthly_chl;
monthly_no3_song = monthly_no3 .* 1000 ./14;
monthly_nh4_song = monthly_nh4 .* 1000 ./14;

% % % cd D:\장기생태\Dynamic\06_river\환경과학원
% % % clearvars yp_w_*_04 yp_w_*_af 
% % % load('sumjin(gure)_polynomial_climate_to2004_3sig.mat','yp_w_*'); % ~2004 : 1st regime, 2004~2018 : 2nd regime
% % % clearvars yp_w_do yp_w_chl yp_w_no3 yp_w_nh4
% % % load('sumjin(gure)_polynomial_climate_to2004_v2_3sig.mat','yp_w_*'); % ~2018 : whole regime

%% cut the songjung to make it same dim. with jinwall's data
%monthly mean and filled missing value
mon_in_do_song_mat = mon_in_do_song(length(mon_in_do_song)-length(mon_in_do_jin)+1:end);
mon_in_chl_song_mat = mon_in_chl_song(length(mon_in_chl_song)-length(mon_in_chl_jin)+1:end);
mon_in_no3_song_mat = mon_in_no3_song(length(mon_in_no3_song)-length(mon_in_no3_jin)+1:end);
mon_in_nh4_song_mat = mon_in_nh4_song(length(mon_in_nh4_song)-length(mon_in_nh4_jin)+1:end);

%monthly mean only
monthly_do_song_mat = monthly_do_song(length(monthly_do_song)-length(monthly_do_jin)+1:end);
monthly_chl_song_mat = monthly_chl_song(length(monthly_chl_song)-length(monthly_chl_jin)+1:end);
monthly_no3_song_mat = monthly_no3_song(length(monthly_no3_song)-length(monthly_no3_jin)+1:end);
monthly_nh4_song_mat = monthly_nh4_song(length(monthly_nh4_song)-length(monthly_nh4_jin)+1:end);

% DO
figure; hold on;
plot(monthly_do_jin,'*','color','b','linew',1)
plot(monthly_do_song_mat,'o','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('DO (mmol/m^3)','fontsize',13)
title('jinwall vs. gure monthly yymm OBS DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:12:length(monthly_do_jin)]);
set(gca,'xticklabel',2007:2019,'fontsize',10);
xlim([1 length(monthly_do_jin)])
legend('jinwall','gure');

clearvars yp_est_do nonan_idx
nonan_idx = find(isnan(monthly_do_song_mat + monthly_do_jin) == 0);
figure; hold on;
plot(monthly_do_song_mat(nonan_idx),monthly_do_jin(nonan_idx),'.','linew',2)
xp = monthly_do_song_mat(nonan_idx);
pf_do = polyfit(monthly_do_song_mat(nonan_idx),monthly_do_jin(nonan_idx),1);
yp_est_do = polyval(pf_do,monthly_do_song_mat(nonan_idx));
plot(monthly_do_song_mat(nonan_idx),yp_est_do,'r','linew',2);
xlabel('gure (mmol/m^3)','fontsize',13)
ylabel('jinwall (mmol/m^3)','fontsize',13)
title('gure vs. jinwall compare raw do','fontsize',13)
grid on;
set(gca,'fontsize',13)
% plot(mon_in_do_song_mat(1:nan_jin_idx(1)-1),yp_est_do_in,'g','linew',2);
% % there are no significant difference

% interped DO
figure; hold on;
plot(mon_in_do_jin,'-','color','b','linew',1)
plot(mon_in_do_song_mat,'-','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('DO (mmol/m^3)','fontsize',13)
title('jinwall vs. gure monthly interp. yymm OBS DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([6 17])
xlim([1 length(monthly_do_jin)])
legend('jinwall','gure');
set(gca,'xtick',[1:12:length(monthly_do_jin)]);
set(gca,'xticklabel',2007:2019,'fontsize',10);
xlim([1 length(monthly_do_jin)])

clearvars yp_est_do_in nan_jin_idx
nan_jin_idx = find(isnan(mon_in_do_jin) == 1)
find(isnan(mon_in_do_song_mat)==1)
figure; hold on;
plot(mon_in_do_song_mat(1:nan_jin_idx(1)-1),mon_in_do_jin(1:nan_jin_idx(1)-1),'.','linew',2)
xp = mon_in_do_song_mat(1:nan_jin_idx(1)-1);
pf_do_in = polyfit(mon_in_do_song_mat(1:nan_jin_idx(1)-1),mon_in_do_jin(1:nan_jin_idx(1)-1),1);
yp_est_do_in = polyval(pf_do_in,mon_in_do_song_mat(1:nan_jin_idx(1)-1));
plot(mon_in_do_song_mat(1:nan_jin_idx(1)-1),yp_est_do_in,'r','linew',2);
xlabel('gure (mmol/m^3)','fontsize',13)
ylabel('jinwall (mmol/m^3)','fontsize',13)
title('gure vs. jinwall compare monthly interp. yymm OBS DO','fontsize',13)
grid on;
set(gca,'fontsize',13)


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
% chl
figure; hold on;
plot(monthly_chl_jin,'*','color','b','linew',1)
plot(monthly_chl_song_mat,'o','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
title('jinwall vs. gure monthly yymm OBS chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:12:length(monthly_chl_jin)]);
set(gca,'xticklabel',2007:2019,'fontsize',10);
xlim([1 length(monthly_chl_jin)])
legend('jinwall','gure');

nan_jin_idx = find(isnan(monthly_chl_jin) == 1)
nan_song_idx = find(isnan(monthly_chl_song_mat)==1)

clearvars yp_est_chl nonan_idx
nonan_idx = find(isnan(monthly_chl_song_mat + monthly_chl_jin) == 0);
figure; hold on;
plot(monthly_chl_song_mat,monthly_chl_jin,'.','linew',2)
xp = monthly_chl_song_mat(nonan_idx);
pf_chl = polyfit(monthly_chl_song_mat(nonan_idx),monthly_chl_jin(nonan_idx),1);
yp_est_chl = polyval(pf_chl,monthly_chl_song_mat(nonan_idx));
plot(monthly_chl_song_mat(nonan_idx),yp_est_chl,'r','linew',2);
xlabel('gure (ug/L)','fontsize',13)
ylabel('jinwall (ug/L)','fontsize',13)
title('gure vs. jinwall compare chl','fontsize',13)
grid on;
set(gca,'fontsize',13)
% plot(mon_in_chl_song_mat(1:nan_jin_idx(1)-1),yp_est_chl_in,'g','linew',2);
% % there are no significant difference

% interped chl
figure; hold on;
plot(mon_in_chl_jin,'-','color','b','linew',1)
plot(mon_in_chl_song_mat,'-','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
title('jinwall vs. gure monthly interp. yymm OBS chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([6 17])
xlim([1 length(monthly_chl_jin)])
legend('jinwall','gure');
set(gca,'xtick',[1:12:length(monthly_chl_jin)]);
set(gca,'xticklabel',2007:2019,'fontsize',10);
xlim([1 length(monthly_chl_jin)])

clearvars yp_est_chl_in nan_jin_idx
nan_jin_idx = find(isnan(mon_in_chl_jin) == 1)
find(isnan(mon_in_chl_song_mat)==1)
figure; hold on;
plot(mon_in_chl_song_mat(1:nan_jin_idx(1)-1),mon_in_chl_jin(1:nan_jin_idx(1)-1),'.','linew',2)
xp = mon_in_chl_song_mat(1:nan_jin_idx(1)-1);
pf_chl_in = polyfit(mon_in_chl_song_mat(1:nan_jin_idx(1)-1),mon_in_chl_jin(1:nan_jin_idx(1)-1),1);
yp_est_chl_in = polyval(pf_chl_in,mon_in_chl_song_mat(1:nan_jin_idx(1)-1));
plot(mon_in_chl_song_mat(1:nan_jin_idx(1)-1),yp_est_chl_in,'r','linew',2);
xlabel('gure (ug/L)','fontsize',13)
ylabel('jinwall (ug/L)','fontsize',13)
title('gure vs. jinwall compare monthly interp. yymm OBS chl','fontsize',13)
grid on;
set(gca,'fontsize',13)

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
plot(monthly_no3_jin,'*','color','b','linew',1)
plot(monthly_no3_song_mat,'o','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('no3 (mmol N /m^3)','fontsize',13)
title('jinwall vs. gure monthly yymm OBS no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([6 17])
legend('jinwall','gure');
set(gca,'xtick',[1:12:length(monthly_no3_jin)]);
set(gca,'xticklabel',2007:2019,'fontsize',10);
xlim([1 length(monthly_no3_jin)])

nan_jin_idx = find(isnan(monthly_no3_jin) == 1)
nan_song_idx = find(isnan(monthly_no3_song_mat)==1)

clearvars yp_est_no3 nonan_idx
nonan_idx = find(isnan(monthly_no3_song_mat + monthly_no3_jin) == 0);
figure; hold on;
plot(monthly_no3_song_mat,monthly_no3_jin,'.','linew',2)
xp = monthly_no3_song_mat(nonan_idx);
pf_no3 = polyfit(monthly_no3_song_mat(nonan_idx),monthly_no3_jin(nonan_idx),1);
yp_est_no3 = polyval(pf_no3,monthly_no3_song_mat(nonan_idx));
plot(monthly_no3_song_mat(nonan_idx),yp_est_no3,'r','linew',2);
xlabel('gure (mmol N /m^3)','fontsize',13)
ylabel('jinwall (mmol N /m^3)','fontsize',13)
title('gure vs. jinwall compare no3','fontsize',13)
grid on;
set(gca,'fontsize',13)
% plot(mon_in_no3_song_mat(1:nan_jin_idx(1)-1),yp_est_no3_in,'g','linew',2);
% % there are no significant difference

% interped no3
figure; hold on;
plot(mon_in_no3_jin,'-','color','b','linew',1)
plot(mon_in_no3_song_mat,'-','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('no3 (mmol N /m^3)','fontsize',13)
title('jinwall vs. gure monthly interp. yymm OBS no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([6 17])
xlim([1 length(monthly_no3_jin)])
legend('jinwall','gure');
set(gca,'xtick',[1:12:length(monthly_no3_jin)]);
set(gca,'xticklabel',2007:2019,'fontsize',10);
xlim([1 length(monthly_no3_jin)])

clearvars yp_est_no3_in nan_jin_idx
nan_jin_idx = find(isnan(mon_in_no3_jin) == 1)
find(isnan(mon_in_no3_song_mat)==1)
figure; hold on;
plot(mon_in_no3_song_mat(1:nan_jin_idx(1)-1),mon_in_no3_jin(1:nan_jin_idx(1)-1),'.','linew',2)
xp = mon_in_no3_song_mat(1:nan_jin_idx(1)-1);
pf_no3_in = polyfit(mon_in_no3_song_mat(1:nan_jin_idx(1)-1),mon_in_no3_jin(1:nan_jin_idx(1)-1),1);
yp_est_no3_in = polyval(pf_no3_in,mon_in_no3_song_mat(1:nan_jin_idx(1)-1));
plot(mon_in_no3_song_mat(1:nan_jin_idx(1)-1),yp_est_no3_in,'r','linew',2);
xlabel('gure (mmol N /m^3)','fontsize',13)
ylabel('jinwall (mmol N /m^3)','fontsize',13)
title('gure vs. jinwall compare monthly interp. yymm OBS no3','fontsize',13)
grid on;
set(gca,'fontsize',13)


% nh4
figure; hold on;
plot(monthly_nh4_jin,'*','color','b','linew',1)
plot(monthly_nh4_song_mat,'o','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('nh4 (mmol N /m^3)','fontsize',13)
title('jinwall vs. gure monthly yymm OBS nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([6 17])
xlim([1 length(monthly_nh4_jin)])
legend('jinwall','gure');
set(gca,'xtick',[1:12:length(monthly_no3_jin)]);
set(gca,'xticklabel',2007:2019,'fontsize',10);
xlim([1 length(monthly_no3_jin)])

nan_jin_idx = find(isnan(monthly_nh4_jin) == 1)
nan_song_idx = find(isnan(monthly_nh4_song_mat)==1)

clearvars yp_est_nh4 nonan_idx
nonan_idx = find(isnan(monthly_nh4_song_mat + monthly_nh4_jin) == 0);
figure; hold on;
plot(monthly_nh4_song_mat,monthly_nh4_jin,'.','linew',2)
xp = monthly_nh4_song_mat(nonan_idx);
pf_nh4 = polyfit(monthly_nh4_song_mat(nonan_idx),monthly_nh4_jin(nonan_idx),1);
yp_est_nh4 = polyval(pf_nh4,monthly_nh4_song_mat(nonan_idx));
plot(monthly_nh4_song_mat(nonan_idx),yp_est_nh4,'r','linew',2);
xlabel('gure (mmol N /m^3)','fontsize',13)
ylabel('jinwall (mmol N /m^3)','fontsize',13)
title('gure vs. jinwall compare nh4','fontsize',13)
grid on;
set(gca,'fontsize',13)
% plot(mon_in_nh4_song_mat(1:nan_jin_idx(1)-1),yp_est_nh4_in,'g','linew',2);
% % there are no significant difference

% interped nh4
figure; hold on;
plot(mon_in_nh4_jin,'-','color','b','linew',1)
plot(mon_in_nh4_song_mat,'-','color','r','linew',1)
xlabel('time(days)','fontsize',13)
ylabel('nh4 (mmol N /m^3)','fontsize',13)
title('jinwall vs. gure monthly interp. yymm OBS nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([6 17])
xlim([1 length(monthly_nh4_jin)])
legend('jinwall','gure');
set(gca,'xtick',[1:12:length(monthly_no3_jin)]);
set(gca,'xticklabel',2007:2019,'fontsize',10);
xlim([1 length(monthly_no3_jin)])

clearvars yp_est_nh4_in nan_jin_idx
nan_jin_idx = find(isnan(mon_in_nh4_jin) == 1)
find(isnan(mon_in_nh4_song_mat)==1)
figure; hold on;
plot(mon_in_nh4_song_mat(1:nan_jin_idx(1)-1),mon_in_nh4_jin(1:nan_jin_idx(1)-1),'.','linew',2)
xp = mon_in_nh4_song_mat(1:nan_jin_idx(1)-1);
pf_nh4_in = polyfit(mon_in_nh4_song_mat(1:nan_jin_idx(1)-1),mon_in_nh4_jin(1:nan_jin_idx(1)-1),1);
yp_est_nh4_in = polyval(pf_nh4_in,mon_in_nh4_song_mat(1:nan_jin_idx(1)-1));
plot(mon_in_nh4_song_mat(1:nan_jin_idx(1)-1),yp_est_nh4_in,'r','linew',2);
xlabel('gure (mmol N /m^3)','fontsize',13)
ylabel('jinwall (mmol N /m^3)','fontsize',13)
title('gure vs. jinwall compare monthly interp. yymm OBS nh4','fontsize',13)
grid on;
set(gca,'fontsize',13)


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

save('jinwall_and_gure_regression_from_monthly_3sig.mat','yp_est_*','pf_*');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reconst.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %monthly mean and filled missing value
% mon_in_do_song = mon_in_do *0.7*44.661;
% mon_in_chl_song = mon_in_chl;
% mon_in_no3_song = mon_in_no3 .* 1000 ./14;
% mon_in_nh4_song = mon_in_nh4 .* 1000 ./14;

% %monthly mean and filled missing value
% mon_in_do_jin = mon_in_do *0.7*44.661;
% mon_in_chl_jin = mon_in_chl;
% mon_in_no3_jin = mon_in_no3 .* 1000 ./14;
% mon_in_nh4_jin = mon_in_nh4 .* 1000 ./14;


figure; hold on;
yp_recon_do = polyval(pf_do,mon_in_do_song);
plot(mon_in_do_song,'b','linew',1);
plot(length(monthly_do_song)-length(monthly_do_jin)+1:length(monthly_do_song)-length(monthly_do_jin)+length(mon_in_do_jin),mon_in_do_jin,'g','linew',1);
plot(yp_recon_do,'r','linew',1);
ylabel('DO (mmol/m^3)')
xlabel('year')
title('DO reconst.from monthly mean')
legend('songjung','jinwal','recon');
xlim([1 length(mon_in_do_song)])
set(gca,'xtick',[1:36:length(mon_in_do_song)]);
set(gca,'xticklabel',1989:3:2019,'fontsize',10);
grid on; set(gca,'fontsize',13)


figure; hold on;
yp_recon_chl = polyval(pf_chl,mon_in_chl_song);
plot(mon_in_chl_song,'b','linew',1);
plot(length(monthly_chl_song)-length(monthly_chl_jin)+1:length(monthly_chl_song)-length(monthly_chl_jin)+length(mon_in_chl_jin),mon_in_chl_jin,'g','linew',1);
plot(yp_recon_chl,'r','linew',1);
ylabel('chl (ug/L)')
xlabel('year')
title('Chl reconst.from monthly mean')
legend('songjung','jinwal','recon');
xlim([36*4 length(mon_in_do_song)])
% xlim([1 length(mon_in_do_song)])
set(gca,'xtick',[1:36:length(mon_in_do_song)]);
set(gca,'xticklabel',1989:3:2019,'fontsize',10);
grid on; set(gca,'fontsize',13)

% plot(length(monthly_chl_song)-length(monthly_chl_jin)+1:length(monthly_chl_song),monthly_chl_jin,'g','linew',2);

figure; hold on;
yp_recon_no3 = polyval(pf_no3,mon_in_no3_song);
plot(mon_in_no3_song,'b','linew',1);
plot(length(monthly_no3_song)-length(monthly_no3_jin)+1:length(monthly_no3_song)-length(monthly_no3_jin)+length(mon_in_no3_jin),mon_in_no3_jin,'g','linew',1);
plot(yp_recon_no3,'r','linew',1);
ylabel('NO3 (mmol N /m^3)')
xlabel('year')
title('NO3 reconst.from monthly mean')
legend('songjung','jinwal','recon');
xlim([36*2 length(mon_in_do_song)])
set(gca,'xtick',[1:36:length(mon_in_do_song)]);
set(gca,'xticklabel',1989:3:2019,'fontsize',10);
grid on; set(gca,'fontsize',13)

figure; hold on;
yp_recon_nh4 = polyval(pf_nh4,mon_in_nh4_song);
plot(mon_in_nh4_song,'b','linew',1);
plot(length(monthly_do_song)-length(monthly_do_jin)+1:length(monthly_do_song)-length(monthly_do_jin)+length(mon_in_nh4_jin),mon_in_nh4_jin,'g','linew',1);
plot(yp_recon_nh4,'r','linew',1);
ylabel('NH4 (mmol N /m^3)')
xlabel('year')
title('NH4 reconst.from monthly mean')
legend('songjung','jinwal','recon');
xlim([36*2 length(mon_in_do_song)])
set(gca,'xtick',[1:36:length(mon_in_do_song)]);
set(gca,'xticklabel',1989:3:2019,'fontsize',10);
grid on; set(gca,'fontsize',13)




%%

load('D:\장기생태\Dynamic\06_river\환경과학원\sumjin(jinwall)_polynomial_climate_advanced(v4)_3sig.mat','yp_w_*_04','reg_clim_+', 'regm_*')
yp_adv_do_jinwal = yp_w_do_04;
yp_adv_chl_jinwal = yp_w_chl_04;
yp_adv_no3_jinwal = yp_w_no3_04;
yp_adv_nh4_jinwal = yp_w_nh4_04;

clearvars yp_w_*_04
load('D:\장기생태\Dynamic\06_river\sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','yp_w_*_04'); 
yp_adv_do = yp_w_do_04;
yp_adv_chl = yp_w_chl_04;
yp_adv_no3 = yp_w_no3_04;
yp_adv_nh4 = yp_w_nh4_04;

clearvars yp_w_*_04 yp_w_*_af 
load('D:\장기생태\Dynamic\06_river\환경과학원\sumjin(songjung)_polynomial_climate_to2004_3sig.mat','yp_w_*'); % ~2004 : 1st regime, 2004~2018 : 2nd regime
clearvars yp_w_do yp_w_chl yp_w_no3 yp_w_nh4
load('D:\장기생태\Dynamic\06_river\환경과학원\sumjin(songjung)_polynomial_climate_to2004_v2_3sig.mat','yp_w_*'); % ~2018 : whole regime

% yp_adv_*_jinwal %jinwall (07~19)
% yp_adv_* % songjung whole(adv)
% yp_w_*_04 % songjung 1st regime only
% yp_w_*_af % songjung 2nd regime only
% yp_w_* % songjung only whole

figure; hold on;
yp_recon_do = polyval(pf_do,yp_adv_do);
plot(yp_recon_do,'r','linew',2);
plot(yp_adv_do,'b','linew',2);
plot(yp_adv_do_jinwal,'g','linew',2);
ylabel('mmol/m^3')
xlabel('days')
title('DO reconst.')
legend('recon','adv-do','adv-do-jinwal');
xlim([1 365])
grid on

figure; hold on;
yp_recon_chl = polyval(pf_chl,yp_adv_chl);
plot(yp_recon_chl,'r','linew',2);
plot(yp_adv_chl,'b*','linew',2);
plot(yp_adv_chl_jinwal,'g','linew',2);

% plot(length(monthly_chl_song)-length(monthly_chl_jin)+1:length(monthly_chl_song),monthly_chl_jin,'g','linew',2);

figure; hold on;
yp_recon_no3 = polyval(pf_no3,yp_adv_no3);
plot(yp_recon_no3,'r','linew',2);
plot(yp_adv_no3,'b*','linew',2);
plot(yp_adv_no3_jinwal,'g','linew',2);

figure; hold on;
yp_recon_nh4 = polyval(pf_nh4,yp_adv_nh4);
plot(yp_recon_nh4,'r','linew',2);
plot(yp_adv_nh4,'b*','linew',2);
plot(yp_adv_nh4_jinwal,'g','linew',2);

%%%%

figure; hold on;
yp_recon_do = polyval(pf_do,monthly_do_song(~isnan(monthly_do_song)));
plot(find(isnan(monthly_do_song)==0),yp_recon_do,'r','linew',2);
plot(monthly_do_song,'b*','linew',2);
plot(length(monthly_do_song)-length(monthly_do_jin)+1:length(monthly_do_song),monthly_do_jin,'g','linew',2);

figure; hold on;
yp_recon_chl = polyval(pf_chl,monthly_chl_song(~isnan(monthly_chl_song)));
plot(find(isnan(monthly_chl_song)==0),yp_recon_chl,'r','linew',2);
plot(monthly_chl_song,'b*','linew',2);
plot(length(monthly_chl_song)-length(monthly_chl_jin)+1:length(monthly_chl_song),monthly_chl_jin,'g','linew',2);

figure; hold on;
yp_recon_no3 = polyval(pf_no3,monthly_no3_song(~isnan(monthly_no3_song)));
plot(find(isnan(monthly_no3_song)==0),yp_recon_no3,'r','linew',2);
plot(monthly_no3_song,'b*','linew',2);
plot(length(monthly_no3_song)-length(monthly_no3_jin)+1:length(monthly_no3_song),monthly_nh4_jin,'g','linew',2);

figure; hold on;
yp_recon_nh4 = polyval(pf_nh4,monthly_nh4_song(~isnan(monthly_nh4_song)));
plot(find(isnan(monthly_nh4_song)==0),yp_recon_nh4,'r','linew',2);
plot(monthly_nh4_song,'b*','linew',2);
plot(length(monthly_nh4_song)-length(monthly_nh4_jin)+1:length(monthly_nh4_song),monthly_nh4_jin,'g','linew',2);


for i=1:12
    do_mon_clim_j(i) = nanmean(monthly_do_jin(i:12:end));
    chl_mon_clim_j(i) = nanmean(monthly_chl_jin(i:12:end));
    nh4_mon_clim_j(i) = nanmean(monthly_nh4_jin(i:12:end));
    no3_mon_clim_j(i) = nanmean(monthly_no3_jin(i:12:end));
    do_mon_clim_s(i) = nanmean(monthly_do_song(i:12:end));
    chl_mon_clim_s(i) = nanmean(monthly_chl_song(i:12:end));
    nh4_mon_clim_s(i) = nanmean(monthly_nh4_song(i:12:end));
    no3_mon_clim_s(i) = nanmean(monthly_no3_song(i:12:end));
end

figure; hold on;
plot(no3_mon_clim_j,'r'); plot(no3_mon_clim_s,'b');xlim([1 12])

figure; hold on;
plot(nh4_mon_clim_j,'r'); plot(nh4_mon_clim_s,'b');xlim([1 12])

%% each month regression

for i=1:12
    figure; hold on;
    plot(monthly_nh4_song_mat(i:12:end), monthly_nh4_jin(i:12:end),'*','linew',2);
    hold off
end