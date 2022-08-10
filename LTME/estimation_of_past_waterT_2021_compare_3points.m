close all; clc; clear all;
cd D:\장기생태\Dynamic\06_river\환경과학원
songjung=load('sumjin_recons_water_temp_present_2021.mat');
cd D:\장기생태\Dynamic\06_river_physics\환경과학원
hadong=load('sumjin_recons_water_temp(hadong)_present_2021.mat');
agyang=load('sumjin_recons_water_temp(agyang)_present_2021.mat');

figure; hold on
% plot(temp_air,'r','linew',1.5);
% plot(songjung.tx(songjung.indx_put),songjung.r_temp,'r','linew',1.5);
plot(songjung.tx,songjung.merg_recon_w,'r','linew',1.5);
plot(hadong.tx,hadong.merg_recon_w,'g','linew',1.5);
plot(agyang.tx,agyang.merg_recon_w,'b','linew',1.5);
% scatter(songjung.tx(indx_put),songjung.r_temp,'r','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp','fontsize',13)
grid on;
xlim([1 length(songjung.temp_air)]);
set(gca,'xtick',songjung.indx_num);
set(gca,'xticklabel',1990:2021);
legend({'sumjin-water(hadong)','recon. water'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_1990~2019']); hold off;
% xlim([1466 length(songjung.temp_air)]);
ylim([-5 35]); xtickangle(45);


%% high_airT

close all; clc; clear all;
cd D:\장기생태\Dynamic\06_river\환경과학원
songjung=load('sumjin_recons_water_temp_present_2021_high_airT.mat');
cd D:\장기생태\Dynamic\06_river_physics\환경과학원
hadong=load('sumjin_recons_water_temp(hadong)_present_2021_high_airT.mat');
agyang=load('sumjin_recons_water_temp(agyang)_present_2021_high_airT.mat');
jinwal=load('sumjin_recons_water_temp(jinwal)_present_2021_high_airT.mat');

figure; hold on
% plot(temp_air,'r','linew',1.5);
% plot(songjung.tx(songjung.indx_put),songjung.r_temp,'r','linew',1.5);
plot(songjung.tx,songjung.merg_recon_w,'r','linew',1.5);
plot(hadong.tx,hadong.merg_recon_w,'g','linew',1.5);
% plot(jinwal.tx,jinwal.merg_recon_w,'y','linew',1.5);
plot(agyang.tx,agyang.merg_recon_w,'b','linew',1.5);
% scatter(songjung.tx(indx_put),songjung.r_temp,'r','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp','fontsize',13)
grid on;
xlim([1 length(songjung.temp_air)]);
set(gca,'xtick',songjung.indx_num);
set(gca,'xticklabel',1990:2021);
legend({'Gure','Hadong','agyang'})
% legend({'Gure','Hadong','Jinwal','Agyang'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_1990~2019']); hold off;
% xlim([1466 length(songjung.temp_air)]);
ylim([-5 35]); xtickangle(45);

figure; hold on
% plot(temp_air,'r','linew',1.5);
% plot(songjung.tx(songjung.indx_put),songjung.r_temp,'r','linew',1.5);
% plot(songjung.tx,songjung.merg_recon_w,'r','linew',1.5);
% plot(hadong.tx,hadong.merg_recon_w,'g','linew',1.5);
plot(jinwal.tx,jinwal.merg_recon_w,'r','linew',1.5);
plot(agyang.tx,agyang.merg_recon_w,'b','linew',1.5);
% scatter(songjung.tx(indx_put),songjung.r_temp,'r','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp','fontsize',13)
grid on;
xlim([1 length(songjung.temp_air)]);
set(gca,'xtick',songjung.indx_num);
set(gca,'xticklabel',1990:2021);
% legend({'Gure','Hadong','agyang'})
legend({'Jinwal','Agyang'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_1990~2019']); hold off;
% xlim([1466 length(songjung.temp_air)]);
ylim([-5 35]); xtickangle(45);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mixed (high_airT to make regression equation and daily mean to put the x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear all;
cd D:\장기생태\Dynamic\06_river\환경과학원
songjung_ori=load('sumjin_recons_water_temp_present_2021.mat');
songjung=load('sumjin_recons_water_temp_present_2021_high_airT.mat');
songjung_mix=load('sumjin_recons_water_temp_present_2021_high_airT_mixed_daily.mat');
cd D:\장기생태\Dynamic\06_river_physics\환경과학원
hadong_ori=load('sumjin_recons_water_temp(hadong)_present_2021.mat');
agyang_ori=load('sumjin_recons_water_temp(agyang)_present_2021.mat');
jinwal_ori=load('sumjin_recons_water_temp(jinwal)_present_2021.mat');
hadong=load('sumjin_recons_water_temp(hadong)_present_2021_high_airT.mat');
agyang=load('sumjin_recons_water_temp(agyang)_present_2021_high_airT.mat');
jinwal=load('sumjin_recons_water_temp(jinwal)_present_2021_high_airT.mat');
hadong_mix=load('sumjin_recons_water_temp(hadong)_present_2021_high_airT_mixed_daily.mat');
agyang_mix=load('sumjin_recons_water_temp(agyang)_present_2021_high_airT_mixed_daily.mat');
jinwal_mix=load('sumjin_recons_water_temp(jinwal)_present_2021_high_airT_mixed_daily.mat');

figure; hold on
% plot(temp_air,'r','linew',1.5);
% plot(songjung.tx(songjung.indx_put),songjung.r_temp,'r','linew',1.5);
% plot(songjung.tx,songjung.merg_recon_w,'r','linew',1.5);
plot(songjung_mix.tx,songjung_mix.merg_recon_w,'r','linew',1.5);
% plot(hadong.tx,hadong.merg_recon_w,'g','linew',1.5);
plot(hadong_mix.tx,hadong_mix.merg_recon_w,'g','linew',1.5);
% plot(jinwal.tx,jinwal.merg_recon_w,'y','linew',1.5);
% plot(agyang.tx,agyang.merg_recon_w,'b','linew',1.5);
plot(agyang_mix.tx,agyang_mix.merg_recon_w,'b','linew',1.5);
plot(jinwal_mix.tx,jinwal_mix.merg_recon_w,'color',[.5 .5 .5],'linew',1.5);
% scatter(songjung.tx(indx_put),songjung.r_temp,'r','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp','fontsize',13)
grid on;
xlim([1 length(songjung.temp_air)]);
set(gca,'xtick',songjung.indx_num);
set(gca,'xticklabel',1990:2021);
% legend({'Gure','Hadong','Jinwal','Agyang'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_1990~2019']); hold off;
legend({'Gure','Hadong','agyang','jinwal'},'fontsize',13,'NumColumns',4)
% xlim([1466 length(songjung.temp_air)]);
ylim([-5 35]); xtickangle(45);

figure; hold on
% plot(temp_air,'r','linew',1.5);
% plot(songjung.tx(songjung.indx_put),songjung.r_temp,'r','linew',1.5);
% plot(songjung.tx,songjung.merg_recon_w,'r','linew',1.5);
plot(songjung_mix.tx,songjung_mix.merg_recon_w,'r','linew',1.5);
% plot(hadong.tx,hadong.merg_recon_w,'g','linew',1.5);
plot(hadong_mix.tx,hadong_mix.merg_recon_w,'g','linew',1.5);
% plot(jinwal.tx,jinwal.merg_recon_w,'y','linew',1.5);
% plot(agyang.tx,agyang.merg_recon_w,'b','linew',1.5);
plot(agyang_mix.tx,agyang_mix.merg_recon_w,'b','linew',1.5);
% plot(jinwal_mix.tx,jinwal_mix.merg_recon_w,'color',[.5 .5 .5],'linew',1.5);
% scatter(songjung.tx(indx_put),songjung.r_temp,'r','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp','fontsize',13)
grid on;
xlim([1 length(songjung.temp_air)]);
set(gca,'xtick',songjung.indx_num);
set(gca,'xticklabel',1990:2021);
% legend({'Gure','Hadong','Jinwal','Agyang'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_1990~2019']); hold off;
legend({'Gure','Hadong','agyang','jinwal'},'fontsize',13,'NumColumns',4)
% xlim([1466 length(songjung.temp_air)]);
ylim([-5 35]); xtickangle(45);



%% Gure
figure; hold on
plot(songjung.tx,songjung.merg_recon_w,'r','linew',1.5);
plot(songjung_ori.tx,songjung_ori.merg_recon_w,'g','linew',1.5);
plot(songjung_mix.tx,songjung_mix.merg_recon_w,'b','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp(Gure)','fontsize',13)
grid on;
xlim([1 length(songjung.temp_air)]);
set(gca,'xtick',songjung.indx_num);
set(gca,'xticklabel',1990:2021);
% legend({'Gure','Hadong','agyang'})
legend({'Jinwal','Agyang'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_1990~2019']); hold off;
% xlim([1466 length(songjung.temp_air)]);
ylim([-5 35]); xtickangle(45);
legend({'일최고','일평균','일최고&일평균'},'fontsize',10,'NumColumns',3)

%% Hadong
figure; hold on
plot(hadong.tx,hadong.merg_recon_w,'r','linew',1.5);
plot(hadong_ori.tx,hadong_ori.merg_recon_w,'g','linew',1.5);
plot(hadong_mix.tx,hadong_mix.merg_recon_w,'b','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp(Hadong)','fontsize',13)
grid on;
xlim([1 length(hadong.temp_air)]);
set(gca,'xtick',hadong.indx_num);
set(gca,'xticklabel',1990:2021);
% legend({'Gure','Hadong','agyang'})
legend({'Jinwal','Agyang'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_1990~2019']); hold off;
% xlim([1466 length(hadong.temp_air)]);
ylim([-5 35]); xtickangle(45);
legend({'일최고','일평균','일최고&일평균'},'fontsize',10,'NumColumns',3)

%% Agyang
figure; hold on
plot(agyang.tx,agyang.merg_recon_w,'r','linew',1.5);
plot(agyang_ori.tx,agyang_ori.merg_recon_w,'g','linew',1.5);
plot(agyang_mix.tx,agyang_mix.merg_recon_w,'b','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp(Agyang)','fontsize',13)
grid on;
xlim([1 length(agyang.temp_air)]);
set(gca,'xtick',agyang.indx_num);
set(gca,'xticklabel',1990:2021);
% legend({'Gure','agyang','agyang'})
legend({'Jinwal','Agyang'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_1990~2019']); hold off;
% xlim([1466 length(agyang.temp_air)]);
ylim([-5 35]); xtickangle(45);
legend({'일최고','일평균','일최고&일평균'},'fontsize',10,'NumColumns',3)

%% Jinwal
figure; hold on
plot(jinwal.tx,jinwal.merg_recon_w,'r','linew',1.5);
plot(jinwal_ori.tx,jinwal_ori.merg_recon_w,'g','linew',1.5);
plot(jinwal_mix.tx,jinwal_mix.merg_recon_w,'b','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp(Jinwal)','fontsize',13)
grid on;
xlim([1 length(jinwal.temp_air)]);
set(gca,'xtick',jinwal.indx_num);
set(gca,'xticklabel',1990:2021);
% legend({'Gure','jinwal','jinwal'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_1990~2019']); hold off;
% xlim([1466 length(jinwal.temp_air)]);
ylim([-5 35]); xtickangle(45);
legend({'일최고','일평균','일최고&일평균'},'fontsize',10,'NumColumns',3)


