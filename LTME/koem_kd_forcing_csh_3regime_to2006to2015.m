close all; clear; clc;
% -------------------------------------------------------------------------
t_year = 1989:2018
cd D:\장기생태\Dynamic\06_river
% from plot_KOEM_data_process_step5_sig_v2_each_st_to06to15.m
yj=load('koem_timeseires_monthly_gy_only_1to3sig_v2_each_st_16p.mat'); % 16 points 

%% 9 point
sp_9p_st = [1:6, 14:16];

koem_secchi_3s=yj.regime_secchi_3s(:,1:end);
koem_secchi_sp=nanmean(koem_secchi_3s(sp_9p_st,:),1);
koem_kd_sp=nanmean(1.7./koem_secchi_3s(sp_9p_st,:),1);

%% cut 3 regime
cut_1 = 120 % 2006-12
cut_2 = 228 % 2015-12

mon_secchi_1 = koem_secchi_3s(:,1:cut_1);
mon_secchi_2 = koem_secchi_3s(:,cut_1+1:cut_2);
mon_secchi_3 = koem_secchi_3s(:,cut_2+1:end);

mon_secchi_sp_1 = koem_secchi_sp(1:cut_1);
mon_secchi_sp_2 = koem_secchi_sp(cut_1+1:cut_2);
mon_secchi_sp_3 = koem_secchi_sp(cut_2+1:end);

mon_kd_1 = 1.7./ koem_secchi_3s(:,1:cut_1);
mon_kd_2 = 1.7./ koem_secchi_3s(:,cut_1+1:cut_2);
mon_kd_3 = 1.7./ koem_secchi_3s(:,cut_2+1:end);

mon_kd_sp_1 = koem_kd_sp(1:cut_1);
mon_kd_sp_2 = koem_kd_sp(cut_1+1:cut_2);
mon_kd_sp_3 = koem_kd_sp(cut_2+1:end);

for i = 1:12 %month   
mon_clim_secchi_1(:,i) = nanmean(mon_secchi_1(:,i:12:end),2);
mon_clim_secchi_2(:,i) = nanmean(mon_secchi_2(:,i:12:end),2);
mon_clim_secchi_3(:,i) = nanmean(mon_secchi_3(:,i:12:end),2); 

mon_clim_secchi_sp_1(i) = nanmean(mon_secchi_sp_1(i:12:end));
mon_clim_secchi_sp_2(i) = nanmean(mon_secchi_sp_2(i:12:end));
mon_clim_secchi_sp_3(i) = nanmean(mon_secchi_sp_3(i:12:end)); 

mon_clim_kd_1(:,i) = nanmean(mon_kd_1(:,i:12:end),2);
mon_clim_kd_2(:,i) = nanmean(mon_kd_2(:,i:12:end),2);
mon_clim_kd_3(:,i) = nanmean(mon_kd_3(:,i:12:end),2); 

mon_clim_kd_sp_1(i) = nanmean(mon_kd_sp_1(i:12:end));
mon_clim_kd_sp_2(i) = nanmean(mon_kd_sp_2(i:12:end));
mon_clim_kd_sp_3(i) = nanmean(mon_kd_sp_3(i:12:end)); 
end

figure; hold on;
plot(mon_clim_secchi_1,'r'); plot(mon_clim_secchi_2,'g'); plot(mon_clim_secchi_3,'b'); grid on;

figure; hold on;
plot(mon_clim_secchi_sp_1,'ro'); plot(mon_clim_secchi_sp_2,'go'); plot(mon_clim_secchi_sp_3,'bo'); grid on;
disp(nanmean(mon_clim_secchi_sp_1)); disp(nanmean(mon_clim_secchi_sp_2)); disp(nanmean(mon_clim_secchi_sp_3))

figure; hold on;
plot(mon_clim_kd_sp_1,'ro'); plot(mon_clim_kd_sp_2,'go'); plot(mon_clim_kd_sp_3,'bo'); grid on;
disp(nanmean(mon_clim_kd_sp_1)); disp(nanmean(mon_clim_kd_sp_2)); disp(nanmean(mon_clim_kd_sp_3))

day_indx=[1, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 365];
day_indx_3 = [day_indx, day_indx+365, day_indx+730];

clearvars day_clim_*
day_clim_kd_1 = NaN(16,365*3);
day_clim_kd_2 = NaN(16,365*3);
day_clim_kd_3 = NaN(16,365*3);

day_clim_kd_sp_1 = NaN(365*3,1);
day_clim_kd_sp_2 = NaN(365*3,1);
day_clim_kd_sp_3 = NaN(365*3,1);

day_clim_secchi_1 = NaN(16,365*3);
day_clim_secchi_2 = NaN(16,365*3);
day_clim_secchi_3 = NaN(16,365*3);

day_clim_secchi_sp_1 = NaN(365*3,1);
day_clim_secchi_sp_2 = NaN(365*3,1);
day_clim_secchi_sp_3 = NaN(365*3,1);


%%
day_clim_kd_1(:,day_indx_3) = [mon_clim_kd_1, mon_clim_kd_1, mon_clim_kd_1]
day_clim_kd_2(:,day_indx_3) = [mon_clim_kd_2, mon_clim_kd_2, mon_clim_kd_2]
day_clim_kd_3(:,day_indx_3) = [mon_clim_kd_3, mon_clim_kd_3, mon_clim_kd_3]

day_clim_kd_sp_1(day_indx_3) = [mon_clim_kd_sp_1, mon_clim_kd_sp_1, mon_clim_kd_sp_1];
day_clim_kd_sp_2(day_indx_3) = [mon_clim_kd_sp_2, mon_clim_kd_sp_2, mon_clim_kd_sp_2];
day_clim_kd_sp_3(day_indx_3) = [mon_clim_kd_sp_3, mon_clim_kd_sp_3, mon_clim_kd_sp_3];

day_clim_secchi_1(:,day_indx_3) = [mon_clim_secchi_1, mon_clim_secchi_1, mon_clim_secchi_1];
day_clim_secchi_2(:,day_indx_3) = [mon_clim_secchi_2, mon_clim_secchi_2, mon_clim_secchi_2];
day_clim_secchi_3(:,day_indx_3) = [mon_clim_secchi_3, mon_clim_secchi_3, mon_clim_secchi_3];

day_clim_secchi_sp_1(day_indx_3) = [mon_clim_secchi_sp_1, mon_clim_secchi_sp_1, mon_clim_secchi_sp_1]
day_clim_secchi_sp_2(day_indx_3) = [mon_clim_secchi_sp_2, mon_clim_secchi_sp_2, mon_clim_secchi_sp_2]
day_clim_secchi_sp_3(day_indx_3) = [mon_clim_secchi_sp_3, mon_clim_secchi_sp_3, mon_clim_secchi_sp_3]

%%
for i = 1:9
    clearvars temp_*
    temp_kd_1 =  day_clim_kd_1(sp_9p_st(i),:)
    temp_kd_2 =  day_clim_kd_2(sp_9p_st(i),:)
    temp_kd_3 =  day_clim_kd_3(sp_9p_st(i),:)
    temp_secchi_1 =  day_clim_secchi_1(sp_9p_st(i),:)
    temp_secchi_2 =  day_clim_secchi_2(sp_9p_st(i),:)
    temp_secchi_3 =  day_clim_secchi_3(sp_9p_st(i),:)

t=1:length(temp_kd_1);
day_clim_kd_1(sp_9p_st(i),isnan(temp_kd_1)) = interp1(t(~isnan(temp_kd_1)),temp_kd_1(~isnan(temp_kd_1)),t(isnan(temp_kd_1)));
t=1:length(temp_kd_2);
day_clim_kd_2(sp_9p_st(i),isnan(temp_kd_2)) = interp1(t(~isnan(temp_kd_2)),temp_kd_2(~isnan(temp_kd_2)),t(isnan(temp_kd_2)));
t=1:length(temp_kd_3);
day_clim_kd_3(sp_9p_st(i),isnan(temp_kd_3)) = interp1(t(~isnan(temp_kd_3)),temp_kd_3(~isnan(temp_kd_3)),t(isnan(temp_kd_3)));

t=1:length(temp_secchi_1);
day_clim_secchi_1(sp_9p_st(i),isnan(temp_secchi_1)) = interp1(t(~isnan(temp_secchi_1)),temp_secchi_1(~isnan(temp_secchi_1)),t(isnan(temp_secchi_1)));
t=1:length(temp_secchi_2);
day_clim_secchi_2(sp_9p_st(i),isnan(temp_secchi_2)) = interp1(t(~isnan(temp_secchi_2)),temp_secchi_2(~isnan(temp_secchi_2)),t(isnan(temp_secchi_2)));
t=1:length(temp_secchi_3);
day_clim_secchi_3(sp_9p_st(i),isnan(temp_secchi_3)) = interp1(t(~isnan(temp_secchi_3)),temp_secchi_3(~isnan(temp_secchi_3)),t(isnan(temp_secchi_3)));
end

t=1:length(day_clim_kd_sp_1);
day_clim_kd_sp_1(isnan(day_clim_kd_sp_1)) = interp1(t(~isnan(day_clim_kd_sp_1)),day_clim_kd_sp_1(~isnan(day_clim_kd_sp_1)),t(isnan(day_clim_kd_sp_1)));
t=1:length(day_clim_kd_sp_2);
day_clim_kd_sp_2(isnan(day_clim_kd_sp_2)) = interp1(t(~isnan(day_clim_kd_sp_2)),day_clim_kd_sp_2(~isnan(day_clim_kd_sp_2)),t(isnan(day_clim_kd_sp_2)));
t=1:length(day_clim_kd_sp_3);
day_clim_kd_sp_3(isnan(day_clim_kd_sp_3)) = interp1(t(~isnan(day_clim_kd_sp_3)),day_clim_kd_sp_3(~isnan(day_clim_kd_sp_3)),t(isnan(day_clim_kd_sp_3)));

t=1:length(day_clim_secchi_sp_1);
day_clim_secchi_sp_1(isnan(day_clim_secchi_sp_1)) = interp1(t(~isnan(day_clim_secchi_sp_1)),day_clim_secchi_sp_1(~isnan(day_clim_secchi_sp_1)),t(isnan(day_clim_secchi_sp_1)));
t=1:length(day_clim_secchi_sp_2);
day_clim_secchi_sp_2(isnan(day_clim_secchi_sp_2)) = interp1(t(~isnan(day_clim_secchi_sp_2)),day_clim_secchi_sp_2(~isnan(day_clim_secchi_sp_2)),t(isnan(day_clim_secchi_sp_2)));
t=1:length(day_clim_secchi_sp_3);
day_clim_secchi_sp_3(isnan(day_clim_secchi_sp_3)) = interp1(t(~isnan(day_clim_secchi_sp_3)),day_clim_secchi_sp_3(~isnan(day_clim_secchi_sp_3)),t(isnan(day_clim_secchi_sp_3)));

day_clim_kd_sp_1_in= day_clim_kd_sp_1(366:730); day_clim_kd_sp_2_in= day_clim_kd_sp_2(366:730);
day_clim_kd_sp_3_in= day_clim_kd_sp_3(366:730);

figure; hold on;
plot(day_clim_kd_sp_1_in,'r'); plot(day_clim_kd_sp_2_in,'g'); plot(day_clim_kd_sp_3_in,'b'); grid on;
disp(nanmean(day_clim_kd_sp_1_in)); disp(nanmean(day_clim_kd_sp_2_in)); disp(nanmean(day_clim_kd_sp_3_in))
xlim([1 365]); grid on; ylabel('Kd coefficient'); xlabel('days')

figure; hold on; xlim([1 365]); grid on; ylabel('Kd coefficient'); xlabel('days')
plot(day_clim_kd_sp_1_in,'r'); grid on; plot(day_clim_kd_1(sp_9p_st,366:730)','color',[.5 .5 .5]); ylim([0 1.8])
figure; hold on; xlim([1 365]); grid on; ylabel('Kd coefficient'); xlabel('days')
plot(day_clim_kd_sp_1_in,'g'); grid on; plot(day_clim_kd_2(sp_9p_st,366:730)','color',[.5 .5 .5]); ylim([0 1.8])
figure; hold on; xlim([1 365]); grid on; ylabel('Kd coefficient'); xlabel('days')
plot(day_clim_kd_sp_1_in,'b'); grid on; plot(day_clim_kd_3(sp_9p_st,366:730)','color',[.5 .5 .5]); ylim([0 1.8])

save('Koem_kd_days_climate_to06to15.mat','day_clim_*');

%% Kd vs. CHL
plt_kd = 1.7 ./yr_secchi;
corrcoef(yr_chl(~isnan(yr_chl)),plt_kd(~isnan(yr_chl)))
fig = figure; hold on;
        plot(plt_kd,'b*','linew',2);
        plot(plt_kd,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. Kd']);
        xlabel('time','fontsize',13)
        ylabel('Kd (m^-1)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:30);
        xlim([1 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
ylim([0 10]);
xtickangle(45)

%% secchi
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_secchi,'b*','linew',2);
        plot(yr_secchi,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 yearly OBS transp vs. secchi']);
        xlabel('time','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:30);
        xlim([1 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)

% Kd
for i = 1:16
clearvars secc r_coef plt_chl coe_* Kd_c idx_*
    secc = yr_secchi_raw(i,:);
    idx_1=find(secc < 2.20);
    idx_2=find(secc >= 2.20);
    Kd_c(idx_1) = 1.16./((secc(idx_1)).^0.62);
    Kd_c(idx_2) = exp((0.15-log(secc(idx_2)).*0.62) .* (1.68-log(secc(idx_2)))./0.89  ...
    + (-0.48-log(secc(idx_2)).*0.72).* (log(secc(idx_2))-0.79)./0.89);
    Kd_c(Kd_c==0)=NaN;
    plt_chl = yr_chl_raw(i,:);
    coe_secc=Kd_c(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_secc)),coe_secc(~isnan(coe_secc)))    
fig = figure; hold on;
        plot(Kd_c,'ko','linew',2);
        plot(Kd_c,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. Kd st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        text(9,2.8, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:30);
        xlim([9 31])
        ylim([0 1.4])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_Kdc_koem_st_',char(name_tag{i})]),'-dpng')
end
