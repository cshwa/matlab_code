fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl,'b','linew',2);
        plot(obm_gy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
        plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl + transport']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
               
trans=ncread('..\06_river\river_2001_realts_biofennel_GY_jinwall.nc','river_transport');

yyaxis right
plot(trans(1,:),'g-','linew',2);

chl_365=gy_chl;
chl_365(60)=[];
corrcoef(trans(1,:),chl_365)

[r,lags] = xcorr(trans(1,:),chl_365,'normalized') 

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL chl vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');


figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL chl vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');