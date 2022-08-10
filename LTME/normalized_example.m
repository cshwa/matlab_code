%% nh4
fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(raw_nh4 .*1000 ./14,'b*','linew',2);
        plot(normalize(interp_nh4,'range'),'r','linew',2);
% %         plot(1:length(obm_gy_nh4),obm_gy_nh4 + obs_std_gy_nh4,'m-','linew',2);
% %         plot(1:length(obm_gy_nh4),obm_gy_nh4 - obs_std_gy_nh4,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL nh4 + transport']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_nh4'),'-dpng') 
yyaxis right
plot(normalize(sj_trans_out,'range'),'g-','linew',2);
xlim([1 length(raw_nh4)])

return

corrcoef(normalize(sj_trans_out(~isnan(interp_nh4)),'range'),normalize(interp_nh4(~isnan(interp_nh4)),'range'))
% corrcoef(sj_trans_out,interp_nh4,'Rows','complete')

figure;
plot(normalize(sj_trans_out(~isnan(interp_nh4)),'range'),normalize(interp_nh4(~isnan(interp_nh4)),'range'),'*')


clearvars r lags
[r,lags] = xcorr(sj_trans_out_x(~isnan(interp_nh4)),interp_nh4_x(~isnan(interp_nh4)),'normalized') 
% [r,lags] = xcorr(sj_trans_out(~isnan(interp_nh4)),interp_nh4(~isnan(interp_nh4)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL nh4 vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL nh4 vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');
