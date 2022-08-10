%% chl
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(koem_chl,'b*','linew',2);
        plot(koem_in_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2018 monthly OBS transp vs. chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([98 length(koem_chl)])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(koem_secchi,'go','linew',2);
plot(koem_in_secchi,'k','linew',2);
yline(1.9)
% for i = 6:12:360; xline(i,'color','m'); end

idx_1=find(koem_secchi < 2.20);
idx_2=find(koem_secchi >= 2.20);

Kd_c = NaN;
Kd_c(idx_1) = 1.16./((koem_secchi(idx_1)).^0.62);
Kd_c(idx_2) = exp((0.15-log(koem_secchi(idx_2)).*0.62) .* (1.68-log(koem_secchi(idx_2)))./0.89  ...
    + (-0.48-log(koem_secchi(idx_2)).*0.72).* (log(koem_secchi(idx_2))-0.79)./0.89);


Kd_c(Kd_c==0)=NaN;
t=find(isnan(Kd_c)==0);
plot(t,Kd_c(~isnan(Kd_c)),'o')
plot(t,Kd_c(~isnan(Kd_c)))