plt_kd = sj.yr_nh4;
corrcoef(yr_trans(~isnan(plt_kd)),plt_kd(~isnan(plt_kd)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS transp. vs. NH4-N']);
        xlabel('time','fontsize',13)
        ylabel('NH4-N (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'b*','linew',2);
plot(yr_trans,'b-','linew',2);
ylabel('transport (m^3/s)');
ylim([0 inf]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)



plt_kd = sj.yr_no3;
corrcoef(yr_trans(~isnan(plt_kd)),plt_kd(~isnan(plt_kd)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS transp. vs. NO3-N']);
        xlabel('time','fontsize',13)
        ylabel('NO3-N (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 120])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'b*','linew',2);
plot(yr_trans,'b-','linew',2);
ylabel('transport (m^3/s)');
ylim([0 inf]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)

plt_kd = sj.yr_po4;
corrcoef(yr_trans(~isnan(plt_kd)),plt_kd(~isnan(plt_kd)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(13:18,repmat(nanmean(plt_kd(1:18)),length(13:18),1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS transp. vs. PO4-P']);
        xlabel('time','fontsize',13)
        ylabel('PO4-P (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 1.1])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'b*','linew',2);
plot(yr_trans,'b-','linew',2);
ylabel('transport (m^3/s)');
ylim([0 inf]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)

%%%%%%%%%%%%%%%%
%% load vs nut
%%%%%%%%%%%%%%%%
% load nut. 
clearvars plt_nh4 plt_no3 plt_po4
plt_nh4=sj.yr_nh4 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 
plt_no3=sj.yr_no3 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 
plt_po4=sj.yr_po4 ./1000 .*30.973762.* yr_trans .* 31536000 ./ 10^6
fig = figure; hold on;
        plot(plt_no3,'ro','linew',2);
        plot(plt_nh4,'ko','linew',2);
         plot(plt_no3,'r','linew',2);
        plot(plt_nh4,'k','linew',2);
        title(['NO3-N, NH4-N, PO4-P load songjung']);
        xlabel('time','fontsize',13)
        ylabel('load (ton/yr)','fontsize',13)
        legend('NO3-N','NH4-N');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([7 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_po4,'bo','linew',2);
plot(plt_po4,'b','linew',2);
ylabel('PO4-P load (ton/yr)','fontsize',13)
ylim([0 inf])
xtickangle(45)
print(fig,strcat(['load_nut_river_yearly']),'-dpng')


% corrcoef(plt_nh4(~isnan(plt_nh4)),sj.yr_nh4(~isnan(plt_nh4)))
corrcoef(plt_no3(~isnan(plt_no3)),sj.yr_no3(~isnan(plt_no3)))
corrcoef(plt_no3(~isnan(plt_no3)),yr_trans(~isnan(plt_no3)))
fig = figure; hold on;
        plot(plt_no3,'ro','linew',2);
%         plot(plt_nh4,'ko','linew',2);
        plot(plt_no3,'r','linew',2);
%         plot(plt_nh4,'k','linew',2);
        title(['NO3-N, NH4-N load songjung vs nutrients concen.']);
        xlabel('time','fontsize',13)
        ylabel('load (ton/yr)','fontsize',13)
        legend('NO3-N load');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([7 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);        
yyaxis right
plot(yr_trans,'bo','linew',2);
% plot(sj.yr_nh4,'bo','linew',2);
plot(yr_trans,'b','linew',2);
% plot(sj.yr_nh4,'b','linew',2);
ylabel('transport (m^3/s)','fontsize',13)
% ylim([0 inf])
xtickangle(45)        


% corrcoef(plt_nh4(~isnan(plt_nh4)),sj.yr_nh4(~isnan(plt_nh4)))
corrcoef(plt_nh4(~isnan(plt_nh4)),sj.yr_nh4(~isnan(plt_nh4)))
corrcoef(plt_nh4(~isnan(plt_nh4)),yr_trans(~isnan(plt_nh4)))
fig = figure; hold on;
        plot(plt_nh4,'ro','linew',2);
%         plot(plt_nh4,'ko','linew',2);
        plot(plt_nh4,'r','linew',2);
%         plot(plt_nh4,'k','linew',2);
        title(['NH4-N load songjung vs nutrients concen.']);
        xlabel('time','fontsize',13)
        ylabel('load (ton/yr)','fontsize',13)
        legend('NH4-N load');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([7 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);        
yyaxis right
plot(yr_trans,'bo','linew',2);
% plot(sj.yr_nh4,'bo','linew',2);
plot(yr_trans,'b','linew',2);
% plot(sj.yr_nh4,'b','linew',2);
ylabel('transport (m^3/s)','fontsize',13)
% ylim([0 inf])
xtickangle(45)        


% corrcoef(plt_nh4(~isnan(plt_nh4)),sj.yr_nh4(~isnan(plt_nh4)))
corrcoef(plt_po4(~isnan(plt_po4)),sj.yr_po4(~isnan(plt_po4)))
corrcoef(plt_po4(~isnan(plt_po4)),yr_trans(~isnan(plt_po4)))
fig = figure; hold on;
        plot(plt_po4,'ro','linew',2);
%         plot(plt_po4,'ko','linew',2);
        plot(plt_po4,'r','linew',2);
%         plot(plt_po4,'k','linew',2);
        title(['PO4-P load songjung vs nutrients concen.']);
        xlabel('time','fontsize',13)
        ylabel('load (ton/yr)','fontsize',13)
        legend('PO4-P load');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([7 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);        
yyaxis right
plot(yr_trans,'bo','linew',2);
% plot(sj.yr_po4,'bo','linew',2);
plot(yr_trans,'b','linew',2);
% plot(sj.yr_po4,'b','linew',2);
ylabel('transport (m^3/s)','fontsize',13)
% ylim([0 inf])
xtickangle(45)      

corrcoef(plt_nh4(~isnan(plt_nh4)),sj.yr_nh4(~isnan(plt_nh4)))
% corrcoef(plt_nh4(~isnan(plt_nh4)),yr_trans(~isnan(plt_nh4)))
fig = figure; hold on;
        plot(plt_nh4,'ro','linew',2);
%         plot(plt_nh4,'ko','linew',2);
        plot(plt_nh4,'r','linew',2);
%         plot(plt_nh4,'k','linew',2);
        title(['NH4-N load songjung vs nutrients concen.']);
        xlabel('time','fontsize',13)
        ylabel('load (ton/yr)','fontsize',13)
        legend('NH4-N load');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([7 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);        
yyaxis right
% plot(yr_trans,'bo','linew',2);
plot(sj.yr_nh4,'bo','linew',2);
% plot(yr_trans,'b','linew',2);
plot(sj.yr_nh4,'b','linew',2);
ylabel('NH4-N (mmol/m^3)','fontsize',13)
% ylim([0 inf])
xtickangle(45)   



plt_kd = sj.yr_nh4;
corrcoef(yr_trans(~isnan(plt_kd)),plt_kd(~isnan(plt_kd)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS transp. vs. NH4-N']);
        xlabel('time','fontsize',13)
        ylabel('NH4-N (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'b*','linew',2);
plot(yr_trans,'b-','linew',2);
ylabel('transport (m^3/s)');
ylim([0 inf]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)



plt_kd = sj.yr_no3;
corrcoef(yr_trans(~isnan(plt_kd)),plt_kd(~isnan(plt_kd)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS transp. vs. NO3-N']);
        xlabel('time','fontsize',13)
        ylabel('NO3-N (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 120])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'b*','linew',2);
plot(yr_trans,'b-','linew',2);
ylabel('transport (m^3/s)');
ylim([0 inf]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)

plt_kd = sj.yr_po4;
corrcoef(yr_trans(~isnan(plt_kd)),plt_kd(~isnan(plt_kd)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(13:18,repmat(nanmean(plt_kd(1:18)),length(13:18),1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS transp. vs. PO4-P']);
        xlabel('time','fontsize',13)
        ylabel('PO4-P (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 1.1])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'b*','linew',2);
plot(yr_trans,'b-','linew',2);
ylabel('transport (m^3/s)');
ylim([0 inf]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)