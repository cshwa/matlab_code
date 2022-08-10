close all; clear all; clc;

load rcp85_timeseriese_plot.mat

%% subplot 4,1
figPos = [0 0 5 4];
fontSizeTick = 12;
fontSizeLab  = 12;
fontSizeCLab = 12;

%arrow for event
sizeline = 2;
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos); hold on;
subplot(5,1,1);
xlim([1 9])
ylim([16 22]); yticks([16:2:22]); grid on; hold on;
xticks(1:9); 
plot(0, 0 ,'.')
plot(t_gy_all_v2(end,2:10)','b','linew',1); hold on;
plot(t_gy_all_v2(1,2:10)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
set(gca,'fontweight','bold'); box on;
ylabel('수온(^oC)','fontsize',7)

subplot(5,1,2);
xlim([1 9])
ylim([28.5 31]); yticks([28.5,29,30,31]); grid on; hold on;
xticks(1:9); 
plot(0, 0 ,'.')
plot(s_gy_all_v2(end,2:10)','b','linew',1); hold on;
plot(s_gy_all_v2(1,2:10)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
set(gca,'fontweight','bold'); box on;
ylabel('염분(g/kg)','fontsize',7)

subplot(5,1,3);
xlim([1 9])
ylim([235 265]); yticks([235:10:265]); grid on; hold on;
xticks(1:9); 
plot(0, 0 ,'.')
plot(o2_gy_all_v2(end,2:10)','b','linew',1); hold on;
plot(o2_gy_all_v2(1,2:10)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
set(gca,'fontweight','bold'); box on;
ylabel('용존산소(uM/L)','fontsize',7)

subplot(5,1,4);
xlim([1 9])
ylim([8 11.5]); yticks([8:12]); grid on; hold on;
xticks(1:9); 
plot(0, 0 ,'.')
plot((no3_gy_all_v2(end,2:10)' + nh4_gy_all_v2(end,2:10)'),'b','linew',1); hold on;
plot((no3_gy_all_v2(1,2:10)' + nh4_gy_all_v2(1,2:10)'),'b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
set(gca,'fontweight','bold'); box on;
ylabel('DIN(uM/L)','fontsize',7)

subplot(5,1,5);
xlim([1 9])
ylim([0.6 2]); yticks([0.6,1,1.5,2]); grid on; hold on;
xticks(1:9); 
plot(0, 0 ,'.')
plot(chl_gy_all_v2(end,2:10)','b','linew',1); hold on;
plot(chl_gy_all_v2(1,2:10)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca,'fontweight','bold'); box on;
ylabel('클로로필(ug/L)','fontsize',7)
xticklabels({'2016~2020','2021~2030','2031~2040','2041~2050','2051~2060', ...
    '2061~2070', '2071~2080', '2081~2090', '2091~2100'});
xtickangle(45);

print(gcf,['yearly_bio_timeseries_RCP','.png'],'-dpng','-r500');


close all; clear all; clc;
load rcp85_timeseriese_plot.mat
%% subplot 4,1
figPos = [0 0 5 4];
fontSizeTick = 12;
fontSizeLab  = 12;
fontSizeCLab = 12;

%arrow for event
sizeline = 2;
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos); hold on;
subplot(4,1,1);
xlim([1 9])
ylim([16 22]); yticks([16:2:22]); grid on; hold on;
xticks(1:9); 
plot(0, 0 ,'.')
plot(t_gy_all_v2(end,2:10)','b','linew',1); hold on;
% plot(t_gy_all_v2(1,2:10)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
% set(gca,'fontweight','bold'); box on;
ylabel('수온(^oC)','fontsize',8)
set(gca,'fontsize',9,'fontweight','bold'); box on;
ax = gca
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

% subplot(5,1,2);
% xlim([1 9])
% ylim([28.5 31]); yticks([28.5,29,30,31]); grid on; hold on;
% xticks(1:9); 
% plot(0, 0 ,'.')
% plot(s_gy_all_v2(end,2:10)','b','linew',1); hold on;
% plot(s_gy_all_v2(1,2:10)','b--','linew',1);
% %set(gca,'YTick',[]); 
% set(gca, 'XTickLabel', [])
% set(gca,'fontweight','bold'); box on;
% ylabel('염분(g/kg)','fontsize',8)

subplot(4,1,2);
xlim([1 9])
ylim([240 270]); yticks([240:10:270]); grid on; hold on;
xticks(1:9); 
plot(0, 0 ,'.')
plot(o2_gy_all_v2(end,2:10)','b','linew',1); hold on;
% plot(o2_gy_all_v2(1,2:10)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
% set(gca,'fontweight','bold'); box on;
ylabel('용존산소(uM/L)','fontsize',8)
set(gca,'fontsize',9,'fontweight','bold'); box on;
ax = gca
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

subplot(4,1,3);
xlim([1 9])
ylim([5 30]); yticks([5:10:30]); grid on; hold on;
xticks(1:9); 
plot(0, 0 ,'.')
plot((no3_gy_all_v2(end,2:10)' + nh4_gy_all_v2(end,2:10)'),'b','linew',1); hold on;
% plot((no3_gy_all_v2(1,2:10)' + nh4_gy_all_v2(1,2:10)'),'b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
% set(gca,'fontweight','bold'); box on;
ylabel('DIN(uM/L)','fontsize',8)
set(gca,'fontsize',9,'fontweight','bold'); box on;
ax = gca
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

subplot(4,1,4);
xlim([1 9])
ylim([1 3]); yticks([1:3]); grid on; hold on;
xticks(1:9); 
plot(0, 0 ,'.')
plot(chl_gy_all_v2(end,2:10)','b','linew',1); hold on;
% plot(chl_gy_all_v2(1,2:10)','b--','linew',1);
%set(gca,'YTick',[]); 
% set(gca,'fontweight','bold'); box on;
ylabel('클로로필(ug/L)','fontsize',8)

xticklabels({'2020','2030','2040','2050','2060', ...
    '2070', '2080', '2090', '2100'});
xlabel('시간(년도)')
% xtickangle(45);
set(gca,'fontsize',9,'fontweight','bold'); box on;
ax = gca
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

print(gcf,['yearly_bio_timeseries_RCP','.png'],'-dpng','-r2000');