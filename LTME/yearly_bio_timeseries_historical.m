close all; clear all; clc;

load past_timeseriese_plot.mat

%% subplot 4,1
figPos = [0 0 5 4];
fontSizeTick = 12;
fontSizeLab  = 12;
fontSizeCLab = 12;

%arrow for event
sizeline = 2;
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos); hold on;
subplot(5,1,1);
xlim([1 5])
ylim([16.5 18]); yticks([16.5,17,17.5,18]); grid on; hold on;
xticks(1:5); 
plot(0, 0 ,'.')
plot(t_gy_all_v3(end,:)','b','linew',1); hold on;
plot(t_gy_all_v3(1,:)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
set(gca,'fontweight','bold'); box on;
ylabel('수온(^oC)','fontsize',7)

subplot(5,1,2);
xlim([1 5])
ylim([24 32]); yticks([24,26,28,30,32]); grid on; hold on;
xticks(1:5); 
plot(0, 0 ,'.')
plot(s_gy_all_v3(end,:)','b','linew',1); hold on;
plot(s_gy_all_v3(1,:)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
set(gca,'fontweight','bold'); box on;
ylabel('염분(g/kg)','fontsize',7)

subplot(5,1,3);
xlim([1 5])
ylim([250 270]); yticks([250:10:270]); grid on; hold on;
xticks(1:5); 
plot(0, 0 ,'.')
plot(o2_gy_all_v3(end,:)','b','linew',1); hold on;
plot(o2_gy_all_v3(1,:)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
set(gca,'fontweight','bold'); box on;
ylabel('용존산소(μM/L)','fontsize',7)

subplot(5,1,4);
xlim([1 5])
ylim([5 30]); yticks([5:10:30]); grid on; hold on;
xticks(1:5); 
plot(0, 0 ,'.')
plot((no3_gy_all_v3(end,:)' + nh4_gy_all_v3(end,:)'),'b','linew',1); hold on;
plot((no3_gy_all_v3(1,:)' + nh4_gy_all_v3(1,:)'),'b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
set(gca,'fontweight','bold'); box on;
ylabel('DIN(μM/L)','fontsize',7)

subplot(5,1,5);
xlim([1 5])
ylim([0 3]); yticks([0:1:3]); grid on; hold on;
xticks(1:5); 
plot(0, 0 ,'.')
plot(chl_gy_all_v3(end,:)','b','linew',1); hold on;
plot(chl_gy_all_v3(1,:)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca,'fontweight','bold'); box on;
ylabel('클로로필(μg/L)','fontsize',7)
xticklabels({'1982~1986','1987~1992','1993~2004','2007~2015','2016~2020'});
xtickangle(45);

print(gcf,['yearly_bio_timeseries_his','.png'],'-dpng','-r500');



close all; clear all; clc;

load past_timeseriese_plot.mat


%% subplot 4,1
figPos = [0 0 5 4];
fontSizeTick = 13;
fontSizeLab  = 13;
fontSizeCLab = 13;

%arrow for event
sizeline = 2;
fig = figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos); hold on;
subplot(4,1,1);
xlim([1 5])
ylim([16 22]); yticks([16,18,20,22]); grid on; hold on;
xticks(1:5); 
plot(0, 0 ,'.')
plot(t_gy_all_v3(end,:)','b','linew',1); hold on;
% plot(t_gy_all_v3(1,:)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
% set(gca,'fontweight','bold'); box on;
ylabel('수온(^oC)','fontsize',8,'fontweight','bold')
set(gca,'fontsize',9,'fontweight','bold'); box on;
ax = gca
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

% subplot(4,1,2);
% xlim([1 5])
% ylim([24 32]); yticks([24,26,28,30,32]); grid on; hold on;
% xticks(1:5); 
% plot(0, 0 ,'.')
% plot(s_gy_all_v3(end,:)','b','linew',1); hold on;
% % plot(s_gy_all_v3(1,:)','b--','linew',1);
% %set(gca,'YTick',[]); 
% set(gca, 'XTickLabel', [])
% set(gca,'fontweight','bold'); box on;
% ylabel('염분(g/kg)','fontsize',8)

subplot(4,1,2);
xlim([1 5])
ylim([240 270]); yticks([240:10:270]); grid on; hold on;
xticks(1:5); 
plot(0, 0 ,'.')
plot(o2_gy_all_v3(end,:)','b','linew',1); hold on;
% plot(o2_gy_all_v3(1,:)','b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
% set(gca,'fontweight','bold'); box on;
ylabel('용존산소(μM/L)','fontsize',8,'fontweight','bold')
set(gca,'fontsize',9,'fontweight','bold'); box on;
ax = gca
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

subplot(4,1,3);
xlim([1 5])
ylim([5 30]); yticks([5:10:30]); grid on; hold on;
xticks(1:5); 
plot(0, 0 ,'.')
plot((no3_gy_all_v3(end,:)' + nh4_gy_all_v3(end,:)'),'b','linew',1); hold on;
% plot((no3_gy_all_v3(1,:)' + nh4_gy_all_v3(1,:)'),'b--','linew',1);
%set(gca,'YTick',[]); 
set(gca, 'XTickLabel', [])
% set(gca,'fontweight','bold'); box on;
ylabel('DIN(μM/L)','fontsize',8,'fontweight','bold')
set(gca,'fontsize',9,'fontweight','bold'); box on;
ax = gca
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

subplot(4,1,4);
xlim([1 5])
ylim([1 3]); yticks([1:3]); grid on; hold on;
xticks(1:5); 
plot(0, 0 ,'.')
plot(chl_gy_all_v3(end,:)','b','linew',1); hold on;
% plot(chl_gy_all_v3(1,:)','b--','linew',1);
%set(gca,'YTick',[]); 

ylabel('클로로필(μg/L)','fontsize',8,'fontweight','bold')
xticklabels({'1982','1992','2004','2015','2020'});
xlabel('시간(년도)')
% xtickangle(45);
set(gca,'fontsize',9,'fontweight','bold'); box on;
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
% ax.XAxis.Color = 'r';

print(gcf,['yearly_bio_timeseries_his','.png'],'-dpng','-r2000');



