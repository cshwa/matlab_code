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

set_ff = ['D:\장기생태\Dynamic\Input\RCP_input\'];

% f_list1 = dir([set_ff, 'river_3th_reg_bio_sewer_det_nosew_mixed_daily.nc']);
% f_list2 = dir([set_ff, 'bndy_csh_RCM_clim_test66_*']);
f_list{1} = [set_ff, 'bndy_auto_NWP_1_10_test6_clim_16to20.nc'];
at_list{1} = [set_ff, 'day_clim_ERA5_to2020_3regime_Tair.nc'];
k=0
for i = [2021:10:2091]
    k = k +1;
f_list{k+1} = [set_ff, 'bndy_csh_RCM_clim_test66_',num2str(i),'to',num2str(i+9),'.nc'];
at_list{k+1} = [set_ff, 'test66_GY_',num2str(i),'to',num2str(i+9),'_Tair.nc'];
end

grd_nwp = [set_ff,'grid_gy_v11_s.nc'];
lon = ncread(grd_nwp,'lon_rho');
lat = ncread(grd_nwp,'lat_rho'); 
mask = ncread(grd_nwp,'mask_rho');
mask_s = repmat(mask(:,1), 1,20);
mask_e = repmat(mask(end,:)', 1,20);

for i = 1:length(f_list)
tempo_s = []; tempo_e=[]; 
tempo_s = mean(ncread(f_list{i},'temp_south'),[3],'omitnan');
tempo_e = mean(ncread(f_list{i},'temp_east'),[3],'omitnan');
temp_s(i) = mean(tempo_s  .* (mask_s./mask_s),[1 2],'omitnan');
temp_e(i) = mean(tempo_e  .* (mask_e./mask_e),[1 2],'omitnan');
at_temp(i) = mean(ncread(at_list{i},'Tair'),[1 2 3],'omitnan');
end

f_list{1} = 'bndy_auto_NWP_1_10_test6_clim_1982to1986.nc';
% f_list{2} = 'bndy_auto_NWP_1_10_test6_clim_16to20.nc';
f_list{2} = 'bndy_auto_NWP_1_10_test6_clim_3regime.nc';

for i = 1:length(f_list)
tempo_s = []; tempo_e=[]; 
tempo_s = mean(ncread(f_list{i},'temp_south'),[3],'omitnan');
tempo_e = mean(ncread(f_list{i},'temp_east'),[3],'omitnan');
temp_s(i) = mean(tempo_s  .* (mask_s./mask_s),[1 2],'omitnan');
temp_e(i) = mean(tempo_e  .* (mask_e./mask_e),[1 2],'omitnan');
% at_temp(i) = mean(ncread(at_list{i},'Tair'),[1 2 3],'omitnan');
end



%% subplot 4,1
figPos = [0 0 5 4];
fontSizeTick = 12;
fontSizeLab  = 12;
fontSizeCLab = 12;

figure; hold on; 
plot((temp_s + temp_e)./2)
plot(at_temp)
plot(t_gy_all_v2(end,2:10))
% plot(t_gy_all_v2(1,2:10))

sizeline = 2;
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos); hold on;
plot((temp_s + temp_e)./2,'r')
plot(at_temp,'b')
plot(t_gy_all_v2(end,2:10),'g')
xlim([1 9])
ylim([14 22]); yticks([14:2:22]); grid on; hold on;
xticks(1:9); 
% plot(0, 0 ,'.')
% plot(t_gy_all_v2(end,2:10)','b','linew',1); hold on;
% plot(t_gy_all_v2(1,2:10)','b--','linew',1);
%set(gca,'YTick',[]); 
% set(gca, 'XTickLabel', [])
% set(gca,'fontweight','bold'); box on;
ylabel('^oC','fontsize',12)
set(gca,'fontsize',12,'fontweight','bold'); box on;
ax = gca
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
xticklabels({'2020','2030','2040','2050','2060', ...
    '2070', '2080', '2090', '2100'});
xlabel('시간(년도)')
% xtickangle(45);
set(gca,'fontsize',12,'fontweight','bold'); box on;
ax = gca
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
legend({'개방경계값','기온','모델결과'},'Location','northwest','NumColumns',2);
print(gcf,['yearly_BC_and_output_compare','.png'],'-dpng','-r2000');



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