close all; clear; clc; 
sig_d=load(['koem_monthly_gyonly_1to6sig_v2_each_st_16p_to06to15.mat']); % 16 points 
load koem_monthly_gyonly_nosig_v2_each_st_16p_to06to15.mat

P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;

sigma = 4;

figure; hold on;
plot([nanmean(mtx_regime_chl_1_s,1) nanmean(mtx_regime_chl_2_s,1) nanmean(mtx_regime_chl_3_s,1)],'.');
plot([nanmean(sig_d.mtx_regime_chl_1_s(:,:,sigma),1) nanmean(sig_d.mtx_regime_chl_2_s(:,:,sigma),1) nanmean(sig_d.mtx_regime_chl_3_s(:,:,sigma),1)],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('chla. (ug/L)')

figure; hold on;
plot([nanmean(mtx_regime_no3_1_s,1)./N_MW nanmean(mtx_regime_no3_2_s,1)./N_MW nanmean(mtx_regime_no3_3_s,1)./N_MW],'.');
plot([nanmean(sig_d.mtx_regime_no3_1_s(:,:,sigma),1)./N_MW nanmean(sig_d.mtx_regime_no3_2_s(:,:,sigma),1)./N_MW nanmean(sig_d.mtx_regime_no3_3_s(:,:,sigma),1)./N_MW],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('NO3 (mmol N /M^3)');

figure; hold on;
plot([nanmean(mtx_regime_nh4_1_s,1)./N_MW nanmean(mtx_regime_nh4_2_s,1)./N_MW nanmean(mtx_regime_nh4_3_s,1)./N_MW],'.');
plot([nanmean(sig_d.mtx_regime_nh4_1_s(:,:,sigma),1)./N_MW nanmean(sig_d.mtx_regime_nh4_2_s(:,:,sigma),1)./N_MW nanmean(sig_d.mtx_regime_nh4_3_s(:,:,sigma),1)./N_MW],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('NH4 (mmol N /M^3)');

figure; hold on;
plot([nanmean(mtx_regime_din_1_s,1)./N_MW nanmean(mtx_regime_din_2_s,1)./N_MW nanmean(mtx_regime_din_3_s,1)./N_MW],'.');
plot([nanmean(sig_d.mtx_regime_din_1_s(:,:,sigma),1)./N_MW nanmean(sig_d.mtx_regime_din_2_s(:,:,sigma),1)./N_MW nanmean(sig_d.mtx_regime_din_3_s(:,:,sigma),1)./N_MW],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('DIN (mmol N /M^3)');


% station plot
for sigma = 1:6
cd(['D:\장기생태\Dynamic\KOEM\to06to15_3regime\',num2str(sigma),'sigma'])
st = [1:6,14:16];
for i=1:length(st)
%   for i=4:4
fig=figure; hold on;
plot([squeeze(mtx_regime_chl_1_s(st(i),:)) squeeze(mtx_regime_chl_2_s(st(i),:)) squeeze(mtx_regime_chl_3_s(st(i),:))],'.');
plot([squeeze(sig_d.mtx_regime_chl_1_s(st(i),:,sigma)) squeeze(sig_d.mtx_regime_chl_2_s(st(i),:,sigma)) squeeze(sig_d.mtx_regime_chl_3_s(st(i),:,sigma))],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('chla. (ug/L)');
ylim([0 40]);
saveas(fig,['chla_',num2str(i),'st.png']);

fig=figure; hold on;
plot([squeeze(mtx_regime_no3_1_s(st(i),:))./N_MW squeeze(mtx_regime_no3_2_s(st(i),:))./N_MW squeeze(mtx_regime_no3_3_s(st(i),:))./N_MW],'.');
plot([squeeze(sig_d.mtx_regime_no3_1_s(st(i),:,sigma))./N_MW squeeze(sig_d.mtx_regime_no3_2_s(st(i),:,sigma))./N_MW squeeze(sig_d.mtx_regime_no3_3_s(st(i),:,sigma))./N_MW],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('NO3 (mmol N /M^3)');
ylim([0 100]);
saveas(fig,['no3_',num2str(i),'st.png']);

fig=figure; hold on;
plot([squeeze(mtx_regime_nh4_1_s(st(i),:))./N_MW squeeze(mtx_regime_nh4_2_s(st(i),:))./N_MW squeeze(mtx_regime_nh4_3_s(st(i),:))./N_MW],'.');
plot([squeeze(sig_d.mtx_regime_nh4_1_s(st(i),:,sigma))./N_MW squeeze(sig_d.mtx_regime_nh4_2_s(st(i),:,sigma))./N_MW squeeze(sig_d.mtx_regime_nh4_3_s(st(i),:,sigma))./N_MW],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('NH4 (mmol N /M^3)');
ylim([0 12]);
saveas(fig,['nh4_',num2str(i),'st.png']);

fig=figure; hold on;
plot([squeeze(mtx_regime_din_1_s(st(i),:))./N_MW squeeze(mtx_regime_din_2_s(st(i),:))./N_MW squeeze(mtx_regime_din_3_s(st(i),:))./N_MW],'.');
plot([squeeze(sig_d.mtx_regime_din_1_s(st(i),:,sigma))./N_MW squeeze(sig_d.mtx_regime_din_2_s(st(i),:,sigma))./N_MW squeeze(sig_d.mtx_regime_din_3_s(st(i),:,sigma))./N_MW],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('DIN (mmol N /M^3)');
ylim([0 120]);
saveas(fig,['din_',num2str(i),'st.png']);
end
close all;
end


% each_month_plot
tmon=[2,5,8,11];
for sigma = 1:6
cd(['D:\장기생태\Dynamic\KOEM\to06to15_3regime\each_month_time\',num2str(sigma),'sigma']);
st = [1:6,14:16];
for i=1:length(st)
    for j = 1:length(tmon)        
%   for i=4:4
fig=figure; hold on;
plot([squeeze(mtx_regime_chl_1_s(st(i),tmon(j):12:end)) squeeze(mtx_regime_chl_2_s(st(i),tmon(j):12:end)) squeeze(mtx_regime_chl_3_s(st(i),tmon(j):12:end))],'o-');
plot([squeeze(sig_d.mtx_regime_chl_1_s(st(i),tmon(j):12:end,sigma)) squeeze(sig_d.mtx_regime_chl_2_s(st(i),tmon(j):12:end,sigma)) squeeze(sig_d.mtx_regime_chl_3_s(st(i),tmon(j):12:end,sigma))],'ro-');
xticks(1:length(1997:2020)); xticklabels(1997:2020); grid on; xlim([1 length(1997:2020)]); xtickangle(45); ylabel('chla. (ug/L)');
ylim([0 40]);
saveas(fig,['chla_',num2str(j),'mon_',num2str(i),'st.png']);

fig=figure; hold on;
plot([squeeze(mtx_regime_no3_1_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_no3_2_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_no3_3_s(st(i),tmon(j):12:end))./N_MW],'o-');
plot([squeeze(sig_d.mtx_regime_no3_1_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_no3_2_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_no3_3_s(st(i),tmon(j):12:end,sigma))./N_MW],'ro-');
xticks(1:length(1997:2020)); xticklabels(1997:2020); grid on; xlim([1 length(1997:2020)]); xtickangle(45); ylabel('NO3 (mmol N /M^3)');
ylim([0 100]);
saveas(fig,['no3_',num2str(j),'mon_',num2str(i),'st.png']);

fig=figure; hold on;
plot([squeeze(mtx_regime_nh4_1_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_nh4_2_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_nh4_3_s(st(i),tmon(j):12:end))./N_MW],'o-');
plot([squeeze(sig_d.mtx_regime_nh4_1_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_nh4_2_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_nh4_3_s(st(i),tmon(j):12:end,sigma))./N_MW],'ro-');
xticks(1:length(1997:2020)); xticklabels(1997:2020); grid on; xlim([1 length(1997:2020)]); xtickangle(45); ylabel('NH4 (mmol N /M^3)');
ylim([0 12]);
saveas(fig,['nh4_',num2str(j),'mon_',num2str(i),'st.png']);

fig=figure; hold on;
plot([squeeze(mtx_regime_din_1_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_din_2_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_din_3_s(st(i),tmon(j):12:end))./N_MW],'o-');
plot([squeeze(sig_d.mtx_regime_din_1_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_din_2_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_din_3_s(st(i),tmon(j):12:end,sigma))./N_MW],'ro-');
xticks(1:length(1997:2020)); xticklabels(1997:2020); grid on; xlim([1 length(1997:2020)]); xtickangle(45); ylabel('DIN (mmol N /M^3)');
ylim([0 120]);
saveas(fig,['din_',num2str(j),'mon_',num2str(i),'st.png']);
close all
    end
end
% close all;
end

%% no cut time series confirm

close all; clear; clc; 

cd  D:\장기생태\Dynamic\06_river

sig_d=load(['koem_monthly_gyonly_1to6sig_v2_each_st_16p_no_regime_cut.mat']); % 16 points 
load koem_monthly_gyonly_nosig_v2_each_st_16p_no_regime_cut.mat

P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;

sigma = 1;
sp_st = [1:6,14:16];

figure; hold on;
plot([nanmean(mtx_regime_chl_s(sp_st,:),1)],'.');
plot([nanmean(sig_d.mtx_regime_chl_s(sp_st,:,sigma),1) ],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('chla. (ug/L)')

figure; hold on;
plot([nanmean(mtx_regime_no3_s(sp_st,:),1)./N_MW],'.');
plot([nanmean(sig_d.mtx_regime_no3_s(sp_st,:,sigma),1)./N_MW],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('NO3 (mmol N /M^3)');

figure; hold on;
plot([nanmean(mtx_regime_nh4_s(sp_st,:),1)./N_MW],'.');
plot([nanmean(sig_d.mtx_regime_nh4_s(sp_st,:,sigma),1)./N_MW],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('NH4 (mmol N /M^3)');

figure; hold on;
plot([nanmean(mtx_regime_din_s(sp_st,:),1)./N_MW],'.');
plot([nanmean(sig_d.mtx_regime_din_s(sp_st,:,sigma),1)./N_MW],'r.');
xticks(1:12:277); xticklabels(1997:2020); grid on; xlim([1 276]); xtickangle(45); ylabel('DIN (mmol N /M^3)');


% each_month_plot
tmon=[2,5,8,11];
for sigma = 1:6
cd(['D:\장기생태\Dynamic\KOEM\to06to15_3regime\each_month_time\',num2str(sigma),'sigma']);
st = [1:6,14:16];
for i=1:length(st)
    for j = 1:length(tmon)        
%   for i=4:4
fig=figure; hold on;
plot([squeeze(mtx_regime_chl_1_s(st(i),tmon(j):12:end)) squeeze(mtx_regime_chl_2_s(st(i),tmon(j):12:end)) squeeze(mtx_regime_chl_3_s(st(i),tmon(j):12:end))],'o-');
plot([squeeze(sig_d.mtx_regime_chl_1_s(st(i),tmon(j):12:end,sigma)) squeeze(sig_d.mtx_regime_chl_2_s(st(i),tmon(j):12:end,sigma)) squeeze(sig_d.mtx_regime_chl_3_s(st(i),tmon(j):12:end,sigma))],'ro-');
xticks(1:length(1997:2020)); xticklabels(1997:2020); grid on; xlim([1 length(1997:2020)]); xtickangle(45); ylabel('chla. (ug/L)');
ylim([0 40]);
saveas(fig,['chla_',num2str(j),'mon_',num2str(i),'st.png']);

fig=figure; hold on;
plot([squeeze(mtx_regime_no3_1_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_no3_2_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_no3_3_s(st(i),tmon(j):12:end))./N_MW],'o-');
plot([squeeze(sig_d.mtx_regime_no3_1_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_no3_2_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_no3_3_s(st(i),tmon(j):12:end,sigma))./N_MW],'ro-');
xticks(1:length(1997:2020)); xticklabels(1997:2020); grid on; xlim([1 length(1997:2020)]); xtickangle(45); ylabel('NO3 (mmol N /M^3)');
ylim([0 100]);
saveas(fig,['no3_',num2str(j),'mon_',num2str(i),'st.png']);

fig=figure; hold on;
plot([squeeze(mtx_regime_nh4_1_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_nh4_2_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_nh4_3_s(st(i),tmon(j):12:end))./N_MW],'o-');
plot([squeeze(sig_d.mtx_regime_nh4_1_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_nh4_2_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_nh4_3_s(st(i),tmon(j):12:end,sigma))./N_MW],'ro-');
xticks(1:length(1997:2020)); xticklabels(1997:2020); grid on; xlim([1 length(1997:2020)]); xtickangle(45); ylabel('NH4 (mmol N /M^3)');
ylim([0 12]);
saveas(fig,['nh4_',num2str(j),'mon_',num2str(i),'st.png']);

fig=figure; hold on;
plot([squeeze(mtx_regime_din_1_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_din_2_s(st(i),tmon(j):12:end))./N_MW squeeze(mtx_regime_din_3_s(st(i),tmon(j):12:end))./N_MW],'o-');
plot([squeeze(sig_d.mtx_regime_din_1_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_din_2_s(st(i),tmon(j):12:end,sigma))./N_MW squeeze(sig_d.mtx_regime_din_3_s(st(i),tmon(j):12:end,sigma))./N_MW],'ro-');
xticks(1:length(1997:2020)); xticklabels(1997:2020); grid on; xlim([1 length(1997:2020)]); xtickangle(45); ylabel('DIN (mmol N /M^3)');
ylim([0 120]);
saveas(fig,['din_',num2str(j),'mon_',num2str(i),'st.png']);
close all
    end
end
% close all;
end


% each_monthly_extract_plot
st = [1:6,14:16];
tmon=[2,5,8,11];
c_set_pick = {'r-','k-','g-','m-','b-','c-'};
cd(['D:\장기생태\Dynamic\KOEM\to06to15_3regime\each_month_time\']);
for i=1:length(st)
    for j = 1:length(tmon)
         if isnan(nanmean(squeeze(chl_sur_clim(st(i),tmon(j):12:end)))) == 0
chl_plt_sur=squeeze(chl_sur_clim(st(i),tmon(j):12:end));

fig=figure; hold on;
plot(find(isnan(chl_plt_sur)==0),chl_plt_sur(~isnan(chl_plt_sur)),'o-');
for sigma = 1:6
    yline(nanmean(chl_plt_sur)+sigma*std(chl_plt_sur,'omitnan'),c_set_pick{sigma});
end
for sigma = 1:6
    yline(nanmean(chl_plt_sur)-sigma*std(chl_plt_sur,'omitnan'),c_set_pick{sigma});
end
xticks(1:length(1997:2019)); xticklabels(1997:2019); grid on; xlim([1 length(1997:2019)]); xtickangle(45); ylabel('chla. (ug/L)');
ylim([0 40]); legend('raw','1sigma','2sigma','3sigma','4sigma','5sigma','6sigma');
saveas(fig,['chla_each_',num2str(j),'mon_extract_',num2str(i),'st.png']);

no3_plt_sur=squeeze(no3_sur_clim(st(i),tmon(j):12:end))./N_MW;
fig=figure; hold on;
plot(find(isnan(no3_plt_sur)==0),no3_plt_sur(~isnan(no3_plt_sur)),'o-');
for sigma = 1:6
    yline(nanmean(no3_plt_sur)+sigma*std(no3_plt_sur,'omitnan'),c_set_pick{sigma});
end
for sigma = 1:6
    yline(nanmean(no3_plt_sur)-sigma*std(no3_plt_sur,'omitnan'),c_set_pick{sigma});
end
xticks(1:length(1997:2019)); xticklabels(1997:2019); grid on; xlim([1 length(1997:2019)]); xtickangle(45); ylabel('NO3 (mmol N /M^3)');
ylim([0 100]); legend('raw','1sigma','2sigma','3sigma','4sigma','5sigma','6sigma');
saveas(fig,['no3_each_',num2str(j),'mon_extract_',num2str(i),'st.png']);

nh4_plt_sur=squeeze(nh4_sur_clim(st(i),tmon(j):12:end))./N_MW;
fig=figure; hold on;
plot(find(isnan(nh4_plt_sur)==0),nh4_plt_sur(~isnan(nh4_plt_sur)),'o-');
for sigma = 1:6
    yline(nanmean(nh4_plt_sur)+sigma*std(nh4_plt_sur,'omitnan'),c_set_pick{sigma});
end
for sigma = 1:6
    yline(nanmean(nh4_plt_sur)-sigma*std(nh4_plt_sur,'omitnan'),c_set_pick{sigma});
end
xticks(1:length(1997:2019)); xticklabels(1997:2019); grid on; xlim([1 length(1997:2019)]); xtickangle(45); ylabel('NH4 (mmol N /M^3)');
ylim([0 12]); legend('raw','1sigma','2sigma','3sigma','4sigma','5sigma','6sigma');
saveas(fig,['nh4_each_',num2str(j),'mon_extract_',num2str(i),'st.png']);

din_plt_sur=squeeze(din_sur_clim(st(i),tmon(j):12:end))./N_MW;
fig=figure; hold on;
plot(find(isnan(din_plt_sur)==0),din_plt_sur(~isnan(din_plt_sur)),'o-');
for sigma = 1:6
    yline(nanmean(din_plt_sur)+sigma*std(din_plt_sur,'omitnan'),c_set_pick{sigma});
end
for sigma = 1:6
    yline(nanmean(din_plt_sur)-sigma*std(din_plt_sur,'omitnan'),c_set_pick{sigma});
end
xticks(1:length(1997:2019)); xticklabels(1997:2019); grid on; xlim([1 length(1997:2019)]); xtickangle(45); ylabel('DIN (mmol N /M^3)');
ylim([0 120]); legend('raw','1sigma','2sigma','3sigma','4sigma','5sigma','6sigma');
saveas(fig,['din_each_',num2str(j),'mon_extract_',num2str(i),'st.png']);
close all
         end
    end
end




% full_month_plot
tmon=[2,5,8,11];
% for sigma = 1:6
cd(['D:\장기생태\Dynamic\KOEM\to06to15_3regime\each_month_time\']);
st = [1:6,14:16];
c_set_pick = {'r-','k-','g-','m-','b-','c-'};
for i=1:length(st)
%     for j = 1:length(tmon)        
%   for i=4:4
chl_plt_sur=squeeze(chl_sur_clim(st(i),:));
fig=figure; hold on;
plot(find(isnan(chl_plt_sur)==0),chl_plt_sur(~isnan(chl_plt_sur)),'o-');
for sigma = 1:6
    yline(nanmean(chl_plt_sur)+sigma*std(chl_plt_sur,'omitnan'),c_set_pick{sigma});
end
for sigma = 1:6
    yline(nanmean(chl_plt_sur)-sigma*std(chl_plt_sur,'omitnan'),c_set_pick{sigma});
end
plot(2:12:length(chl_plt_sur),chl_plt_sur(2:12:end),'r.');
xticks(1:12:length(1997:2019)*12); xticklabels(1997:2019); grid on; xlim([1 length(1997:2019)*12]); xtickangle(45); ylabel('chla. (ug/L)');
ylim([0 40]); legend('raw','1sigma','2sigma','3sigma','4sigma','5sigma','6sigma');
saveas(fig,['chla_full_',num2str(i),'st.png']);

no3_plt_sur=squeeze(no3_sur_clim(st(i),:))./N_MW;
fig=figure; hold on;
plot(find(isnan(no3_plt_sur)==0),no3_plt_sur(~isnan(no3_plt_sur)),'o-');
for sigma = 1:6
    yline(nanmean(no3_plt_sur)+sigma*std(no3_plt_sur,'omitnan'),c_set_pick{sigma});
end
for sigma = 1:6
    yline(nanmean(no3_plt_sur)-sigma*std(no3_plt_sur,'omitnan'),c_set_pick{sigma});
end
xticks(1:12:length(1997:2019)*12); xticklabels(1997:2019); grid on; xlim([1 length(1997:2019)*12]); xtickangle(45); ylabel('NO3 (mmol N /M^3)');
ylim([0 100]); legend('raw','1sigma','2sigma','3sigma','4sigma','5sigma','6sigma');
saveas(fig,['no3_full_',num2str(i),'st.png']);

nh4_plt_sur=squeeze(nh4_sur_clim(st(i),:))./N_MW;
fig=figure; hold on;
plot(find(isnan(nh4_plt_sur)==0),nh4_plt_sur(~isnan(nh4_plt_sur)),'o-');
for sigma = 1:6
    yline(nanmean(nh4_plt_sur)+sigma*std(nh4_plt_sur,'omitnan'),c_set_pick{sigma});
end
for sigma = 1:6
    yline(nanmean(nh4_plt_sur)-sigma*std(nh4_plt_sur,'omitnan'),c_set_pick{sigma});
end
xticks(1:12:length(1997:2019)*12); xticklabels(1997:2019); grid on; xlim([1 length(1997:2019)*12]); xtickangle(45); ylabel('NH4 (mmol N /M^3)');
ylim([0 12]); legend('raw','1sigma','2sigma','3sigma','4sigma','5sigma','6sigma');
saveas(fig,['nh4_full_',num2str(i),'st.png']);

din_plt_sur=squeeze(din_sur_clim(st(i),:))./N_MW;
fig=figure; hold on;
plot(find(isnan(din_plt_sur)==0),din_plt_sur(~isnan(din_plt_sur)),'o-');
for sigma = 1:6
    yline(nanmean(din_plt_sur)+sigma*std(din_plt_sur,'omitnan'),c_set_pick{sigma});
end
for sigma = 1:6
    yline(nanmean(din_plt_sur)-sigma*std(din_plt_sur,'omitnan'),c_set_pick{sigma});
end
xticks(1:12:length(1997:2019)*12); xticklabels(1997:2019); grid on; xlim([1 length(1997:2019)*12]); xtickangle(45); ylabel('DIN (mmol N /M^3)');
ylim([0 120]); legend('raw','1sigma','2sigma','3sigma','4sigma','5sigma','6sigma');
saveas(fig,['din_full_',num2str(i),'st.png']);
close all
end

% spatial mean chla. time series
chl_plt_sur=squeeze(chl_sur_clim(st,2:12:end));
fig=figure; hold on;
plot(nanmean(chl_plt_sur,1),'o-');
for i = 1:length(chl_plt_sur)
% errorbar(i,nanmean(chl_plt_sur(:,i)),'r.');
errorbar(i,nanmean(chl_plt_sur(:,i)),nanstd(chl_plt_sur(:,i)),'r');
end
xticks(1:length(1997:2019)); xticklabels(1997:2019); grid on; xlim([1 length(1997:2019)]); xtickangle(45); ylabel('chla. (ug/L)');
ylim([0 21]); title('spatial mean KOEM chlorophyll')
    

