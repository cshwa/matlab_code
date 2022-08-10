close all; clear; clc;
% -------------------------------------------------------------------------
t_year = 1989:2018
cd D:\장기생태\Dynamic\06_river
% from plot_KOEM_data_process_step5_sig_v2_each_st_to06to15.m
% yj=load(['koem_monthly_gyonly_1to6sig_v2_each_st_16p_to06to15.mat']); % 16 points 
yj=load(['koem_monthly_gyonly_1to6sig_v2_each_st_16p_to06to15_2020.mat']); % 16 points 

%% 9 point
sp_9p_st = [1:6, 14:16];
sig = 3; % 3sigma

koem_secchi_1_3s=squeeze(yj.mtx_regime_secchi_1_s(:,1:end, sig));
koem_secchi_1_sp=nanmean(koem_secchi_1_3s(sp_9p_st,:),1);
koem_kd_1_sp=nanmean(1.7./koem_secchi_1_3s(sp_9p_st,:),1);

koem_secchi_2_3s=squeeze(yj.mtx_regime_secchi_2_s(:,1:end, sig));
koem_secchi_2_sp=nanmean(koem_secchi_2_3s(sp_9p_st,:),1);
koem_kd_2_sp=nanmean(1.7./koem_secchi_2_3s(sp_9p_st,:),1);

koem_secchi_3_3s=squeeze(yj.mtx_regime_secchi_3_s(:,1:end, sig));
koem_secchi_3_sp=nanmean(koem_secchi_3_3s(sp_9p_st,:),1);
koem_kd_3_sp=nanmean(1.7./koem_secchi_3_3s(sp_9p_st,:),1);


%% Make Kd_c form ref.
for i = 1:3 %regime 
clearvars temp_secchi_all kd_c idx_* kd_c_s
eval(['temp_secchi_all',' = koem_secchi_',num2str(i),'_3s(sp_9p_st,:);']);
for j = 1:length(sp_9p_st)
    clearvars temp_secchi kd_c
    temp_secchi=temp_secchi_all(j,:);
    idx_1=find(temp_secchi < 2.20);
    idx_2=find(temp_secchi >= 2.20);
    Kd_c = NaN(length(temp_secchi),1);
    Kd_c(idx_1) = 1.16./((temp_secchi(idx_1)).^0.62);
    Kd_c(idx_2) = exp((0.15-log(temp_secchi(idx_2)).*0.62) .* (1.68-log(temp_secchi(idx_2)))./0.89  ...
        + (-0.48-log(temp_secchi(idx_2)).*0.72).* (log(temp_secchi(idx_2))-0.79)./0.89);
    Kd_c(Kd_c==0)=NaN;
    eval(['kd_c_s(j,:)',' = Kd_c;']);
if j == length(sp_9p_st)
    eval(['kd_c_',num2str(i),'_s',' = nanmean(kd_c_s,1);']);
end
end
end
t=find(isnan(kd_c_1_s)==0); t2=find(isnan(kd_c_2_s)==0); t3=find(isnan(kd_c_3_s)==0);
plot(t,kd_c_1_s(~isnan(kd_c_1_s)),'o'); hold on;
plot(t,kd_c_1_s(~isnan(kd_c_1_s)))

figure; hold on; xlim([1 length(kd_c_1_s)]); grid on; ylabel('Kd coefficient'); xlabel('month')
plot(t,kd_c_1_s(~isnan(kd_c_1_s)),'r'); plot(t2,kd_c_2_s(~isnan(kd_c_2_s)),'g');  plot(t3,kd_c_3_s(~isnan(kd_c_3_s)),'b');

for i = 1:12 %month   
mon_clim_kd_c_1(i) = nanmean(kd_c_1_s(:,i:12:end),2);
mon_clim_kd_c_2(i) = nanmean(kd_c_2_s(:,i:12:end),2);
mon_clim_kd_c_3(i) = nanmean(kd_c_3_s(:,i:12:end),2);  
end

figure; hold on;
plot(mon_clim_kd_c_1,'ro'); plot(mon_clim_kd_c_2,'go'); plot(mon_clim_kd_c_3,'bo'); grid on;
disp(nanmean(mon_clim_kd_c_1)); disp(nanmean(mon_clim_kd_c_2)); disp(nanmean(mon_clim_kd_c_3))


%% cut 3 regime

mon_secchi_1 = koem_secchi_1_3s;
mon_secchi_2 = koem_secchi_2_3s;
mon_secchi_3 = koem_secchi_3_3s;

mon_secchi_sp_1 = koem_secchi_1_sp;
mon_secchi_sp_2 = koem_secchi_2_sp;
mon_secchi_sp_3 = koem_secchi_3_sp;

mon_kd_1 = 1.7./ koem_secchi_1_3s;
mon_kd_2 = 1.7./ koem_secchi_2_3s;
mon_kd_3 = 1.7./ koem_secchi_3_3s;

mon_kd_sp_1 = koem_kd_1_sp;
mon_kd_sp_2 = koem_kd_2_sp;
mon_kd_sp_3 = koem_kd_3_sp;

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

day_clim_kdc_sp_1 = NaN(365*3,1);
day_clim_kdc_sp_2 = NaN(365*3,1);
day_clim_kdc_sp_3 = NaN(365*3,1);

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

day_clim_kdc_sp_1(day_indx_3) = [mon_clim_kd_c_1, mon_clim_kd_c_1, mon_clim_kd_c_1];
day_clim_kdc_sp_2(day_indx_3) = [mon_clim_kd_c_2, mon_clim_kd_c_2, mon_clim_kd_c_2];
day_clim_kdc_sp_3(day_indx_3) = [mon_clim_kd_c_3, mon_clim_kd_c_3, mon_clim_kd_c_3];

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

t=1:length(day_clim_kdc_sp_1);
day_clim_kdc_sp_1(isnan(day_clim_kdc_sp_1)) = interp1(t(~isnan(day_clim_kdc_sp_1)),day_clim_kdc_sp_1(~isnan(day_clim_kdc_sp_1)),t(isnan(day_clim_kdc_sp_1)));
t=1:length(day_clim_kdc_sp_2);
day_clim_kdc_sp_2(isnan(day_clim_kdc_sp_2)) = interp1(t(~isnan(day_clim_kdc_sp_2)),day_clim_kdc_sp_2(~isnan(day_clim_kdc_sp_2)),t(isnan(day_clim_kdc_sp_2)));
t=1:length(day_clim_kdc_sp_3);
day_clim_kdc_sp_3(isnan(day_clim_kdc_sp_3)) = interp1(t(~isnan(day_clim_kdc_sp_3)),day_clim_kdc_sp_3(~isnan(day_clim_kdc_sp_3)),t(isnan(day_clim_kdc_sp_3)));

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

day_clim_kdc_sp_1_in= day_clim_kdc_sp_1(366:730); day_clim_kdc_sp_2_in= day_clim_kdc_sp_2(366:730);
day_clim_kdc_sp_3_in= day_clim_kdc_sp_3(366:730);

figure; hold on;
plot(day_clim_kd_sp_1_in,'r'); plot(day_clim_kd_sp_2_in,'g'); plot(day_clim_kd_sp_3_in,'b'); grid on;
disp(nanmean(day_clim_kd_sp_1_in)); disp(nanmean(day_clim_kd_sp_2_in)); disp(nanmean(day_clim_kd_sp_3_in))
xlim([1 365]); grid on; ylabel('Kd coefficient'); xlabel('days')
plot(day_clim_kdc_sp_1_in,'r--'); plot(day_clim_kdc_sp_2_in,'g--'); plot(day_clim_kdc_sp_3_in,'b--'); grid on;
disp(nanmean(day_clim_kdc_sp_1_in)); disp(nanmean(day_clim_kdc_sp_2_in)); disp(nanmean(day_clim_kdc_sp_3_in))
legend('Kd(1997~2006)','Kd(2007~2015)','Kd(2016~2020)','Kd_c(1997~2006)','Kd_c(2007~2015)','Kd_c(2016~2020)');
ylim([0.2 1.2])


figure; hold on; xlim([1 365]); grid on; ylabel('Kd coefficient'); xlabel('days')
plot(day_clim_kd_sp_1_in,'r'); grid on; plot(day_clim_kd_1(sp_9p_st,366:730)','color',[.5 .5 .5]); ylim([0 1.8])
figure; hold on; xlim([1 365]); grid on; ylabel('Kd coefficient'); xlabel('days')
plot(day_clim_kd_sp_2_in,'g'); grid on; plot(day_clim_kd_2(sp_9p_st,366:730)','color',[.5 .5 .5]); ylim([0 1.8])
figure; hold on; xlim([1 365]); grid on; ylabel('Kd coefficient'); xlabel('days')
plot(day_clim_kd_sp_3_in,'b'); grid on; plot(day_clim_kd_3(sp_9p_st,366:730)','color',[.5 .5 .5]); ylim([0 1.8])

save('Koem_kd_days_climate_to06to15_v2_2020.mat','day_clim_*');

%% Kd vs. CHL
kd_times = 4;
kd1=(mean(day_clim_kd_sp_1_in)); kd2=(mean(day_clim_kd_sp_2_in)); kd3=(mean(day_clim_kd_sp_3_in));
kdc1=(mean(day_clim_kdc_sp_1_in)); kdc2=(mean(day_clim_kdc_sp_2_in)); kdc3=(mean(day_clim_kdc_sp_3_in));
kdt2=(kd_times.* mean(day_clim_kd_sp_1_in)); kdt2_2=(kd_times.*mean(day_clim_kd_sp_2_in)); kdt2_3=(kd_times.* mean(day_clim_kd_sp_3_in));
% Kd
for k = 1:20
if k==1
PAR(k) = (1-exp(-kd1))/kd1
else
PAR(k) = PAR(k-1) * (1-exp(-kd1))/kd1
end
end

for k = 1:20
if k==1
PAR2(k) = (1-exp(-kd2))/kd2
else
PAR2(k) = PAR2(k-1) * (1-exp(-kd2))/kd2
end
end

for k = 1:20
if k==1
PAR3(k) = (1-exp(-kd3))/kd3
else
PAR3(k) = PAR3(k-1) * (1-exp(-kd3))/kd3
end
end

% Kdc
for k = 1:20
if k==1
PARC(k) = (1-exp(-kdc1))/kdc1
else
PARC(k) = PARC(k-1) * (1-exp(-kdc1))/kdc1
end
end

for k = 1:20
if k==1
PARC2(k) = (1-exp(-kdc2))/kdc2
else
PARC2(k) = PARC2(k-1) * (1-exp(-kdc2))/kdc2
end
end

for k = 1:20
if k==1
PARC3(k) = (1-exp(-kdc3))/kdc3
else
PARC3(k) = PARC3(k-1) * (1-exp(-kdc3))/kdc3
end
end

% Kd
for k = 1:20
if k==1
PARR(k) = (1-exp(-0.04))/0.04
else
PARR(k) = PARR(k-1) * (1-exp(-0.04))/0.04
end
end

% kd * 2
for k = 1:20
if k==1
PART(k) = (1-exp(-kdt2))/kdt2
else
PART(k) = PART(k-1) * (1-exp(-kdt2))/kdt2
end
end

for k = 1:20
if k==1
PART2(k) = (1-exp(-kdt2_2))/kdt2_2
else
PART2(k) = PART2(k-1) * (1-exp(-kdt2_2))/kdt2_2
end
end

for k = 1:20
if k==1
PART3(k) = (1-exp(-kdt2_3))/kdt2_3
else
PART3(k) = PART3(k-1) * (1-exp(-kdt2_3))/kdt2_3
end
end


figure; hold on;
plot([1 PARR],'k'); 
plot([1 PAR],'r');  plot([1 PAR2],'g'); plot([1 PAR3],'b'); grid on;
plot([1 PARC],'r--');  plot([1 PARC2],'g--'); plot([1 PARC3],'b--'); grid on;
legend('Kd(fennel)','Kd(1st)','Kd(2nd)','Kd(3rd)','KdC(1st)','KdC(2nd)','KdC(3rd)')
xlim([1 21]);
xlabel('num. of vertical layer'); ylabel('Normalizaed irradiance')
xticks(1:21); xticklabels(0:20);

figure; hold on;
plot([1 PARR],'k'); 
plot([1 PAR],'r');  plot([1 PAR2],'g'); plot([1 PAR3],'b'); grid on;
% plot([1 PARC],'r--');  plot([1 PARC2],'g--'); plot([1 PARC3],'b--'); 
plot([1 PART],'r--');  plot([1 PART2],'g--'); plot([1 PART3],'b--'); 
grid on;
legend('Kd(fennel)','Kd(1st)','Kd(2nd)','Kd(3rd)',[num2str(kd_times),'*Kd(1st)'],[num2str(kd_times),'*Kd(2nd)'],[num2str(kd_times),'*Kd(3rd)'])
xlim([1 21]);
xlabel('num. of vertical layer'); ylabel('Normalizaed irradiance')
xticks(1:21); xticklabels(0:20);


figure; hold on;
plot([1 PARR],'k'); 
plot([1 PAR],'r');  plot([1 PAR2],'g'); plot([1 PAR3],'b'); grid on;
% plot([1 PARC],'r--');  plot([1 PARC2],'g--'); plot([1 PARC3],'b--'); 
% plot([1 PART],'r--');  plot([1 PART2],'g--'); plot([1 PART3],'b--'); 
grid on;
legend('Kd(fennel)','Kd(1st)','Kd(2nd)','Kd(3rd)')
xlim([1 21]);
xlabel('Depth (m)'); ylabel('Normalizaed irradiance')
xticks(1:21); xticklabels(0:20);




%% + add Chlorophyll
clim1=load('D:\장기생태\Dynamic\06_river\koem_climate_to06to15_v5_sig_gy_9points.mat');
clearvars *PAR* biof*
% Kd
biof=0.02486*nanmean(clim1.obm_gy_chl_ft_1(:,3));
biof2=0.02486*nanmean(clim1.obm_gy_chl_ft_2(:,3));
biof3=0.02486*nanmean(clim1.obm_gy_chl_ft_3(:,3));
for k = 1:20
if k==1
BPAR(k) = (1-exp(-(kd1+biof)))./(kd1+biof)
else
BPAR(k) = BPAR(k-1) * (1-exp(-(kd1+biof)))./(kd1+biof)
end
end

for k = 1:20
if k==1
BPAR2(k) = (1-exp(-(kd2+biof2)))/(kd2+biof2)
else
BPAR2(k) = BPAR2(k-1) * (1-exp(-(kd2+biof2)))/(kd2+biof2)
end
end

for k = 1:20
if k==1
BPAR3(k) = (1-exp(-(kd3+biof3)))/(kd3+biof3)
else
BPAR3(k) = BPAR3(k-1) * (1-exp(-(kd3+biof3)))/(kd3+biof3)
end
end

% Kdc
for k = 1:20
if k==1
BPARC(k) = (1-exp(-(kdc1+biof)))/(kdc1+biof)
else
BPARC(k) = BPARC(k-1) * (1-exp(-(kdc1+biof)))/(kdc1+biof)
end
end

for k = 1:20
if k==1
BPARC2(k) = (1-exp(-(kdc2+biof2)))/(kdc2+biof2)
else
BPARC2(k) = BPARC2(k-1) * (1-exp(-(kdc2+biof2)))/(kdc2+biof)
end
end

for k = 1:20
if k==1
BPARC3(k) = (1-exp(-(kdc3+biof3)))/(kdc3+biof3)
else
BPARC3(k) = BPARC3(k-1) * (1-exp(-(kdc3+biof3)))/(kdc3+biof3)
end
end

% Kd
for k = 1:20
if k==1
BPARR(k) = (1-exp(-(0.04+biof)))/(0.04+biof)
else
BPARR(k) = BPARR(k-1) * (1-exp(-(0.04+biof)))/(0.04+biof)
end
end

for k = 1:20
if k==1
BPARR2(k) = (1-exp(-(0.04+biof2)))/(0.04+biof2)
else
BPARR2(k) = BPARR2(k-1) * (1-exp(-(0.04+biof2)))/(0.04+biof2)
end
end

% kd * 2
for k = 1:20
if k==1
BPART(k) = (1-exp(-(kdt2+biof)))/(kdt2+biof)
else
BPART(k) = BPART(k-1) * (1-exp(-(kdt2+biof)))/(kdt2+biof)
end
end

for k = 1:20
if k==1
BPART2(k) = (1-exp(-(kdt2_2+biof2)))/(kdt2_2+biof2)
else
BPART2(k) = BPART2(k-1) * (1-exp(-(kdt2_2+biof2)))/(kdt2_2+biof)
end
end

for k = 1:20
if k==1
BPART3(k) = (1-exp(-(kdt2_3+biof3)))/(kdt2_3+biof3)
else
BPART3(k) = BPART3(k-1) * (1-exp(-(kdt2_3+biof3)))/(kdt2_3+biof3)
end
end

% Kd
for k = 1:20
if k==1
BPARR(k) = (1-exp(-(0.04+biof)))/(0.04+biof)
else
BPARR(k) = BPARR(k-1) * (1-exp(-(0.04+biof)))/(0.04+biof)
end
end

for k = 1:20
if k==1
BPARR2(k) = (1-exp(-(0.04+biof2)))/(0.04+biof2)
else
BPARR2(k) = BPARR2(k-1) * (1-exp(-(0.04+biof2)))/(0.04+biof2)
end
end

for k = 1:20
if k==1
BPARR3(k) = (1-exp(-(0.04+biof3)))/(0.04+biof3)
else
BPARR3(k) = BPARR3(k-1) * (1-exp(-(0.04+biof3)))/(0.04+biof3)
end
end


figure; hold on;
plot([1 BPARR],'r.-'); plot([1 BPARR2],'g.-'); plot([1 BPARR3],'b.-'); grid on;
plot([1 BPAR],'r');  plot([1 BPAR2],'g'); plot([1 BPAR3],'b'); grid on;
plot([1 BPARC],'r--');  plot([1 BPARC2],'g--'); plot([1 BPARC3],'b--'); grid on;
legend('Kd(1st-fennel)','Kd(2nd-fennel)','Kd(3rd-fennel)','Kd(1st)','Kd(2nd)','Kd(3rd)','KdC(1st)','KdC(2nd)','KdC(3rd)')
xlim([1 21]);
xlabel('num. of vertical layer'); ylabel('Normalizaed irradiance')
xticks(1:21); xticklabels(0:20);
text(5,0.5,[num2str(nanmean(clim1.obm_gy_chl_ft_1(:,3)),'%.2f'),' ug/L'],'color','r')
text(5,0.45,[num2str(nanmean(clim1.obm_gy_chl_ft_2(:,3)),'%.2f'),' ug/L'],'color','g')
text(5,0.4,[num2str(nanmean(clim1.obm_gy_chl_ft_3(:,3)),'%.2f'),' ug/L'],'color','b')


figure; hold on;
plot([1 BPARR],'r.-'); plot([1 BPARR2],'g.-'); plot([1 BPARR3],'b.-'); grid on;
plot([1 BPAR],'r');  plot([1 BPAR2],'g'); plot([1 BPAR3],'b'); grid on;
plot([1 BPART],'r--');  plot([1 BPART2],'g--'); plot([1 BPART3],'b--'); grid on;
legend('Kd(1st-fennel)','Kd(2nd-fennel)','Kd(3rd-fennel)','Kd(1st)','Kd(2nd)','Kd(3rd)',[num2str(kd_times),'*Kd(1st)'],[num2str(kd_times),'*Kd(2nd)'],[num2str(kd_times),'*Kd(3rd)'])
xlim([1 21]);
xlabel('num. of vertical layer'); ylabel('Normalizaed irradiance')
xticks(1:21); xticklabels(0:20);
text(5,0.5,[num2str(nanmean(clim1.obm_gy_chl_ft_1(:,3)),'%.2f'),' ug/L'],'color','r')
text(5,0.45,[num2str(nanmean(clim1.obm_gy_chl_ft_2(:,3)),'%.2f'),' ug/L'],'color','g')
text(5,0.4,[num2str(nanmean(clim1.obm_gy_chl_ft_3(:,3)),'%.2f'),' ug/L'],'color','b')


figure; hold on;
plot([1 BPARR],'r.-'); plot([1 BPARR2],'g.-'); plot([1 BPARR3],'b.-'); grid on;
plot([1 BPAR],'r');  plot([1 BPAR2],'g'); plot([1 BPAR3],'b'); grid on;
% plot([1 BPARC],'r--');  plot([1 BPARC2],'g--'); plot([1 BPARC3],'b--'); grid on;
legend('Kd(1st-fennel)','Kd(2nd-fennel)','Kd(3rd-fennel)','Kd(1st)','Kd(2nd)','Kd(3rd)')
xlim([1 21]);
xlabel('Depth (m)'); ylabel('Normalizaed irradiance')
xticks(1:21); xticklabels(0:20);
text(5,0.5,[num2str(nanmean(clim1.obm_gy_chl_ft_1(:,3)),'%.2f'),' ug/L'],'color','r')
text(5,0.45,[num2str(nanmean(clim1.obm_gy_chl_ft_2(:,3)),'%.2f'),' ug/L'],'color','g')
text(5,0.4,[num2str(nanmean(clim1.obm_gy_chl_ft_3(:,3)),'%.2f'),' ug/L'],'color','b')


% Kd
figure; hold on;
plot(clim1.obm_gy_chl_ft_1(:,3),'r','linew',2); plot(clim1.obm_gy_chl_ft_2(:,3),'g','linew',2); 
plot(clim1.obm_gy_chl_ft_3(:,3),'b','linew',2); grid on;
legend('1997~2006','2007~2015','2016~2019');
ylabel('Chl.a (ug/L)'); xlabel('days'); title('KOEM OBS Climate');
xlim([1 365]);



