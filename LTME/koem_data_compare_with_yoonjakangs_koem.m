close all; clear; clc;

load koem_timeseires_monthly_gy_only_raw.mat %not extract outlier (raw data)

% name_tag{sp_gy}
% ±¤¾çÇ×, ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5, ¿©¼ö2, ¿©¼ö3, ¿©¼ö1
% [res I]=sort([4,3,2,1,5]);
P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;

% DIP
gy_csh=regime_po4(sp_gy(2:6),:)./P_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
gy_csh_b=regime_po4_b(sp_gy(2:6),:)./P_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5

gy_csh_all=regime_po4(sp_gy,:)./P_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
gy_csh_b_all=regime_po4_b(sp_gy,:)./P_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5

% gy_csh=regime_po4(sp_gy(2:6),:)./PO4_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
% gy_csh_b=regime_po4_b(sp_gy(2:6),:)./PO4_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
% 
% gy_csh_all=regime_po4(sp_gy,:)./PO4_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
% gy_csh_b_all=regime_po4_b(sp_gy,:)./PO4_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5

%NO3-N
no3_csh=regime_no3(sp_gy(2:6),:)./N_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
no3_csh_b=regime_no3_b(sp_gy(2:6),:)./N_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5

no3_csh_all=regime_no3(sp_gy,:)./N_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
no3_csh_b_all=regime_no3_b(sp_gy,:)./N_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5

%NH4-N
nh4_csh=regime_nh4(sp_gy(2:6),:)./N_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
nh4_csh_b=regime_nh4_b(sp_gy(2:6),:)./N_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5

nh4_csh_all=regime_nh4(sp_gy,:)./N_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
nh4_csh_b_all=regime_nh4_b(sp_gy,:)./N_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5

% nh4_csh=regime_nh4(sp_gy(2:6),:)./NH4_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
% nh4_csh_b=regime_nh4_b(sp_gy(2:6),:)./NH4_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
% 
% nh4_csh_all=regime_nh4(sp_gy,:)./NH4_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5
% nh4_csh_b_all=regime_nh4_b(sp_gy,:)./NH4_MW;  % ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5

% DIP
for i = 1:size(gy_csh,1)
clearvars data_tmp tmp_x
data_tmp=gy_csh(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
gy_csh_in(i,:)=data_tmp;
end

for i = 1:size(gy_csh_b,1)
clearvars data_tmp tmp_x
data_tmp=gy_csh_b(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
gy_csh_b_in(i,:)=data_tmp;
end

for i = 1:size(gy_csh,1)
clearvars data_tmp tmp_x
data_tmp=gy_csh_all(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
gy_csh_all_in(i,:)=data_tmp;
end

for i = 1:size(gy_csh_b,1)
clearvars data_tmp tmp_x
data_tmp=gy_csh_b_all(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
gy_csh_b_all_in(i,:)=data_tmp;
end

%NH4-N
for i = 1:size(nh4_csh,1)
clearvars data_tmp tmp_x
data_tmp=nh4_csh(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
nh4_csh_in(i,:)=data_tmp;
end

for i = 1:size(nh4_csh_b,1)
clearvars data_tmp tmp_x
data_tmp=nh4_csh_b(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
nh4_csh_b_in(i,:)=data_tmp;
end

for i = 1:size(nh4_csh,1)
clearvars data_tmp tmp_x
data_tmp=nh4_csh_all(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
nh4_csh_all_in(i,:)=data_tmp;
end

for i = 1:size(nh4_csh_b,1)
clearvars data_tmp tmp_x
data_tmp=nh4_csh_b_all(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
nh4_csh_b_all_in(i,:)=data_tmp;
end

%NO3-N
for i = 1:size(no3_csh,1)
clearvars data_tmp tmp_x
data_tmp=no3_csh(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
no3_csh_in(i,:)=data_tmp;
end

for i = 1:size(no3_csh_b,1)
clearvars data_tmp tmp_x
data_tmp=no3_csh_b(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
no3_csh_b_in(i,:)=data_tmp;
end

for i = 1:size(no3_csh,1)
clearvars data_tmp tmp_x
data_tmp=no3_csh_all(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
no3_csh_all_in(i,:)=data_tmp;
end

for i = 1:size(no3_csh_b,1)
clearvars data_tmp tmp_x
data_tmp=no3_csh_b_all(i,:);
tmp_x = 1:length(data_tmp);
data_tmp(isnan(data_tmp)) = interp1(tmp_x(~isnan(data_tmp)),data_tmp(~isnan(data_tmp)),tmp_x(isnan(data_tmp))); 
no3_csh_b_all_in(i,:)=data_tmp;
end

[res I]=sort([4,3,2,1,5]);
gy_csh_in_re = gy_csh_in(I,:); %sort to ±¤¾ç1~±¤¾ç5
gy_csh_b_in_re = gy_csh_b_in(I,:); %sort to ±¤¾ç1~±¤¾ç5

no3_csh_in_re = no3_csh_in(I,:); %sort to ±¤¾ç1~±¤¾ç5
no3_csh_b_in_re = no3_csh_b_in(I,:); %sort to ±¤¾ç1~±¤¾ç5

nh4_csh_in_re = nh4_csh_in(I,:); %sort to ±¤¾ç1~±¤¾ç5
nh4_csh_b_in_re = nh4_csh_b_in(I,:); %sort to ±¤¾ç1~±¤¾ç5

figure; 
plot((gy_csh_in_re)');
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');

figure; 
plot((gy_csh_b_in_re)');
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');

% 5-point spatial mean vs. 9-point spatial mean
figure; hold on;
plot(nanmean(gy_csh_in_re,1),'b');
plot(nanmean(gy_csh_all_in,1),'r');
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

% 5-point spatial mean
figure; hold on;
plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(nanmean(gy_csh_in_re,1),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

% yoonja_kang's KOEM data
% 1: year,  2: Month, 3: St., 4: depth,  5: temp, 6: salt, 7: pH, 8: DO, 9:
% COD, 10: NH4, 11: NO2, 12: NO3, 13: NO2NO3, 14: DIN, 15: TN, 16: DIP, 17:
% TP, 18: SIO2, 19: SS, 20: Chla, 21: Sechhi, 22: NP
[raw txt]=xlsread('D:\Àå±â»ýÅÂ\Dynamic\KOEM\koem_yoonja_kang\GB_Environment_Selected_Converted_csh.xlsx','surf');
[rawb txtb]=xlsread('D:\Àå±â»ýÅÂ\Dynamic\KOEM\koem_yoonja_kang\GB_Environment_Selected_Converted_csh.xlsx','bot');

for i = 1:5 %'±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5'
ind(i,:)=find(raw(:,3) == i );
indb(i,:)=find(rawb(:,3) == i );
end

clearvars po4_yjk po4_b_yjk
for i = 1:5
po4_yjk(i,:) = raw(ind(i,:),16);
po4_b_yjk(i,:) = rawb(indb(i,:),16);
no3_yjk(i,:) = raw(ind(i,:),12);
no3_b_yjk(i,:) = rawb(indb(i,:),12);
nh4_yjk(i,:) = raw(ind(i,:),10);
nh4_b_yjk(i,:) = rawb(indb(i,:),10);
end

% ref_date(145:252); %2007~2017 (yoonja kang's time period);
ref_date(146:3:251); % in situ. time for yoonja kang's data


%% compare 5point csh vs. yoonja kang's
figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(nanmean(gy_csh_in_re,1),'color','b','linew',2);
plot(146:3:251,nanmean(po4_yjk,1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
ylabel('PO4-P (mmol P/l)')
xlim([0 264]); grid on;
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(nanmean(no3_csh_in_re,1),'color','b','linew',2);
plot(146:3:251,nanmean(no3_yjk,1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
ylabel('NO3-N (mmol N/l)')
xlim([0 264]); grid on;
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(nanmean(nh4_csh_in_re,1),'color','b','linew',2);
plot(146:3:251,nanmean(nh4_yjk,1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
ylabel('NH4-N (mmol N/l)')
xlim([0 264]); grid on;
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);


%% compare 5point csh vs. yoonja kang's (metching it)
figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(nanmean(gy_csh_in_re,1).*(P_MW/PO4_MW),'color','b','linew',2);
plot(146:3:251,nanmean(po4_yjk,1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylabel('PO4-P (mmol P/l)')
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(nanmean(no3_csh_in_re,1).*(N_MW/NO3_MW),'color','b','linew',2);
plot(146:3:251,nanmean(no3_yjk,1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylabel('NO3-N (mmol N/l)')
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(nanmean(nh4_csh_in_re,1),'color','b','linew',2);
plot(nanmean(nh4_csh_in_re,1).*(N_MW/NH4_MW),'color','b','linew',2);
plot(146:3:251,nanmean(nh4_yjk,1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylabel('NH4-N (mmol N/l)')
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);


figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(nanmean(gy_csh_b_in_re,1),'color','b','linew',2);
plot(146:3:251,nanmean(po4_b_yjk,1),'color','r','linew',2);

xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

%% raw data from yoonjakang

yoonja=load('yoonjakangs_koem_data_monthly.mat');


%% compare 5point yoonja kang's converting data vs. yoonja kang's raw(metching it)
figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(yoonja.po4_sur(2:6,2:3:264),1).*(1/PO4_MW),'color','b','linew',2);
plot(146:3:251,nanmean(po4_yjk,1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylabel('PO4-P (mmol P/l)')
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(yoonja.no3_sur(2:6,2:3:264),1).*(1/NO3_MW),'color','b','linew',2);
plot(146:3:251,nanmean(no3_yjk,1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylabel('NO3-N (mmol N/l)')
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(yoonja.nh4_sur(2:6,2:3:264),1).*(1/NH4_MW),'color','b','linew',2);
plot(146:3:251,nanmean(nh4_yjk,1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylabel('NH4-N (mmol N/l)')
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(yoonja.nh4_sur(2:6,2:3:264),1).*(1/NH4_MW),'color','b','linew',2);
plot(146:3:251,nanmean(nh4_yjk,1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylabel('NH4-N (mmol N/l)')
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);

nanmean(yoonja.po4_sur(2:6,146:3:251),1) ./ nanmean(po4_yjk,1)
nanmean(yoonja.no3_sur(2:6,146:3:251),1) ./ nanmean(no3_yjk,1)
nanmean(yoonja.nh4_sur(2:6,146:3:251),1) ./ nanmean(nh4_yjk,1)




cd D:\Àå±â»ýÅÂ\Dynamic\KOEM\koem_yoonja_kang\pic_compare_raw
%% 9point- mean
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_temp(sp_gy(I),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.temp_sur(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
legend('cshwa','Pf.Kang','fontsize',15);
print(fig,strcat('label.png'),'-dpng')

% temp 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_temp(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.temp_sur(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
print(fig,strcat('temp_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_temp_b(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.temp_bot(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title('temp bottom')
xtickangle(45);
print(fig,strcat('temp_b.png'),'-dpng')

% salt 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_salt(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.salt_sur(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([20 35])
ylabel('salt (psu)')
title('salt surface')
xtickangle(45);
print(fig,strcat('salt_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_salt_b(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.salt_bot(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([20 35])
ylabel('salt (psu)')
title('salt bottom')
xtickangle(45);
print(fig,strcat('salt_b.png'),'-dpng')

% chl 
fig=figure; hold on;
plot(nanmean(regime_chl(sp_gy,:),1),'.','color','r','linew',2);
plot(nanmean(yoonja.chl_sur(:,:),1),'.','color','b','linew',2);
plot(2:3:264,nanmean(regime_chl(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.chl_sur(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 22]);
ylabel('chl (ug/L)')
title('chl surface')
xtickangle(45);
print(fig,strcat('chl_s.png'),'-dpng')

fig=figure; hold on;
plot(nanmean(regime_chl_b(sp_gy,:),1),'.','color','r','linew',2);
plot(nanmean(yoonja.chl_bot(:,:),1),'.','color','b','linew',2);
plot(2:3:264,nanmean(regime_chl_b(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.chl_bot(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 22]);
ylabel('chl (ug/L)')
title('chl bottom')
xtickangle(45);
print(fig,strcat('chl_b.png'),'-dpng')

% nh4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_nh4(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.nh4_sur(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 120])
ylabel('NH4-N (ug/L)')
title('NH4-N surface')
xtickangle(45);
print(fig,strcat('nh4_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_nh4_b(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.nh4_bot(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 120])
ylabel('NH4-N (ug/L)')
title('NH4-N bottom')
xtickangle(45);
print(fig,strcat('nh4_b.png'),'-dpng')

% no3 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_no3(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.no3_sur(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 450])
ylabel('NO3-N (ug/L)')
title('NO3-N surface')
xtickangle(45);
print(fig,strcat('no3_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_no3_b(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.no3_bot(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 450])
ylabel('NO3-N (ug/L)')
title('NO3-N bottom')
xtickangle(45);
print(fig,strcat('no3_b.png'),'-dpng')

% DO
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_do(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.do_sur(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([2 14])
ylabel('DO (mg/L)')
title('DO surface')
xtickangle(45);
print(fig,strcat('do_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_do_b(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.do_bot(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([2 14])
ylabel('DO (mg/L)')
title('DO bottom')
xtickangle(45);
print(fig,strcat('do_b.png'),'-dpng')

% po4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_po4(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.po4_sur(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP surface')
xtickangle(45);
print(fig,strcat('DIP_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_po4_b(sp_gy,2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.po4_bot(:,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP bottom')
xtickangle(45);
print(fig,strcat('DIP_b.png'),'-dpng')

% DIN (ug/L)
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(yoonja.DIN_sur(:,2:3:264),1),'color','b','linew',2);
plot(2:3:264,nanmean(yoonja.DIN_bot(:,2:3:264),1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
% ylim([0 120])
ylabel('DIN (ug/L)')
title('DIN surface & bottom')
legend('surface','bottom')
xtickangle(45);
print(fig,strcat('DIN_raw.png'),'-dpng')

% DIN (umol)
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(yoonja.DIN_sur(:,2:3:264),1) ./14,'color','b','linew',2);
plot(2:3:264,nanmean(yoonja.DIN_bot(:,2:3:264),1) ./14,'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
% ylim([0 120])
ylabel('DIN (umol/L)')
title('DIN surface & bottom')
legend('surface','bottom')
xtickangle(45);
print(fig,strcat('DIN_raw_mol.png'),'-dpng')






%% raw data each st. confirm
% ±¤¾çÇ×, ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5, ¿©¼ö2, ¿©¼ö3, ¿©¼ö1
cd D:\Àå±â»ýÅÂ\Dynamic\KOEM\koem_yoonja_kang\pic_compare_raw\each_st
[res I]=sort([1,5,4,3,2,6,8,9,7]);
% sorting confirm
nn=['±¤¾çÇ×'; '±¤¾ç4'; '±¤¾ç3'; '±¤¾ç2'; '±¤¾ç1'; '±¤¾ç5'; '¿©¼ö2'; '¿©¼ö3'; '¿©¼ö1'];
nn(I,:)

fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_temp(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.temp_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
legend('cshwa','Pf.Kang','fontsize',15);
print(fig,strcat('label.png'),'-dpng')

for j = 1:1 %st.
% temp 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_temp(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.temp_sur(j,2:3:264),'.','color','b','linew',2);
plot(2:3:264,regime_temp(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.temp_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title(['st-',num2str(j),'temp surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_temp_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_temp_b(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.temp_bot(j,2:3:264),'.','color','b','linew',2);
plot(2:3:264,regime_temp_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.temp_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title(['st-',num2str(j),'temp bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_temp_b.png']),'-dpng')

% salt 
fig=figure; hold on;
plot(2:3:264,regime_salt(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.salt_sur(j,2:3:264),'.','color','b','linew',2);
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_salt(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.salt_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([5 35])
ylabel('salt (psu)')
title(['st-',num2str(j),'salt surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_salt_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_salt_b(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.salt_bot(j,2:3:264),'.','color','b','linew',2);
plot(2:3:264,regime_salt_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.salt_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([5 35])
ylabel('salt (psu)')
title(['st-',num2str(j),'salt bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_salt_b.png']),'-dpng')

% chl 
fig=figure; hold on;
plot(regime_chl(sp_gy(I(j)),:),'.','color','r','linew',2);
plot(yoonja.chl_sur(j,:),'.','color','b','linew',2);
plot(2:3:264,regime_chl(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.chl_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 38]);
ylabel('chl (ug/L)')
title(['st-',num2str(j),'chl surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_chl_s.png']),'-dpng')

fig=figure; hold on;
plot(regime_chl_b(sp_gy(I(j)),:),'.','color','r','linew',2);
plot(yoonja.chl_bot(j,:),'.','color','b','linew',2);
plot(2:3:264,regime_chl_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.chl_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 38]);
ylabel('chl (ug/L)')
title(['st-',num2str(j),'chl bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_chl_b.png']),'-dpng')

% nh4
fig=figure; hold on;
plot(2:3:264,regime_nh4(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.nh4_sur(j,2:3:264),'.','color','b','linew',2);
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_nh4(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.nh4_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 240])
ylabel('NH4-N (ug/L)')
title(['st-',num2str(j),'NH4-N surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_nh4_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_nh4_b(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.nh4_bot(j,2:3:264),'.','color','b','linew',2);
plot(2:3:264,regime_nh4_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.nh4_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 240])
ylabel('NH4-N (ug/L)')
title(['st-',num2str(j),'NH4-N bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_nh4_b.png']),'-dpng')

% no3 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_no3(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.no3_sur(j,2:3:264),'.','color','b','linew',2);
plot(2:3:264,regime_no3(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.no3_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 450])
ylabel('NO3-N (ug/L)')
title(['st-',num2str(j),'NO3-N surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_no3_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_no3_b(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.no3_bot(j,2:3:264),'.','color','b','linew',2);
plot(2:3:264,regime_no3_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.no3_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 450])
ylabel('NO3-N (ug/L)')
title(['st-',num2str(j),'NO3-N bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_no3_b.png']),'-dpng')

% DO
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_do(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.do_sur(j,2:3:264),'.','color','b','linew',2);
plot(2:3:264,regime_do(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.do_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 20])
ylabel('DO (mg/L)')
title(['st-',num2str(j),'DO surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_do_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_do_b(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.do_bot(j,2:3:264),'.','color','b','linew',2);
plot(2:3:264,regime_do_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.do_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 20])
ylabel('DO (mg/L)')
title(['st-',num2str(j),'DO bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_do_b.png']),'-dpng')

% po4
fig=figure; hold on;
plot(2:3:264,regime_po4(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.po4_sur(j,2:3:264),'.','color','b','linew',2);
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_po4(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 200])
ylabel('DIP (ug/L)')
title(['st-',num2str(j),'DIP surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_DIP_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_po4_b(sp_gy(I(j)),2:3:264),'.','color','r','linew',2);
plot(2:3:264,yoonja.po4_bot(j,2:3:264),'.','color','b','linew',2);
plot(2:3:264,regime_po4_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 200])
ylabel('DIP (ug/L)')
title(['st-',num2str(j),'DIP bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_DIP_b.png']),'-dpng')

close all;
end

for j = 2:9 %st.
% temp 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_temp(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.temp_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title(['st-',num2str(j),'temp surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_temp_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_temp_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.temp_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title(['st-',num2str(j),'temp bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_temp_b.png']),'-dpng')

% salt 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_salt(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.salt_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([5 35])
ylabel('salt (psu)')
title(['st-',num2str(j),'salt surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_salt_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_salt_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.salt_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([5 35])
ylabel('salt (psu)')
title(['st-',num2str(j),'salt bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_salt_b.png']),'-dpng')

% chl 
fig=figure; hold on;
plot(regime_chl(sp_gy(I(j)),:),'.','color','r','linew',2);
plot(yoonja.chl_sur(j,:),'.','color','b','linew',2);
plot(2:3:264,regime_chl(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.chl_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 38]);
ylabel('chl (ug/L)')
title(['st-',num2str(j),'chl surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_chl_s.png']),'-dpng')

fig=figure; hold on;
plot(regime_chl_b(sp_gy(I(j)),:),'.','color','r','linew',2);
plot(yoonja.chl_bot(j,:),'.','color','b','linew',2);
plot(2:3:264,regime_chl_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.chl_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 38]);
ylabel('chl (ug/L)')
title(['st-',num2str(j),'chl bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_chl_b.png']),'-dpng')

% nh4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_nh4(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.nh4_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 240])
ylabel('NH4-N (ug/L)')
title(['st-',num2str(j),'NH4-N surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_nh4_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_nh4_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.nh4_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 240])
ylabel('NH4-N (ug/L)')
title(['st-',num2str(j),'NH4-N bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_nh4_b.png']),'-dpng')

% no3 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_no3(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.no3_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 450])
ylabel('NO3-N (ug/L)')
title(['st-',num2str(j),'NO3-N surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_no3_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_no3_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.no3_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 450])
ylabel('NO3-N (ug/L)')
title(['st-',num2str(j),'NO3-N bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_no3_b.png']),'-dpng')

% DO
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_do(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.do_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 20])
ylabel('DO (mg/L)')
title(['st-',num2str(j),'DO surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_do_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_do_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.do_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 20])
ylabel('DO (mg/L)')
title(['st-',num2str(j),'DO bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_do_b.png']),'-dpng')

% po4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_po4(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 200])
ylabel('DIP (ug/L)')
title(['st-',num2str(j),'DIP surface'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_DIP_s.png']),'-dpng')

fig=figure; hold on;
plot(2:3:264,regime_po4_b(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_bot(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 200])
ylabel('DIP (ug/L)')
title(['st-',num2str(j),'DIP bottom'])
xtickangle(45);
print(fig,strcat(['st_',num2str(j),'_DIP_b.png']),'-dpng')

close all;
end

%% 5point-mean
cd D:\Àå±â»ýÅÂ\Dynamic\KOEM\koem_yoonja_kang\pic_compare_raw\5-point-mean
% temp 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_temp(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.temp_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 32])
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
print(fig,strcat('temp_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_temp_b(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.temp_bot(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 32])
ylabel('temp (^oC)')
title('temp bottom')
xtickangle(45);
print(fig,strcat('temp_b.png'),'-dpng')

% salt 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_salt(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.salt_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([16 35])
ylabel('salt (psu)')
title('salt surface')
xtickangle(45);
print(fig,strcat('salt_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_salt_b(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.salt_bot(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([16 35])
ylabel('salt (psu)')
title('salt bottom')
xtickangle(45);
print(fig,strcat('salt_b.png'),'-dpng')

% chl 
fig=figure; hold on;
plot(nanmean(regime_chl(sp_gy(2:6),:),1),'.','color','r','linew',2);
plot(nanmean(yoonja.chl_sur(2:6,:),1),'.','color','b','linew',2);
plot(2:3:264,nanmean(regime_chl(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 22]);
ylabel('chl (ug/L)')
title('chl surface')
xtickangle(45);
print(fig,strcat('chl_s.png'),'-dpng')

fig=figure; hold on;
plot(nanmean(regime_chl_b(sp_gy(2:6),:),1),'.','color','r','linew',2);
plot(nanmean(yoonja.chl_bot(2:6,:),1),'.','color','b','linew',2);
plot(2:3:264,nanmean(regime_chl_b(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.chl_bot(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 22]);
ylabel('chl (ug/L)')
title('chl bottom')
xtickangle(45);
print(fig,strcat('chl_b.png'),'-dpng')

% nh4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_nh4(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.nh4_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 120])
ylabel('NH4-N (ug/L)')
title('NH4-N surface')
xtickangle(45);
print(fig,strcat('nh4_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_nh4_b(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.nh4_bot(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 120])
ylabel('NH4-N (ug/L)')
title('NH4-N bottom')
xtickangle(45);
print(fig,strcat('nh4_b.png'),'-dpng')

% no3 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_no3(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.no3_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 475])
ylabel('NO3-N (ug/L)')
title('NO3-N surface')
xtickangle(45);
print(fig,strcat('no3_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_no3_b(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.no3_bot(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 475])
ylabel('NO3-N (ug/L)')
title('NO3-N bottom')
xtickangle(45);
print(fig,strcat('no3_b.png'),'-dpng')

% DO
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_do(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.do_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([2 15])
ylabel('DO (mg/L)')
title('DO surface')
xtickangle(45);
print(fig,strcat('do_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_do_b(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.do_bot(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([2 15])
ylabel('DO (mg/L)')
title('DO bottom')
xtickangle(45);
print(fig,strcat('do_b.png'),'-dpng')

% po4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,nanmean(regime_po4(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.po4_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP surface')
xtickangle(45);
print(fig,strcat('DIP_s.png'),'-dpng')

fig=figure; hold on;
plot(2:3:264,nanmean(regime_po4_b(sp_gy(2:6),2:3:264),1),'color','r','linew',2);
plot(2:3:264,nanmean(yoonja.po4_bot(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP bottom')
xtickangle(45);
print(fig,strcat('DIP_b.png'),'-dpng')

%% st. compare 
cd D:\Àå±â»ýÅÂ\Dynamic\KOEM\koem_yoonja_kang\pic_compare_raw\st_compare
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_temp(sp_gy(I),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.temp_sur(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
% legend('GY1','Pf.Kang','fontsize',15);
print(fig,strcat('label.png'),'-dpng')

% temp 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_temp(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.temp_sur(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
print(fig,strcat('temp_s.png'),'-dpng')

fig=figure; hold on;
% plot(2:3:264,regime_temp_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.temp_bot(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title('temp bottom')
xtickangle(45);
print(fig,strcat('temp_b.png'),'-dpng')

% salt 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_salt(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.salt_sur(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([20 35])
ylabel('salt (psu)')
title('salt surface')
xtickangle(45);
print(fig,strcat('salt_s.png'),'-dpng')

fig=figure; hold on;
% plot(2:3:264,regime_salt_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.salt_bot(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([20 35])
ylabel('salt (psu)')
title('salt bottom')
xtickangle(45);
print(fig,strcat('salt_b.png'),'-dpng')

% chl 
fig=figure; hold on;
% plot(regime_chl(sp_gy,:),'.','color','r','linew',2);
plot(yoonja.chl_sur(:,:)','.','linew',2);
% plot(2:3:264,regime_chl(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.chl_sur(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 22]);
ylabel('chl (ug/L)')
title('chl surface')
xtickangle(45);
print(fig,strcat('chl_s.png'),'-dpng')

fig=figure; hold on;
% plot(regime_chl_b(sp_gy,:),'.','color','r','linew',2);
plot(yoonja.chl_bot(:,:)','.','linew',2);
% plot(2:3:264,regime_chl_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.chl_bot(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 22]);
ylabel('chl (ug/L)')
title('chl bottom')
xtickangle(45);
print(fig,strcat('chl_b.png'),'-dpng')

% nh4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_nh4(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.nh4_sur(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 120])
ylabel('NH4-N (ug/L)')
title('NH4-N surface')
xtickangle(45);
print(fig,strcat('nh4_s.png'),'-dpng')

fig=figure; hold on;
% plot(2:3:264,regime_nh4_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.nh4_bot(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 120])
ylabel('NH4-N (ug/L)')
title('NH4-N bottom')
xtickangle(45);
print(fig,strcat('nh4_b.png'),'-dpng')

% no3 
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_no3(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.no3_sur(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 450])
ylabel('NO3-N (ug/L)')
title('NO3-N surface')
xtickangle(45);
print(fig,strcat('no3_s.png'),'-dpng')

fig=figure; hold on;
% plot(2:3:264,regime_no3_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.no3_bot(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 450])
ylabel('NO3-N (ug/L)')
title('NO3-N bottom')
xtickangle(45);
print(fig,strcat('no3_b.png'),'-dpng')

% DO
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_do(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.do_sur(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([2 14])
ylabel('DO (mg/L)')
title('DO surface')
xtickangle(45);
print(fig,strcat('do_s.png'),'-dpng')

fig=figure; hold on;
% plot(2:3:264,regime_do_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.do_bot(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([2 14])
ylabel('DO (mg/L)')
title('DO bottom')
xtickangle(45);
print(fig,strcat('do_b.png'),'-dpng')

% po4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_po4(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_sur(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP surface')
xtickangle(45);
print(fig,strcat('DIP_s.png'),'-dpng')

fig=figure; hold on;
% plot(2:3:264,regime_po4_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_bot(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP bottom')
xtickangle(45);
print(fig,strcat('DIP_b.png'),'-dpng')


% po4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_po4(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_sur(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 72]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP surface')
xtickangle(45);
print(fig,strcat('DIP_72_s.png'),'-dpng')

fig=figure; hold on;
% plot(2:3:264,regime_po4_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_bot(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 72]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP bottom')
xtickangle(45);
print(fig,strcat('DIP_72_b.png'),'-dpng')


% po4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_po4(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_sur(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 72]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP surface')
xtickangle(45);
print(fig,strcat('DIP_72_s.png'),'-dpng')

fig=figure; hold on;
% plot(2:3:264,regime_po4_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_bot(:,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 72]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP bottom')
xtickangle(45);
print(fig,strcat('DIP_72_b.png'),'-dpng')

fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_po4(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_sur(:,2:3:264),'color',[.5,.5,.5],'linew',2);
plot(2:3:264,nanmean(yoonja.po4_sur(:,2:3:264),1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 72]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP surface')
xtickangle(45);
print(fig,strcat('DIP_72black_s.png'),'-dpng')

fig=figure; hold on;
% plot(2:3:264,regime_po4_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_bot(:,2:3:264),'color',[.5,.5,.5],'linew',2);
plot(2:3:264,nanmean(yoonja.po4_bot(:,2:3:264),1),'color','r','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 72]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP bottom')
xtickangle(45);
print(fig,strcat('DIP_72black_b.png'),'-dpng')



% po4
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
% plot(2:3:264,regime_po4(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_sur(1,2:3:264),'.','linew',2);
plot(2:3:264,yoonja.po4_sur(2:9,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 72]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP surface')
xtickangle(45);
% print(fig,strcat('DIP_72_s.png'),'-dpng')

fig=figure; hold on;
% plot(2:3:264,regime_po4_b(sp_gy,2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.po4_bot(1,2:3:264),'.','linew',2);
plot(2:3:264,yoonja.po4_bot(2:9,2:3:264),'linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 72]); grid on;
ylim([0 140])
ylabel('DIP (ug/L)')
title('DIP bottom')
xtickangle(45);
% print(fig,strcat('DIP_72_b.png'),'-dpng')


%% KODC 
close all; clear; clc

% 20501 : PO4 is no regime on surf and bot.
kodc205 = load('KODC_data_monthly_v5_20501.mat') %(umol / L)
x = 1:length(kodc205.ref_yymm);
fig=figure; hold on;
plot(x(~isnan(kodc205.po4_sur)),kodc205.po4_sur(~isnan(kodc205.po4_sur)),'linew',1.5);
plot(kodc205.po4_sur,'r.','linew',2);
xticklabels(1987:2:2019);
xticks(85:24:length(kodc205.ref_yymm));
xlim([85 length(kodc205.ref_yymm)]); grid on;
ylim([0 1.3])
ylabel('PO4-P (umol/L)')
title('Á¤¼±ÇØ¾ç°üÃø 205-01 PO4-P surface')
xtickangle(45);
set(gca,'fontweight','bold','fontsize',13);set(gca,'linewidth',1.5); 
print(fig,strcat('po4_kodc_20501_sur.png'),'-dpng')

% 40016 : PO4 is no regime for surf but bot had it (shift 2011-03) 
kodc400 = load('KODC_data_monthly_v5_40016.mat') %(umol / L)
x = 1:length(kodc400.ref_yymm);
fig=figure; hold on;
plot(x(~isnan(kodc400.po4_sur)),kodc400.po4_sur(~isnan(kodc400.po4_sur)),'linew',1.5);
plot(kodc400.po4_sur,'r.','linew',2);
xticklabels(1987:2:2019);
xticks(85:24:length(kodc400.ref_yymm));
xlim([85 length(kodc400.ref_yymm)]); grid on;
ylim([0 1.3])
ylabel('PO4-P (umol/L)')
title('Á¤¼±ÇØ¾ç°üÃø 400-16 PO4-P surface')
xtickangle(45);
set(gca,'fontweight','bold','fontsize',13);set(gca,'linewidth',1.5); 
print(fig,strcat('po4_kodc_40016_sur.png'),'-dpng')

% 20501 : PO4 is no regime on botf and bot.
kodc205 = load('KODC_data_monthly_v5_20501.mat') %(umol / L)
x = 1:length(kodc205.ref_yymm);
fig=figure; hold on;
plot(x(~isnan(kodc205.po4_bot)),kodc205.po4_bot(~isnan(kodc205.po4_bot)),'linew',1.5);
plot(kodc205.po4_bot,'r.','linew',2);
xticklabels(1987:2:2019);
xticks(85:24:length(kodc205.ref_yymm));
xlim([85 length(kodc205.ref_yymm)]); grid on;
ylim([0 1.3])
ylabel('PO4-P (umol/L)')
title('Á¤¼±ÇØ¾ç°üÃø 205-01 PO4-P bottom')
xtickangle(45);
set(gca,'fontweight','bold','fontsize',13);set(gca,'linewidth',1.5); 
print(fig,strcat('po4_kodc_20501_bot.png'),'-dpng')

% 40016 : PO4 is no regime for botf but bot had it (shift 2011-03) 
kodc400 = load('KODC_data_monthly_v5_40016.mat') %(umol / L)
x = 1:length(kodc400.ref_yymm);
fig=figure; hold on;
plot(x(~isnan(kodc400.po4_bot)),kodc400.po4_bot(~isnan(kodc400.po4_bot)),'linew',1.5);
plot(kodc400.po4_bot,'r.','linew',2);
xticklabels(1987:2:2019);
xticks(85:24:length(kodc400.ref_yymm));
xlim([85 length(kodc400.ref_yymm)]); grid on;
ylim([0 1.3])
ylabel('PO4-P (umol/L)')
title('Á¤¼±ÇØ¾ç°üÃø 400-16 PO4-P bottom')
xtickangle(45);
set(gca,'fontweight','bold','fontsize',13);set(gca,'linewidth',1.5); 
print(fig,strcat('po4_kodc_40016_bot.png'),'-dpng')

% no3
% 20501 : no3 is no regime on surf and bot.
kodc205 = load('KODC_data_monthly_v5_20501.mat') %(umol / L)
x = 1:length(kodc205.ref_yymm);
fig=figure; hold on;
plot(x(~isnan(kodc205.no3_sur)),kodc205.no3_sur(~isnan(kodc205.no3_sur)),'linew',1.5);
plot(kodc205.no3_sur,'r.','linew',2);
xticklabels(1987:2:2019);
xticks(85:24:length(kodc205.ref_yymm));
xlim([85 length(kodc205.ref_yymm)]); grid on;
%ylim([0 1.3])
ylabel('NO3-N (umol/L)')
title('Á¤¼±ÇØ¾ç°üÃø 205-01 NO3-N surface')
xtickangle(45);
set(gca,'fontweight','bold','fontsize',13);set(gca,'linewidth',1.5); 
print(fig,strcat('no3_kodc_20501_sur.png'),'-dpng')

% 40016 : no3 is no regime for surf but bot had it (shift 2011-03) 
kodc400 = load('KODC_data_monthly_v5_40016.mat') %(umol / L)
x = 1:length(kodc400.ref_yymm);
fig=figure; hold on;
plot(x(~isnan(kodc400.no3_sur)),kodc400.no3_sur(~isnan(kodc400.no3_sur)),'linew',1.5);
plot(kodc400.no3_sur,'r.','linew',2);
xticklabels(1987:2:2019);
xticks(85:24:length(kodc400.ref_yymm));
xlim([85 length(kodc400.ref_yymm)]); grid on;
%ylim([0 1.3])
ylabel('NO3-N (umol/L)')
title('Á¤¼±ÇØ¾ç°üÃø 400-16 NO3-N surface')
xtickangle(45);
set(gca,'fontweight','bold','fontsize',13);set(gca,'linewidth',1.5); 
print(fig,strcat('no3_kodc_40016_sur.png'),'-dpng')

% 20501 : no3 is no regime on botf and bot.
kodc205 = load('KODC_data_monthly_v5_20501.mat') %(umol / L)
x = 1:length(kodc205.ref_yymm);
fig=figure; hold on;
plot(x(~isnan(kodc205.no3_bot)),kodc205.no3_bot(~isnan(kodc205.no3_bot)),'linew',1.5);
plot(kodc205.no3_bot,'r.','linew',2);
xticklabels(1987:2:2019);
xticks(85:24:length(kodc205.ref_yymm));
xlim([85 length(kodc205.ref_yymm)]); grid on;
ylim([0 35])
ylabel('NO3-N (umol/L)')
title('Á¤¼±ÇØ¾ç°üÃø 205-01 NO3-N bottom')
xtickangle(45);
set(gca,'fontweight','bold','fontsize',13);set(gca,'linewidth',1.5); 
print(fig,strcat('no3_kodc_20501_bot.png'),'-dpng')

% 40016 : no3 is no regime for botf but bot had it (shift 2011-03) 
kodc400 = load('KODC_data_monthly_v5_40016.mat') %(umol / L)
x = 1:length(kodc400.ref_yymm);
fig=figure; hold on;
plot(x(~isnan(kodc400.no3_bot)),kodc400.no3_bot(~isnan(kodc400.no3_bot)),'linew',1.5);
plot(kodc400.no3_bot,'r.','linew',2);
xticklabels(1987:2:2019);
xticks(85:24:length(kodc400.ref_yymm));
xlim([85 length(kodc400.ref_yymm)]); grid on;
ylim([0 35])
ylabel('NO3-N (umol/L)')
title('Á¤¼±ÇØ¾ç°üÃø 400-16 NO3-N bottom')
xtickangle(45);
set(gca,'fontweight','bold','fontsize',13);set(gca,'linewidth',1.5); 
print(fig,strcat('no3_kodc_40016_bot.png'),'-dpng')


