close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\환경과학원
[raw txt]=xlsread('수질측정지점_구례_fix.xls','검색결과','');
% [raw txt]=xlsread('수질측정지점_사천천.xls','검색결과','');

dash_c = '.';
% r_txt_ud = flipud(txt);
r_txt_ud = txt;
r_date_txt=[char(r_txt_ud(2:end,5))];

r_do_txt=[r_txt_ud(2:end,9)];
for i = 1:length(r_do_txt)
    if strcmp(r_do_txt{i,1},'') == 1
       r_do(i) = NaN;
    elseif strcmp(r_do_txt{i,1},'') == 0
       r_do(i) = str2num(char(r_do_txt{i,1}));
    end
end

r_chl_txt=[r_txt_ud(2:end,16)];
for i = 1:length(r_chl_txt)
    if strcmp(r_chl_txt{i,1},'') == 1
       r_chl(i) = NaN;
    elseif strcmp(r_chl_txt{i,1},'') == 0
       r_chl(i) = str2num(char(r_chl_txt{i,1}));
    end
end


r_no3_txt=[r_txt_ud(2:end,20)];
for i = 1:length(r_no3_txt)
    if strcmp(r_no3_txt{i,1},'') == 1
       r_no3(i) = NaN;
    elseif strcmp(r_no3_txt{i,1},'') == 0
       r_no3(i) = str2num(char(r_no3_txt{i,1}));
    end
end

r_nh4_txt=[r_txt_ud(2:end,21)];
for i = 1:length(r_nh4_txt)
    if strcmp(r_nh4_txt{i,1},'') == 1
       r_nh4(i) = NaN;
    elseif strcmp(r_nh4_txt{i,1},'') == 0
       r_nh4(i) = str2num(char(r_nh4_txt{i,1}));
    end
end


for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
    r_raw_txt_c{i,1} = r_date_txt(i,:); % raw
end

%% make 1989~present
k=0
for i = 1989:2019
            k=k+1;
    for j = 1:12
        eom_d_raw(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

%make t-axis for yearly 
temp_indx_pre=[];
for i = 1:size(eom_d_raw,1); temp_indx_pre = [temp_indx_pre eom_d_raw(i,:)]; end;

for i = 1:length(temp_indx_pre)
    if i ==1
         t_raw_indx(i) = temp_indx_pre(i);
    else
        t_raw_indx(i) = sum(temp_indx_pre(1:i));
    end
end

k=0; m=0;
for i = 1:31
    l=0
    for n = 1:12
        m = m+1;
    for j = 1:eom_d_raw(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        yymmdd_txt(k,:)=[num2str(i+1988) '.' num2str(n,'%02d') '.'  num2str(j,'%02d')];
    end
    end
end


% make it to cell-array
for i = 1:length(yymmdd_txt)
    yymmdd_txt_c{i,1} = yymmdd_txt(i,:); % raw
end

% pick matched date from water temp date
for i = 1:length(yymmdd_txt_c)
       indx_raw{i} = find([strcmp(yymmdd_txt_c{i}, r_raw_txt_c)] == 1);
end



%% make 366 mm-dd
for i = 1:12
  eom_d(i) = eomday(1980,i); % 1980 is leap-yr
end

k=0
for i = 1:12
    l=0
    for j = 1:eom_d(i)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mth_d_txt(k,:)=[num2str(i,'%02d') '.'  num2str(l,'%02d')]
    end
end

% make it to cell-array
for i = 1:length(mth_d_txt)
    mth_d_txt_c{i,1} = mth_d_txt(i,:); % delete year
end


%raw
for i = 1:length(yymmdd_txt_c)
    if size(indx_raw{i},1) == 0
        raw_do(i) = NaN;
        raw_chl(i) = NaN; 
        raw_no3(i) = NaN; 
        raw_nh4(i) = NaN; 
        r_raw_date_c{i} = '';
    else
        raw_do(i) = nanmean(r_do(indx_raw{i}));
        raw_chl(i) = nanmean(r_chl(indx_raw{i}));
        raw_no3(i) = nanmean(r_no3(indx_raw{i}));
        raw_nh4(i) = nanmean(r_nh4(indx_raw{i}));
        r_raw_date_c{i} = r_date_txt(indx_raw{i},6:end);  %remove year
    end
end

% make 1989~present
k=0
for i = 1989:2019
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

%make t-axis for yearly 
t_indx_pre=[];
for i = 1:size(eom_d,1); t_indx_pre = [t_indx_pre eom_d(i,:)]; end;

for i = 1:length(t_indx_pre)
    if i ==1
         t_indx(i) = t_indx_pre(i);
    else
        t_indx(i) = sum(t_indx_pre(1:i));
    end
end

%make t-axis for yearly 
t_tick_pre=sum(eom_d,2);
for i = 1:15
    if i ==1
         t_tick(i) = t_tick_pre(i);
    else
        t_tick(i) = sum(t_tick_pre(1:i));
    end
end
t_tick = [1 t_tick];

%% divide the regime (in case, 2004)

%divide the data
% 1989 ~ 2004
regime_do_pre = raw_do(1:t_tick(end));
regime_chl_pre = raw_chl(1:t_tick(end));
regime_no3_pre = raw_no3(1:t_tick(end));
regime_nh4_pre = raw_nh4(1:t_tick(end));
% 2004 ~ 2019
regime_do_af_pre = raw_do(t_tick(end)+1:end);
regime_chl_af_pre = raw_chl(t_tick(end)+1:end);
regime_no3_af_pre = raw_no3(t_tick(end)+1:end);
regime_nh4_af_pre = raw_nh4(t_tick(end)+1:end);

%% extract over 3sig
sig = 3; %% sigma
% ~2004
clearvars idx1
idx1 = find(isnan(regime_do_pre) == 0);
regime_do=regime_do_pre;
regime_do(find(regime_do_pre > mean(regime_do_pre(idx1)) + sig*std(regime_do_pre(idx1))))=NaN;
regime_do(find(regime_do_pre < mean(regime_do_pre(idx1)) - sig*std(regime_do_pre(idx1))))=NaN;
regm_do = nanmean(regime_do)

clearvars idx1
idx1 = find(isnan(regime_chl_pre) == 0);
regime_chl=regime_chl_pre;
regime_chl(find(regime_chl_pre > mean(regime_chl_pre(idx1)) + sig*std(regime_chl_pre(idx1))))=NaN;
regime_chl(find(regime_chl_pre < mean(regime_chl_pre(idx1)) - sig*std(regime_chl_pre(idx1))))=NaN;
regm_chl = nanmean(regime_chl)

clearvars idx1
idx1 = find(isnan(regime_no3_pre) == 0);
regime_no3=regime_no3_pre;
regime_no3(find(regime_no3_pre > mean(regime_no3_pre(idx1)) + sig*std(regime_no3_pre(idx1))))=NaN;
regime_no3(find(regime_no3_pre < mean(regime_no3_pre(idx1)) - sig*std(regime_no3_pre(idx1))))=NaN;
regm_no3 = nanmean(regime_no3)

clearvars idx1
idx1 = find(isnan(regime_nh4_pre) == 0);
regime_nh4=regime_nh4_pre;
regime_nh4(find(regime_nh4_pre > mean(regime_nh4_pre(idx1)) + sig*std(regime_nh4_pre(idx1))))=NaN;
regime_nh4(find(regime_nh4_pre < mean(regime_nh4_pre(idx1)) - sig*std(regime_nh4_pre(idx1))))=NaN;
regm_nh4 = nanmean(regime_nh4)

% 2004~
clearvars idx1
idx1 = find(isnan(regime_do_af_pre) == 0);
regime_do_af=regime_do_af_pre;
regime_do_af(find(regime_do_af_pre > mean(regime_do_af_pre(idx1)) + sig*std(regime_do_af_pre(idx1))))=NaN;
regime_do_af(find(regime_do_af_pre < mean(regime_do_af_pre(idx1)) - sig*std(regime_do_af_pre(idx1))))=NaN;
regm_do_af = nanmean(regime_do_af)

clearvars idx1
idx1 = find(isnan(regime_chl_af_pre) == 0);
regime_chl_af=regime_chl_af_pre;
regime_chl_af(find(regime_chl_af_pre > mean(regime_chl_af_pre(idx1)) + sig*std(regime_chl_af_pre(idx1))))=NaN;
regime_chl_af(find(regime_chl_af_pre < mean(regime_chl_af_pre(idx1)) - sig*std(regime_chl_af_pre(idx1))))=NaN;
regm_chl_af = nanmean(regime_chl_af)

clearvars idx1
idx1 = find(isnan(regime_no3_af_pre) == 0);
regime_no3_af=regime_no3_af_pre;
regime_no3_af(find(regime_no3_af_pre > mean(regime_no3_af_pre(idx1)) + sig*std(regime_no3_af_pre(idx1))))=NaN;
regime_no3_af(find(regime_no3_af_pre < mean(regime_no3_af_pre(idx1)) - sig*std(regime_no3_af_pre(idx1))))=NaN;
regm_no3_af = nanmean(regime_no3_af)

clearvars idx1
idx1 = find(isnan(regime_nh4_af_pre) == 0);
regime_nh4_af=regime_nh4_af_pre;
regime_nh4_af(find(regime_nh4_af_pre > mean(regime_nh4_af_pre(idx1)) + sig*std(regime_nh4_af_pre(idx1))))=NaN;
regime_nh4_af(find(regime_nh4_af_pre < mean(regime_nh4_af_pre(idx1)) - sig*std(regime_nh4_af_pre(idx1))))=NaN;
regm_nh4_af = nanmean(regime_nh4_af)

%% extract regime mean
regime_do = regime_do - regm_do;
regime_chl = regime_chl - regm_chl;
regime_no3 = regime_no3 - regm_no3;
regime_nh4 = regime_nh4 - regm_nh4;

regime_do_af = regime_do_af - (regm_do_af - regm_do) - regm_do;
regime_chl_af = regime_chl_af - (regm_chl_af - regm_chl) - regm_chl;
regime_no3_af = regime_no3_af - (regm_no3_af - regm_no3) - regm_no3;
regime_nh4_af = regime_nh4_af - (regm_nh4_af - regm_nh4) - regm_nh4;

%% combine
com_do = [regime_do'; regime_do_af';];
com_chl = [regime_chl'; regime_chl_af';];
com_no3 = [regime_no3'; regime_no3_af';];
com_nh4 = [regime_nh4'; regime_nh4_af';];


%divide the date
for i = 1:length(com_do)
regime_date_c{i} = r_raw_date_c{i}; % 1989 ~ 2004
end

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       regime_indx{i} = find([strcmp(mth_d_txt_c{i}, regime_date_c)] == 1)
end


%confirm how may data on each the days
size_c = [];
for i = 1:366; size_c(i)=length(regime_indx{i}); end
bar(size_c)



% make quasi_climate 1989~2004
for i = 1:length(regime_indx)
    if size(regime_indx{i},1) == 0     
        reg_clim_do(i) = NaN;
        reg_clim_chl(i) = NaN; 
        reg_clim_no3(i) = NaN; 
        reg_clim_nh4(i) = NaN; 
    else
        reg_clim_do(i) = nanmean(com_do(regime_indx{i}));
        reg_clim_chl(i) = nanmean(com_chl(regime_indx{i}));
        reg_clim_no3(i) = nanmean(com_no3(regime_indx{i}));
        reg_clim_nh4(i) = nanmean(com_nh4(regime_indx{i}));
    end
end



%DO
clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_do = find(isnan(reg_clim_do)==0);
pf_w_do = polyfit(xp_w_do,reg_clim_do(xp_w_do),3);
yp_w_do_04 = polyval(pf_w_do,xp);
yp_w_do_04 = yp_w_do_04 + regm_do;

%CHL
clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_chl = find(isnan(reg_clim_chl)==0);
pf_w_chl = polyfit(xp_w_chl,reg_clim_chl(xp_w_chl),3);
yp_w_chl_04 = polyval(pf_w_chl,xp);
yp_w_chl_04 = yp_w_chl_04 + regm_chl;

%NO3
clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_no3 = find(isnan(reg_clim_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_clim_no3(xp_w_no3),3);
yp_w_no3_04 = polyval(pf_w_no3,xp);
yp_w_no3_04 = yp_w_no3_04 + regm_no3;

%NH4
clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_nh4 = find(isnan(reg_clim_nh4)==0);
pf_w_nh4 = polyfit(xp_w_nh4,reg_clim_nh4(xp_w_nh4),3);
yp_w_nh4_04 = polyval(pf_w_nh4,xp);
yp_w_nh4_04 = yp_w_nh4_04 + regm_nh4;

return


%% DO
figure;
scatter(1:366,reg_clim_do + regm_do);
hold on
plot(1:366, yp_w_do_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([2 18]) 
xlim([1 366])

%% chl
figure;
scatter(1:366,reg_clim_chl + regm_chl);
hold on
plot(1:366, yp_w_chl_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title('sumjin(songjung)-polynomial chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 40])
xlim([1 366])


%% no3
figure;
scatter(1:366,reg_clim_no3 + regm_no3);
hold on
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 3])
xlim([1 366])

%% nh4
figure;
scatter(1:366,reg_clim_nh4 + regm_nh4);
hold on
plot(1:366, yp_w_nh4_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 1.2])
xlim([1 366])

%% ALL
figure;
hold on;
plot(1:366, yp_w_do_04,'--','color','b','linew',2)
plot(1:366, yp_w_chl_04./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
plot(1:366, yp_w_nh4_04,'--','color','c','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Concentration (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial.','fontsize',13)
grid on;
legend('DO','Chla','no3','nh4');
set(gca,'fontsize',13)
xlim([1 366])

save('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat'); 

%% vs.

yp_adv_do = yp_w_do_04;
yp_adv_chl = yp_w_chl_04;
yp_adv_no3 = yp_w_no3_04;
yp_adv_nh4 = yp_w_nh4_04;

clearvars yp_w_*_04 yp_w_*_af 
load('sumjin(songjung)_polynomial_climate_to2004_3sig.mat','yp_w_*'); % ~2004 : 1st regime, 2004~2018 : 2nd regime
clearvars yp_w_do yp_w_chl yp_w_no3 yp_w_nh4
load('sumjin(songjung)_polynomial_climate_to2004_v2_3sig.mat','yp_w_*'); % ~2018 : whole regime


% DO
figure;
plot(1:366, yp_adv_do,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_do_04,'--','color','k','linew',2)
plot(1:366, yp_w_do_af,'--','color','r','linew',2)
plot(1:366, yp_w_do,'--','color','g','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([6 17])
xlim([1 366])
legend('whole(adv)','1st','2nd','whole');

% CHL
figure;
plot(1:366, yp_adv_chl,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_chl_04,'--','color','k','linew',2)
plot(1:366, yp_w_chl_af,'--','color','r','linew',2)
plot(1:366, yp_w_chl,'--','color','g','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('chl (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 20])
xlim([1 366])
legend('whole(adv)','1st','2nd','whole');


% no3
figure;
plot(1:366, yp_adv_no3,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_no3_04,'--','color','k','linew',2)
plot(1:366, yp_w_no3_af,'--','color','r','linew',2)
plot(1:366, yp_w_no3,'--','color','g','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('no3 (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0.4 2.2])
xlim([1 366])
legend('whole(adv)','1st','2nd','whole');


% nh4
figure;
plot(1:366, yp_adv_nh4,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_nh4_04,'--','color','k','linew',2)
plot(1:366, yp_w_nh4_af,'--','color','r','linew',2)
plot(1:366, yp_w_nh4,'--','color','g','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('nh4 (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 0.5])
xlim([1 366])
legend('whole(adv)','1st','2nd','whole');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NAMGANG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\환경과학원
% [raw txt]=xlsread('수질측정지점_구례_fix.xls','검색결과','');
[raw txt]=xlsread('수질측정지점_남강댐1.xls','검색결과','');

dash_c = '.';
% r_txt_ud = flipud(txt);
r_txt_ud = txt;
r_date_txt=[char(r_txt_ud(2:end,5))];

r_do_txt=[r_txt_ud(2:end,9)];
for i = 1:length(r_do_txt)
    if strcmp(r_do_txt{i,1},'') == 1
       r_do(i) = NaN;
    elseif strcmp(r_do_txt{i,1},'') == 0
       r_do(i) = str2num(char(r_do_txt{i,1}));
    end
end

r_chl_txt=[r_txt_ud(2:end,16)];
for i = 1:length(r_chl_txt)
    if strcmp(r_chl_txt{i,1},'') == 1
       r_chl(i) = NaN;
    elseif strcmp(r_chl_txt{i,1},'') == 0
       r_chl(i) = str2num(char(r_chl_txt{i,1}));
    end
end


r_no3_txt=[r_txt_ud(2:end,20)];
for i = 1:length(r_no3_txt)
    if strcmp(r_no3_txt{i,1},'') == 1
       r_no3(i) = NaN;
    elseif strcmp(r_no3_txt{i,1},'') == 0
       r_no3(i) = str2num(char(r_no3_txt{i,1}));
    end
end

r_nh4_txt=[r_txt_ud(2:end,21)];
for i = 1:length(r_nh4_txt)
    if strcmp(r_nh4_txt{i,1},'') == 1
       r_nh4(i) = NaN;
    elseif strcmp(r_nh4_txt{i,1},'') == 0
       r_nh4(i) = str2num(char(r_nh4_txt{i,1}));
    end
end


for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
    r_raw_txt_c{i,1} = r_date_txt(i,:); % raw
end

%% make 1989~present
k=0
for i = 1989:2019
            k=k+1;
    for j = 1:12
        eom_d_raw(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

%make t-axis for yearly 
temp_indx_pre=[];
for i = 1:size(eom_d_raw,1); temp_indx_pre = [temp_indx_pre eom_d_raw(i,:)]; end;

for i = 1:length(temp_indx_pre)
    if i ==1
         t_raw_indx(i) = temp_indx_pre(i);
    else
        t_raw_indx(i) = sum(temp_indx_pre(1:i));
    end
end

k=0; m=0;
for i = 1:31
    l=0
    for n = 1:12
        m = m+1;
    for j = 1:eom_d_raw(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        yymmdd_txt(k,:)=[num2str(i+1988) '.' num2str(n,'%02d') '.'  num2str(j,'%02d')];
    end
    end
end


% make it to cell-array
for i = 1:length(yymmdd_txt)
    yymmdd_txt_c{i,1} = yymmdd_txt(i,:); % raw
end

% pick matched date from water temp date
for i = 1:length(yymmdd_txt_c)
       indx_raw{i} = find([strcmp(yymmdd_txt_c{i}, r_raw_txt_c)] == 1);
end



%% make 366 mm-dd
for i = 1:12
  eom_d(i) = eomday(1980,i); % 1980 is leap-yr
end

k=0
for i = 1:12
    l=0
    for j = 1:eom_d(i)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mth_d_txt(k,:)=[num2str(i,'%02d') '.'  num2str(l,'%02d')]
    end
end

% make it to cell-array
for i = 1:length(mth_d_txt)
    mth_d_txt_c{i,1} = mth_d_txt(i,:); % delete year
end


%raw
for i = 1:length(yymmdd_txt_c)
    if size(indx_raw{i},1) == 0
        raw_do(i) = NaN;
        raw_chl(i) = NaN; 
        raw_no3(i) = NaN; 
        raw_nh4(i) = NaN; 
        r_raw_date_c{i} = '';
    else
        raw_do(i) = nanmean(r_do(indx_raw{i}));
        raw_chl(i) = nanmean(r_chl(indx_raw{i}));
        raw_no3(i) = nanmean(r_no3(indx_raw{i}));
        raw_nh4(i) = nanmean(r_nh4(indx_raw{i}));
        t_indx_temp_raw = indx_raw{i}; %remove same days (because already average it)
        r_raw_date_c{i} = r_date_txt(t_indx_temp_raw(1),6:end);  %remove year
    end
end

% make 1989~present
k=0
for i = 1989:2019
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

%make t-axis for yearly 
t_indx_pre=[];
for i = 1:size(eom_d,1); t_indx_pre = [t_indx_pre eom_d(i,:)]; end;

for i = 1:length(t_indx_pre)
    if i ==1
         t_indx(i) = t_indx_pre(i);
    else
        t_indx(i) = sum(t_indx_pre(1:i));
    end
end

%make t-axis for yearly 
t_tick_pre=sum(eom_d,2);
for i = 1:15
    if i ==1
         t_tick(i) = t_tick_pre(i);
    else
        t_tick(i) = sum(t_tick_pre(1:i));
    end
end
t_tick = [1 t_tick];

%% divide the regime (in case, 2004)

%divide the data
% 1989 ~ 2004
regime_do_pre = raw_do(1:t_tick(end));
regime_chl_pre = raw_chl(1:t_tick(end));
regime_no3_pre = raw_no3(1:t_tick(end));
regime_nh4_pre = raw_nh4(1:t_tick(end));
% 2004 ~ 2019
regime_do_af_pre = raw_do(t_tick(end)+1:end);
regime_chl_af_pre = raw_chl(t_tick(end)+1:end);
regime_no3_af_pre = raw_no3(t_tick(end)+1:end);
regime_nh4_af_pre = raw_nh4(t_tick(end)+1:end);

%% extract over 3sig
sig = 3; %% sigma
% ~2004
clearvars idx1
idx1 = find(isnan(regime_do_pre) == 0);
regime_do=regime_do_pre;
regime_do(find(regime_do_pre > mean(regime_do_pre(idx1)) + sig*std(regime_do_pre(idx1))))=NaN;
regime_do(find(regime_do_pre < mean(regime_do_pre(idx1)) - sig*std(regime_do_pre(idx1))))=NaN;
regm_do = nanmean(regime_do)

clearvars idx1
idx1 = find(isnan(regime_chl_pre) == 0);
regime_chl=regime_chl_pre;
regime_chl(find(regime_chl_pre > mean(regime_chl_pre(idx1)) + sig*std(regime_chl_pre(idx1))))=NaN;
regime_chl(find(regime_chl_pre < mean(regime_chl_pre(idx1)) - sig*std(regime_chl_pre(idx1))))=NaN;
regm_chl = nanmean(regime_chl)

clearvars idx1
idx1 = find(isnan(regime_no3_pre) == 0);
regime_no3=regime_no3_pre;
regime_no3(find(regime_no3_pre > mean(regime_no3_pre(idx1)) + sig*std(regime_no3_pre(idx1))))=NaN;
regime_no3(find(regime_no3_pre < mean(regime_no3_pre(idx1)) - sig*std(regime_no3_pre(idx1))))=NaN;
regm_no3 = nanmean(regime_no3)

clearvars idx1
idx1 = find(isnan(regime_nh4_pre) == 0);
regime_nh4=regime_nh4_pre;
regime_nh4(find(regime_nh4_pre > mean(regime_nh4_pre(idx1)) + sig*std(regime_nh4_pre(idx1))))=NaN;
regime_nh4(find(regime_nh4_pre < mean(regime_nh4_pre(idx1)) - sig*std(regime_nh4_pre(idx1))))=NaN;
regm_nh4 = nanmean(regime_nh4)

% 2004~
clearvars idx1
idx1 = find(isnan(regime_do_af_pre) == 0);
regime_do_af=regime_do_af_pre;
regime_do_af(find(regime_do_af_pre > mean(regime_do_af_pre(idx1)) + sig*std(regime_do_af_pre(idx1))))=NaN;
regime_do_af(find(regime_do_af_pre < mean(regime_do_af_pre(idx1)) - sig*std(regime_do_af_pre(idx1))))=NaN;
regm_do_af = nanmean(regime_do_af)

clearvars idx1
idx1 = find(isnan(regime_chl_af_pre) == 0);
regime_chl_af=regime_chl_af_pre;
regime_chl_af(find(regime_chl_af_pre > mean(regime_chl_af_pre(idx1)) + sig*std(regime_chl_af_pre(idx1))))=NaN;
regime_chl_af(find(regime_chl_af_pre < mean(regime_chl_af_pre(idx1)) - sig*std(regime_chl_af_pre(idx1))))=NaN;
regm_chl_af = nanmean(regime_chl_af)

clearvars idx1
idx1 = find(isnan(regime_no3_af_pre) == 0);
regime_no3_af=regime_no3_af_pre;
regime_no3_af(find(regime_no3_af_pre > mean(regime_no3_af_pre(idx1)) + sig*std(regime_no3_af_pre(idx1))))=NaN;
regime_no3_af(find(regime_no3_af_pre < mean(regime_no3_af_pre(idx1)) - sig*std(regime_no3_af_pre(idx1))))=NaN;
regm_no3_af = nanmean(regime_no3_af)

clearvars idx1
idx1 = find(isnan(regime_nh4_af_pre) == 0);
regime_nh4_af=regime_nh4_af_pre;
regime_nh4_af(find(regime_nh4_af_pre > mean(regime_nh4_af_pre(idx1)) + sig*std(regime_nh4_af_pre(idx1))))=NaN;
regime_nh4_af(find(regime_nh4_af_pre < mean(regime_nh4_af_pre(idx1)) - sig*std(regime_nh4_af_pre(idx1))))=NaN;
regm_nh4_af = nanmean(regime_nh4_af)

%% extract regime mean
regime_do = regime_do - regm_do;
regime_chl = regime_chl - regm_chl;
regime_no3 = regime_no3 - regm_no3;
regime_nh4 = regime_nh4 - regm_nh4;

regime_do_af = regime_do_af - (regm_do_af - regm_do) - regm_do;
regime_chl_af = regime_chl_af - (regm_chl_af - regm_chl) - regm_chl;
regime_no3_af = regime_no3_af - (regm_no3_af - regm_no3) - regm_no3;
regime_nh4_af = regime_nh4_af - (regm_nh4_af - regm_nh4) - regm_nh4;

%% combine
com_do = [regime_do'; regime_do_af';];
com_chl = [regime_chl'; regime_chl_af';];
com_no3 = [regime_no3'; regime_no3_af';];
com_nh4 = [regime_nh4'; regime_nh4_af';];

%divide the date
for i = 1:length(com_do)
regime_date_c{i} = r_raw_date_c{i}; % 1989 ~
end

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       regime_indx{i} = find([strcmp(mth_d_txt_c{i}, regime_date_c)] == 1)
end


%confirm how may data on each the days
size_c = [];
for i = 1:366; size_c(i)=length(regime_indx{i}); end
bar(size_c)



% make quasi_climate 1989~2004
for i = 1:length(regime_indx)
    if size(regime_indx{i},1) == 0     
        reg_clim_do(i) = NaN;
        reg_clim_chl(i) = NaN; 
        reg_clim_no3(i) = NaN; 
        reg_clim_nh4(i) = NaN; 
    else
        reg_clim_do(i) = nanmean(com_do(regime_indx{i}));
        reg_clim_chl(i) = nanmean(com_chl(regime_indx{i}));
        reg_clim_no3(i) = nanmean(com_no3(regime_indx{i}));
        reg_clim_nh4(i) = nanmean(com_nh4(regime_indx{i}));
    end
end


%DO
clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_do = find(isnan(reg_clim_do)==0);
pf_w_do = polyfit(xp_w_do,reg_clim_do(xp_w_do),3);
yp_w_do_04 = polyval(pf_w_do,xp);
yp_w_do_04 = yp_w_do_04 + regm_do;

%CHL
clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_chl = find(isnan(reg_clim_chl)==0);
pf_w_chl = polyfit(xp_w_chl,reg_clim_chl(xp_w_chl),3);
yp_w_chl_04 = polyval(pf_w_chl,xp);
yp_w_chl_04 = yp_w_chl_04 + regm_chl;

%NO3
clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_no3 = find(isnan(reg_clim_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_clim_no3(xp_w_no3),3);
yp_w_no3_04 = polyval(pf_w_no3,xp);
yp_w_no3_04 = yp_w_no3_04 + regm_no3;

%NH4
clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_nh4 = find(isnan(reg_clim_nh4)==0);
pf_w_nh4 = polyfit(xp_w_nh4,reg_clim_nh4(xp_w_nh4),3);
yp_w_nh4_04 = polyval(pf_w_nh4,xp);
yp_w_nh4_04 = yp_w_nh4_04 + regm_nh4;


%% DO
figure;
scatter(1:366,reg_clim_do + regm_do);
hold on
plot(1:366, yp_w_do_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('Namgang-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([2 18]) 
xlim([1 366])

%% chl
figure;
scatter(1:366,reg_clim_chl + regm_chl);
hold on
plot(1:366, yp_w_chl_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title('Namgang-polynomial chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 40])
xlim([1 366])


%% no3
figure;
scatter(1:366,reg_clim_no3 + regm_no3);
hold on
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title('Namgang-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 3])
xlim([1 366])

%% nh4
figure;
scatter(1:366,reg_clim_nh4 + regm_nh4);
hold on
plot(1:366, yp_w_nh4_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title('Namgang-polynomial nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 1.2])
xlim([1 366])

%% ALL
figure;
hold on;
plot(1:366, yp_w_do_04,'--','color','b','linew',2)
plot(1:366, yp_w_chl_04./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
plot(1:366, yp_w_nh4_04,'--','color','c','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Concentration (mg/L)','fontsize',13)
title('Namgang-polynomial.','fontsize',13)
grid on;
legend('DO','Chla','no3','nh4');
set(gca,'fontsize',13)
xlim([1 366])

save('Namgang_polynomial_climate_to2004_advanced(v4)_3sig.mat'); 


%% vs.

yp_adv_do = yp_w_do_04;
yp_adv_chl = yp_w_chl_04;
yp_adv_no3 = yp_w_no3_04;
yp_adv_nh4 = yp_w_nh4_04;

clearvars yp_w_*_04 yp_w_*_af 
load('Namgang_polynomial_climate_to2004_3sig.mat','yp_w_*'); % ~2004 : 1st regime, 2004~2018 : 2nd regime
clearvars yp_w_do yp_w_chl yp_w_no3 yp_w_nh4
load('Namgang_polynomial_climate_to2004_v2_3sig.mat','yp_w_*'); % ~2018 : whole regime


% DO
figure;
plot(1:366, yp_adv_do,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_do_04,'--','color','k','linew',2)
plot(1:366, yp_w_do_af,'--','color','r','linew',2)
plot(1:366, yp_w_do,'--','color','g','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('Namgang-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([6 17])
xlim([1 366])
legend('whole(adv)','1st','2nd','whole');

% CHL
figure;
plot(1:366, yp_adv_chl,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_chl_04,'--','color','k','linew',2)
plot(1:366, yp_w_chl_af,'--','color','r','linew',2)
plot(1:366, yp_w_chl,'--','color','g','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('chl (mg/L)','fontsize',13)
title('Namgang-polynomial chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 20])
xlim([1 366])
legend('whole(adv)','1st','2nd','whole');


% no3
figure;
plot(1:366, yp_adv_no3,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_no3_04,'--','color','k','linew',2)
plot(1:366, yp_w_no3_af,'--','color','r','linew',2)
plot(1:366, yp_w_no3,'--','color','g','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('no3 (mg/L)','fontsize',13)
title('Namgang-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0.4 2.2])
xlim([1 366])
legend('whole(adv)','1st','2nd','whole');


% nh4
figure;
plot(1:366, yp_adv_nh4,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_nh4_04,'--','color','k','linew',2)
plot(1:366, yp_w_nh4_af,'--','color','r','linew',2)
plot(1:366, yp_w_nh4,'--','color','g','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('nh4 (mg/L)','fontsize',13)
title('Namgang-polynomial nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 0.5])
xlim([1 366])
legend('whole(adv)','1st','2nd','whole');