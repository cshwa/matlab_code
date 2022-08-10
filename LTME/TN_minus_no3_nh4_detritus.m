close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\환경과학원
[raw txt]=xlsread('수질측정지점_구례_fix.xls','검색결과','');
% [raw txt]=xlsread('수질측정지점_사천천.xls','검색결과','');

% 11 col :TN, 20 col :NO3-N, 21 col :NH3-N (mg / L)

dash_c = '.';
% r_txt_ud = flipud(txt);
r_txt_ud = txt;
r_date_txt=[char(r_txt_ud(2:end,5))];

r_tn_txt=[r_txt_ud(2:end,11)];
for i = 1:length(r_tn_txt)
    if strcmp(r_tn_txt{i,1},'') == 1
       r_tn(i) = NaN;
    elseif strcmp(r_tn_txt{i,1},'') == 0
       r_tn(i) = str2num(char(r_tn_txt{i,1}));
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
        raw_tn(i) = NaN;
        raw_no3(i) = NaN; 
        raw_nh4(i) = NaN; 
        r_raw_date_c{i} = '';
    else
        raw_tn(i) = nanmean(r_tn(indx_raw{i}));
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
for i = 1:31
    if i ==1
         t_tick(i) = t_tick_pre(i);
    else
        t_tick(i) = sum(t_tick_pre(1:i));
    end
end
t_tick_all = [1 (t_tick+1)];

%make t-axis for yearly (~2003) for seperating regime
clearvars t_tick
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

raw_det = raw_tn - (raw_no3+raw_nh4);

raw_det(find(raw_det < 0)) = NaN;

%divide the data
% 1989 ~ 2004
regime_tn_pre = raw_tn(1:t_tick(end));
regime_no3_pre = raw_no3(1:t_tick(end));
regime_nh4_pre = raw_nh4(1:t_tick(end));
regime_det_pre = raw_det(1:t_tick(end));
% 2004 ~ 2019
regime_tn_af_pre = raw_tn(t_tick(end)+1:end);
regime_no3_af_pre = raw_no3(t_tick(end)+1:end);
regime_nh4_af_pre = raw_nh4(t_tick(end)+1:end);
regime_det_af_pre = raw_det(t_tick(end)+1:end);


%% extract over 3sig
sig = 3; %% sigma
% ~2004
clearvars idx1
idx1 = find(isnan(regime_tn_pre) == 0);
regime_tn=regime_tn_pre;
regime_tn(find(regime_tn_pre > mean(regime_tn_pre(idx1)) + sig*std(regime_tn_pre(idx1))))=NaN;
regime_tn(find(regime_tn_pre < mean(regime_tn_pre(idx1)) - sig*std(regime_tn_pre(idx1))))=NaN;
regm_tn = nanmean(regime_tn)

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

clearvars idx1
idx1 = find(isnan(regime_det_pre) == 0);
regime_det=regime_det_pre;
regime_det(find(regime_det_pre > mean(regime_det_pre(idx1)) + sig*std(regime_det_pre(idx1))))=NaN;
regime_det(find(regime_det_pre < mean(regime_det_pre(idx1)) - sig*std(regime_det_pre(idx1))))=NaN;
regm_det = nanmean(regime_det)


% 2004~
clearvars idx1
idx1 = find(isnan(regime_tn_af_pre) == 0);
regime_tn_af=regime_tn_af_pre;
regime_tn_af(find(regime_tn_af_pre > mean(regime_tn_af_pre(idx1)) + sig*std(regime_tn_af_pre(idx1))))=NaN;
regime_tn_af(find(regime_tn_af_pre < mean(regime_tn_af_pre(idx1)) - sig*std(regime_tn_af_pre(idx1))))=NaN;
regm_tn_af = nanmean(regime_tn_af)

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


clearvars idx1
idx1 = find(isnan(regime_det_af_pre) == 0);
regime_det_af=regime_det_af_pre;
regime_det_af(find(regime_det_af_pre > mean(regime_det_af_pre(idx1)) + sig*std(regime_det_af_pre(idx1))))=NaN;
regime_det_af(find(regime_det_af_pre < mean(regime_det_af_pre(idx1)) - sig*std(regime_det_af_pre(idx1))))=NaN;
regm_det_af = nanmean(regime_det_af)

%% extract regime mean
regime_tn = regime_tn - regm_tn;
regime_no3 = regime_no3 - regm_no3;
regime_nh4 = regime_nh4 - regm_nh4;
regime_det = regime_det - regm_det;

regime_tn_af = regime_tn_af - (regm_tn_af - regm_tn) - regm_tn;
regime_no3_af = regime_no3_af - (regm_no3_af - regm_no3) - regm_no3;
regime_nh4_af = regime_nh4_af - (regm_nh4_af - regm_nh4) - regm_nh4;
regime_det_af = regime_det_af - (regm_det_af - regm_det) - regm_det;

%% combine
com_tn = [regime_tn'; regime_tn_af';];
com_no3 = [regime_no3'; regime_no3_af';];
com_nh4 = [regime_nh4'; regime_nh4_af';];
com_det = [regime_det'; regime_det_af';];

%divide the date
for i = 1:length(com_tn)
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
        reg_clim_tn(i) = NaN;
        reg_clim_no3(i) = NaN; 
        reg_clim_nh4(i) = NaN; 
        reg_clim_det(i) = NaN; 
    else
        reg_clim_tn(i) = nanmean(com_tn(regime_indx{i}));
        reg_clim_no3(i) = nanmean(com_no3(regime_indx{i}));
        reg_clim_nh4(i) = nanmean(com_nh4(regime_indx{i}));
        reg_clim_det(i) = nanmean(com_det(regime_indx{i}));
    end
end



%tn
clearvars xp_w_tn pf_w_tn
xp = 1:366;
xp_w_tn = find(isnan(reg_clim_tn)==0);
pf_w_tn = polyfit(xp_w_tn,reg_clim_tn(xp_w_tn),3);
yp_w_tn_04 = polyval(pf_w_tn,xp);
yp_w_tn_04 = yp_w_tn_04 + regm_tn;

%NO3
clearvars xp_w_tn pf_w_tn
xp = 1:366;
xp_w_no3 = find(isnan(reg_clim_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_clim_no3(xp_w_no3),3);
yp_w_no3_04 = polyval(pf_w_no3,xp);
yp_w_no3_04 = yp_w_no3_04 + regm_no3;

%NH4
clearvars xp_w_tn pf_w_tn
xp = 1:366;
xp_w_nh4 = find(isnan(reg_clim_nh4)==0);
pf_w_nh4 = polyfit(xp_w_nh4,reg_clim_nh4(xp_w_nh4),3);
yp_w_nh4_04 = polyval(pf_w_nh4,xp);
yp_w_nh4_04 = yp_w_nh4_04 + regm_nh4;

%det
clearvars xp_w_tn pf_w_tn
xp = 1:366;
xp_w_det = find(isnan(reg_clim_det)==0);
pf_w_det = polyfit(xp_w_det,reg_clim_det(xp_w_det),3);
yp_w_det_04 = polyval(pf_w_det,xp);
yp_w_det_04 = yp_w_det_04 + regm_det;

return


%% tn
figure;
scatter(1:366,reg_clim_tn + regm_tn);
hold on
plot(1:366, yp_w_tn_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('tn (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial tn.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 4]) 
xlim([1 366])


%% no3
figure;
scatter(1:366,reg_clim_no3 + regm_no3);
hold on
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('NO3 (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial no3.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 3])
xlim([1 366])

%% nh4
figure;
scatter(1:366,reg_clim_nh4 + regm_nh4);
hold on
plot(1:366, yp_w_nh4_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('NH4 (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial nh4.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 1.2])
xlim([1 366])

%% det
figure;
scatter(1:366,reg_clim_det + regm_det);
hold on
plot(1:366, yp_w_det_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('det (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial det.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 1.2])
xlim([1 366])

%% ALL
figure;
hold on;
plot(1:366, yp_w_tn_04,'--','color','b','linew',2)
plot(1:366, yp_w_tn_04 - (yp_w_no3_04 + yp_w_nh4_04),'--','color','k','linew',2)
% plot(1:366, yp_w_det_04,'--','color','g','linew',2)
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
plot(1:366, yp_w_nh4_04,'--','color','c','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('Concentration (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial.','fontsize',16)
grid on;
% legend('tn','diff(det.)','det','no3','nh4');
legend('tn','diff(det.)','no3','nh4');
set(gca,'fontsize',16)
xlim([1 366])



figure;
hold on;
plot(1:366, yp_w_tn_04,'--','color','b','linew',2)
plot(1:366, yp_w_tn_04 - (yp_w_no3_04 + yp_w_nh4_04),'--','color','k','linew',2)
plot(1:366, yp_w_det_04,'--','color','g','linew',2)
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
plot(1:366, yp_w_nh4_04,'--','color','c','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('Concentration (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial.','fontsize',16)
grid on;
legend('tn','diff(det.)','det','no3','nh4');
set(gca,'fontsize',16)
xlim([1 366])

figure;
plot(raw_tn,'.'); hold on
plot(raw_no3+raw_nh4,'r.');
plot(raw_tn-(raw_no3+raw_nh4),'g.');
xlabel('time','fontsize',16)
ylabel('Concentration (mg/L)','fontsize',16)
yline(0)
ylim([-inf 10]) 
grid on;
xlim([1 length(raw_tn)])
title('sumjin(songjung)- TN, NO3+NH4, diff.','fontsize',16)
legend('TN', 'no3+nh4','diff(det.)')
set(gca,'xtick',[t_tick_all(1:3:end)]);
set(gca,'xticklabel',1989:3:2019,'fontsize',10);
set(gca,'fontsize',16)

save('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig_TN.mat'); 

close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\환경과학원
[raw txt]=xlsread('수질측정지점_구례_fix.xls','검색결과','');
% [raw txt]=xlsread('수질측정지점_사천천.xls','검색결과','');

dash_c = '.';
% r_txt_ud = flipud(txt);
r_txt_ud = txt;
r_date_txt=[char(r_txt_ud(2:end,5))];

r_tn_txt=[r_txt_ud(2:end,11)];
for i = 1:length(r_tn_txt)
    if strcmp(r_tn_txt{i,1},'') == 1
       r_tn(i) = NaN;
    elseif strcmp(r_tn_txt{i,1},'') == 0
       r_tn(i) = str2num(char(r_tn_txt{i,1}));
    end
end

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
        raw_tn(i) = NaN;
        raw_do(i) = NaN;
        raw_chl(i) = NaN; 
        raw_no3(i) = NaN; 
        raw_nh4(i) = NaN; 
        r_raw_date_c{i} = '';
    else
        raw_tn(i) = nanmean(r_tn(indx_raw{i}));
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
for i = 1:31
    if i ==1
         t_tick(i) = t_tick_pre(i);
    else
        t_tick(i) = sum(t_tick_pre(1:i));
    end
end
t_tick_all = [1 (t_tick+1)];

%make t-axis for yearly (~2003) for seperating regime
clearvars t_tick
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

raw_det = raw_tn - (raw_no3+raw_nh4);
raw_det(find(raw_det < 0)) = NaN;

%divide the data
% 1989 ~ 2004
regime_tn = raw_tn(1:t_tick(end));
regime_do = raw_do(1:t_tick(end));
regime_chl = raw_chl(1:t_tick(end));
regime_no3 = raw_no3(1:t_tick(end));
regime_nh4 = raw_nh4(1:t_tick(end));
regime_det = raw_det(1:t_tick(end));
% 2004 ~ 2019
regime_tn_af = raw_tn(t_tick(end)+1:end);
regime_do_af = raw_do(t_tick(end)+1:end);
regime_chl_af = raw_chl(t_tick(end)+1:end);
regime_no3_af = raw_no3(t_tick(end)+1:end);
regime_nh4_af = raw_nh4(t_tick(end)+1:end);
regime_det_af = raw_det(t_tick(end)+1:end);

%divide the date
for i = 1:t_tick(end)
regime_date_c{i} = r_raw_date_c{i}; % 1989 ~ 2004
end

j=0;
for i = t_tick(end)+1:t_tick(end)+length(regime_do_af)
    j=j+1;
regime_date_af_c{j} = r_raw_date_c{i}; %  2004 ~ 2019
end

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       regime_indx{i} = find([strcmp(mth_d_txt_c{i}, regime_date_c)] == 1)
end
for i = 1:length(mth_d_txt_c)
       regime_indx_af{i} = find([strcmp(mth_d_txt_c{i}, regime_date_af_c)] == 1)
end

%confirm how may data on each the days
size_c = [];
for i = 1:366; size_c(i)=length(regime_indx{i}); end
bar(size_c)

size_c = [];
for i = 1:366; size_c(i)=length(regime_indx_af{i}); end
bar(size_c)

% % 1989 ~ 2004
% regime_tn = raw_tn(1:t_tick(end));
% regime_do = raw_do(1:t_tick(end));
% regime_chl = raw_chl(1:t_tick(end));
% regime_no3 = raw_no3(1:t_tick(end));
% regime_nh4 = raw_nh4(1:t_tick(end));
% regime_det = raw_det(1:t_tick(end));
% % 2004 ~ 2019
% regime_tn_af = raw_tn(t_tick(end)+1:end);
% regime_do_af = raw_do(t_tick(end)+1:end);
% regime_chl_af = raw_chl(t_tick(end)+1:end);
% regime_no3_af = raw_no3(t_tick(end)+1:end);
% regime_nh4_af = raw_nh4(t_tick(end)+1:end);
% regime_det_af = raw_det(t_tick(end)+1:end);

% make quasi_climate 1989~2004
for i = 1:length(regime_indx)
    if size(regime_indx{i},1) == 0     
        reg_clim_tn(i) = NaN;
        reg_clim_do(i) = NaN;
        reg_clim_chl(i) = NaN; 
        reg_clim_no3(i) = NaN; 
        reg_clim_nh4(i) = NaN; 
        reg_clim_det(i) = NaN;
    else
        reg_clim_tn(i) = nanmean(regime_tn(regime_indx{i}));
        reg_clim_do(i) = nanmean(regime_do(regime_indx{i}));
        reg_clim_chl(i) = nanmean(regime_chl(regime_indx{i}));
        reg_clim_no3(i) = nanmean(regime_no3(regime_indx{i}));
        reg_clim_nh4(i) = nanmean(regime_nh4(regime_indx{i}));
        reg_clim_det(i) = nanmean(regime_det(regime_indx{i}));
    end
end

% make quasi_climate 2004~2019
for i = 1:length(regime_indx_af)
    if size(regime_indx_af{i},1) == 0  
        reg_af_clim_tn(i) = NaN;
        reg_af_clim_do(i) = NaN;
        reg_af_clim_chl(i) = NaN; 
        reg_af_clim_no3(i) = NaN; 
        reg_af_clim_nh4(i) = NaN; 
        reg_af_clim_det(i) = NaN; 
    else
        reg_af_clim_tn(i) = nanmean(regime_tn_af(regime_indx_af{i}));
        reg_af_clim_do(i) = nanmean(regime_do_af(regime_indx_af{i}));
        reg_af_clim_chl(i) = nanmean(regime_chl_af(regime_indx_af{i}));
        reg_af_clim_no3(i) = nanmean(regime_no3_af(regime_indx_af{i}));
        reg_af_clim_nh4(i) = nanmean(regime_nh4_af(regime_indx_af{i}));
        reg_af_clim_det(i) = nanmean(regime_det_af(regime_indx_af{i}));
    end
end

sig = 3; %% sigma
%det
clearvars idx1
idx1 = find(isnan(reg_clim_det) == 0);
reg_clim_det_te=reg_clim_det;
upper_bc_det(sig) = mean(reg_clim_det(idx1)) + sig*std(reg_clim_det(idx1));
lower_bc_det(sig) = mean(reg_clim_det(idx1)) - sig*std(reg_clim_det(idx1));
reg_clim_det_te(find(reg_clim_det > mean(reg_clim_det(idx1)) + sig*std(reg_clim_det(idx1))))=NaN;
reg_clim_det_te(find(reg_clim_det < mean(reg_clim_det(idx1)) - sig*std(reg_clim_det(idx1))))=NaN;
nanmean(reg_clim_det_te(idx1))

clearvars xp_w_det pf_w_det
xp = 1:366;
xp_w_det = find(isnan(reg_clim_det_te)==0);
pf_w_det = polyfit(xp_w_det,reg_clim_det_te(xp_w_det),3);
yp_w_det_04 = polyval(pf_w_det,xp);

clearvars idx1
idx1 = find(isnan(reg_af_clim_det) == 0);
reg_af_clim_det_te=reg_af_clim_det;
upper_bc_det_af(sig) = mean(reg_af_clim_det(idx1)) + sig*std(reg_af_clim_det(idx1));
lower_bc_det_af(sig) = mean(reg_af_clim_det(idx1)) - sig*std(reg_af_clim_det(idx1));
reg_af_clim_det_te(find(reg_af_clim_det > mean(reg_af_clim_det(idx1)) + sig*std(reg_af_clim_det(idx1))))=NaN;
reg_af_clim_det_te(find(reg_af_clim_det < mean(reg_af_clim_det(idx1)) - sig*std(reg_af_clim_det(idx1))))=NaN;
nanmean(reg_af_clim_det_te(idx1))

clearvars xp_w_det pf_w_det
xp = 1:366;
xp_w_det = find(isnan(reg_af_clim_det_te)==0);
pf_w_det = polyfit(xp_w_det,reg_af_clim_det_te(xp_w_det),3);
yp_w_det_af = polyval(pf_w_det,xp);

%tn
clearvars idx1
idx1 = find(isnan(reg_clim_tn) == 0);
reg_clim_tn_te=reg_clim_tn;
upper_bc_tn(sig) = mean(reg_clim_tn(idx1)) + sig*std(reg_clim_tn(idx1));
lower_bc_tn(sig) = mean(reg_clim_tn(idx1)) - sig*std(reg_clim_tn(idx1));
reg_clim_tn_te(find(reg_clim_tn > mean(reg_clim_tn(idx1)) + sig*std(reg_clim_tn(idx1))))=NaN;
reg_clim_tn_te(find(reg_clim_tn < mean(reg_clim_tn(idx1)) - sig*std(reg_clim_tn(idx1))))=NaN;
nanmean(reg_clim_tn_te(idx1))

clearvars xp_w_tn pf_w_tn
xp = 1:366;
xp_w_tn = find(isnan(reg_clim_tn_te)==0);
pf_w_tn = polyfit(xp_w_tn,reg_clim_tn_te(xp_w_tn),3);
yp_w_tn_04 = polyval(pf_w_tn,xp);

clearvars idx1
idx1 = find(isnan(reg_af_clim_tn) == 0);
reg_af_clim_tn_te=reg_af_clim_tn;
upper_bc_tn_af(sig) = mean(reg_af_clim_tn(idx1)) + sig*std(reg_af_clim_tn(idx1));
lower_bc_tn_af(sig) = mean(reg_af_clim_tn(idx1)) - sig*std(reg_af_clim_tn(idx1));
reg_af_clim_tn_te(find(reg_af_clim_tn > mean(reg_af_clim_tn(idx1)) + sig*std(reg_af_clim_tn(idx1))))=NaN;
reg_af_clim_tn_te(find(reg_af_clim_tn < mean(reg_af_clim_tn(idx1)) - sig*std(reg_af_clim_tn(idx1))))=NaN;
nanmean(reg_af_clim_tn_te(idx1))

clearvars xp_w_tn pf_w_tn
xp = 1:366;
xp_w_tn = find(isnan(reg_af_clim_tn_te)==0);
pf_w_tn = polyfit(xp_w_tn,reg_af_clim_tn_te(xp_w_tn),3);
yp_w_tn_af = polyval(pf_w_tn,xp);

%DO
clearvars idx1
idx1 = find(isnan(reg_clim_do) == 0);
reg_clim_do_te=reg_clim_do;
upper_bc_do(sig) = mean(reg_clim_do(idx1)) + sig*std(reg_clim_do(idx1));
lower_bc_do(sig) = mean(reg_clim_do(idx1)) - sig*std(reg_clim_do(idx1));
reg_clim_do_te(find(reg_clim_do > mean(reg_clim_do(idx1)) + sig*std(reg_clim_do(idx1))))=NaN;
reg_clim_do_te(find(reg_clim_do < mean(reg_clim_do(idx1)) - sig*std(reg_clim_do(idx1))))=NaN;
nanmean(reg_clim_do_te(idx1))

clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_do = find(isnan(reg_clim_do_te)==0);
pf_w_do = polyfit(xp_w_do,reg_clim_do_te(xp_w_do),3);
yp_w_do_04 = polyval(pf_w_do,xp);

clearvars idx1
idx1 = find(isnan(reg_af_clim_do) == 0);
reg_af_clim_do_te=reg_af_clim_do;
upper_bc_do_af(sig) = mean(reg_af_clim_do(idx1)) + sig*std(reg_af_clim_do(idx1));
lower_bc_do_af(sig) = mean(reg_af_clim_do(idx1)) - sig*std(reg_af_clim_do(idx1));
reg_af_clim_do_te(find(reg_af_clim_do > mean(reg_af_clim_do(idx1)) + sig*std(reg_af_clim_do(idx1))))=NaN;
reg_af_clim_do_te(find(reg_af_clim_do < mean(reg_af_clim_do(idx1)) - sig*std(reg_af_clim_do(idx1))))=NaN;
nanmean(reg_af_clim_do_te(idx1))

clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_do = find(isnan(reg_af_clim_do_te)==0);
pf_w_do = polyfit(xp_w_do,reg_af_clim_do_te(xp_w_do),3);
yp_w_do_af = polyval(pf_w_do,xp);

%CHL
clearvars idx1
idx1 = find(isnan(reg_clim_chl) == 0);
reg_clim_chl_te=reg_clim_chl;
upper_bc_chl(sig) = mean(reg_clim_chl(idx1)) + sig*std(reg_clim_chl(idx1));
lower_bc_chl(sig) = mean(reg_clim_chl(idx1)) - sig*std(reg_clim_chl(idx1));
reg_clim_chl_te(find(reg_clim_chl > mean(reg_clim_chl(idx1)) + sig*std(reg_clim_chl(idx1))))=NaN;
reg_clim_chl_te(find(reg_clim_chl < mean(reg_clim_chl(idx1)) - sig*std(reg_clim_chl(idx1))))=NaN;
nanmean(reg_clim_chl_te(idx1))

clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_chl = find(isnan(reg_clim_chl_te)==0);
pf_w_chl = polyfit(xp_w_chl,reg_clim_chl_te(xp_w_chl),3);
yp_w_chl_04 = polyval(pf_w_chl,xp);

clearvars idx1
idx1 = find(isnan(reg_af_clim_chl) == 0);
reg_af_clim_chl_te=reg_af_clim_chl;
upper_bc_chl_af(sig) = mean(reg_af_clim_chl(idx1)) + sig*std(reg_af_clim_chl(idx1));
lower_bc_chl_af(sig) = mean(reg_af_clim_chl(idx1)) - sig*std(reg_af_clim_chl(idx1));
reg_af_clim_chl_te(find(reg_af_clim_chl > mean(reg_af_clim_chl(idx1)) + sig*std(reg_af_clim_chl(idx1))))=NaN;
reg_af_clim_chl_te(find(reg_af_clim_chl < mean(reg_af_clim_chl(idx1)) - sig*std(reg_af_clim_chl(idx1))))=NaN;
nanmean(reg_af_clim_chl_te(idx1))

clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_chl = find(isnan(reg_af_clim_chl_te)==0);
pf_w_chl = polyfit(xp_w_chl,reg_af_clim_chl_te(xp_w_chl),3);
yp_w_chl_af = polyval(pf_w_chl,xp);

%NO3
clearvars idx1
idx1 = find(isnan(reg_clim_no3) == 0);
reg_clim_no3_te=reg_clim_no3;
upper_bc_no3(sig) = mean(reg_clim_no3(idx1)) + sig*std(reg_clim_no3(idx1));
lower_bc_no3(sig) = mean(reg_clim_no3(idx1)) - sig*std(reg_clim_no3(idx1));
reg_clim_no3_te(find(reg_clim_no3 > mean(reg_clim_no3(idx1)) + sig*std(reg_clim_no3(idx1))))=NaN;
reg_clim_no3_te(find(reg_clim_no3 < mean(reg_clim_no3(idx1)) - sig*std(reg_clim_no3(idx1))))=NaN;
nanmean(reg_clim_no3_te(idx1))

clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_no3 = find(isnan(reg_clim_no3_te)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_clim_no3_te(xp_w_no3),3);
yp_w_no3_04 = polyval(pf_w_no3,xp);

clearvars idx1
idx1 = find(isnan(reg_af_clim_no3) == 0);
reg_af_clim_no3_te=reg_af_clim_no3;
upper_bc_no3_af(sig) = mean(reg_af_clim_no3(idx1)) + sig*std(reg_af_clim_no3(idx1));
lower_bc_no3_af(sig) = mean(reg_af_clim_no3(idx1)) - sig*std(reg_af_clim_no3(idx1));
reg_af_clim_no3_te(find(reg_af_clim_no3 > mean(reg_af_clim_no3(idx1)) + sig*std(reg_af_clim_no3(idx1))))=NaN;
reg_af_clim_no3_te(find(reg_af_clim_no3 < mean(reg_af_clim_no3(idx1)) - sig*std(reg_af_clim_no3(idx1))))=NaN;
nanmean(reg_af_clim_no3_te(idx1))

clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_no3 = find(isnan(reg_af_clim_no3_te)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_af_clim_no3_te(xp_w_no3),3);
yp_w_no3_af = polyval(pf_w_no3,xp);

%NH4
clearvars idx1
idx1 = find(isnan(reg_clim_nh4) == 0);
reg_clim_nh4_te=reg_clim_nh4;
upper_bc_nh4(sig) = mean(reg_clim_nh4(idx1)) + sig*std(reg_clim_nh4(idx1));
lower_bc_nh4(sig) = mean(reg_clim_nh4(idx1)) - sig*std(reg_clim_nh4(idx1));
reg_clim_nh4_te(find(reg_clim_nh4 > mean(reg_clim_nh4(idx1)) + sig*std(reg_clim_nh4(idx1))))=NaN;
reg_clim_nh4_te(find(reg_clim_nh4 < mean(reg_clim_nh4(idx1)) - sig*std(reg_clim_nh4(idx1))))=NaN;
nanmean(reg_clim_nh4_te(idx1))

clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_nh4 = find(isnan(reg_clim_nh4_te)==0);
pf_w_nh4 = polyfit(xp_w_nh4,reg_clim_nh4_te(xp_w_nh4),3);
yp_w_nh4_04 = polyval(pf_w_nh4,xp);

clearvars idx1
idx1 = find(isnan(reg_af_clim_nh4) == 0);
reg_af_clim_nh4_te=reg_af_clim_nh4;
upper_bc_nh4_af(sig) = mean(reg_af_clim_nh4(idx1)) + sig*std(reg_af_clim_nh4(idx1));
lower_bc_nh4_af(sig) = mean(reg_af_clim_nh4(idx1)) - sig*std(reg_af_clim_nh4(idx1));
reg_af_clim_nh4_te(find(reg_af_clim_nh4 > mean(reg_af_clim_nh4(idx1)) + sig*std(reg_af_clim_nh4(idx1))))=NaN;
reg_af_clim_nh4_te(find(reg_af_clim_nh4 < mean(reg_af_clim_nh4(idx1)) - sig*std(reg_af_clim_nh4(idx1))))=NaN;
nanmean(reg_af_clim_nh4_te(idx1))

clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_nh4 = find(isnan(reg_af_clim_nh4_te)==0);
pf_w_nh4 = polyfit(xp_w_nh4,reg_af_clim_nh4_te(xp_w_nh4),3);
yp_w_nh4_af = polyval(pf_w_nh4,xp);

return

figure;
hold on;
plot(1:366, yp_w_tn_04,'--','color','b','linew',2)
plot(1:366, yp_w_tn_04 - (yp_w_no3_04 + yp_w_nh4_04),'--','color','k','linew',2)
% plot(1:366, yp_w_det_04,'--','color','g','linew',2)
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
plot(1:366, yp_w_nh4_04,'--','color','c','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('Concentration (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial.','fontsize',16)
grid on;
% legend('tn','diff(det.)','det','no3','nh4');
legend('tn','diff(det.)','no3','nh4');
set(gca,'fontsize',16,'fontweight','bold')
xlim([1 366])
ylim([0 inf])

figure;
hold on;
plot(1:366, yp_w_tn_af,'--','color','b','linew',2)
plot(1:366, yp_w_tn_af - (yp_w_no3_af + yp_w_nh4_af),'--','color','k','linew',2)
% plot(1:366, yp_w_det_af,'--','color','g','linew',2)
plot(1:366, yp_w_no3_af,'--','color','r','linew',2)
plot(1:366, yp_w_nh4_af,'--','color','c','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('Concentration (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial.','fontsize',16)
grid on;
% legend('tn','diff(det.)','det','no3','nh4');
legend('tn','diff(det.)','no3','nh4');
set(gca,'fontsize',16,'fontweight','bold')
xlim([1 366])



%% DO
figure;
scatter(1:366,reg_clim_do_te);
hold on
plot(1:366, yp_w_do_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('DO (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial DO.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([2 18]) 
xlim([1 366])

%% chl
figure;
scatter(1:366,reg_clim_chl_te);
hold on
plot(1:366, yp_w_chl_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('Chla (mg/m^3)','fontsize',16)
title('sumjin(songjung)-polynomial chl.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 40])
xlim([1 366])


%% no3
figure;
scatter(1:366,reg_clim_no3_te);
hold on
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('NO3 (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial no3.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 3])
xlim([1 366])

%% nh4
figure;
scatter(1:366,reg_clim_nh4_te);
hold on
plot(1:366, yp_w_nh4_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('NH4 (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial nh4.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 2.5])
xlim([1 366])

%% ALL
figure;
hold on;
plot(1:366, yp_w_do_04,'--','color','b','linew',2)
plot(1:366, yp_w_chl_04./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
plot(1:366, yp_w_nh4_04,'--','color','c','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('Concentration (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial.','fontsize',16)
grid on;
legend('DO','Chla','no3','nh4');
set(gca,'fontsize',16)
xlim([1 366])


%% 2004~2019

%% DO
figure;
scatter(1:366,reg_af_clim_do_te);
hold on
plot(1:366, yp_w_do_af,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('DO (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial DO.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([2 18])
xlim([1 366])

%% chl
figure;
scatter(1:366,reg_af_clim_chl_te);
hold on
plot(1:366, yp_w_chl_af,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('Chla (mg/m^3)','fontsize',16)
title('sumjin(songjung)-polynomial chl.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 40])
xlim([1 366])


%% no3
figure;
scatter(1:366,reg_af_clim_no3_te);
hold on
plot(1:366, yp_w_no3_af,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('NO3 (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial no3.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 3])
xlim([1 366])

%% nh4
figure;
scatter(1:366,reg_af_clim_nh4_te);
hold on
plot(1:366, yp_w_nh4_af,'--','color','r','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('NH4 (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial nh4.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 2.5])
xlim([1 366])

%% ALL
figure;
hold on;
plot(1:366, yp_w_do_af,'--','color','b','linew',2)
plot(1:366, yp_w_chl_af./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3_af,'--','color','r','linew',2)
plot(1:366, yp_w_nh4_af,'--','color','c','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('Concentration (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial.','fontsize',16)
grid on;
legend('DO','Chla','no3','nh4');
set(gca,'fontsize',16)
xlim([1 366])


save('sumjin(songjung)_polynomial_climate_to2004_3sig_det.mat'); 

%% vs.
load('songjung_polynomial_climate.mat','yp_w_*')


% DO
figure;
plot(1:366, yp_w_do_04,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_do_af,'--','color','r','linew',2)
plot(1:366, yp_w_do,'--','color','k','linew',2)
xlabel('time(days)','fontsize',16)
ylabel('DO (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial DO.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([6 17])
xlim([1 366])
legend('1st','2nd','whole');

% CHL
figure;
plot(1:366, yp_w_chl_04,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_chl_af,'--','color','r','linew',2)
plot(1:366, yp_w_chl,'--','color','k','linew',2)
legend('1st','2nd','whole');
xlabel('time(days)','fontsize',16)
ylabel('Chla (mg/m^3)','fontsize',16)
title('sumjin(songjung)-polynomial chl.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 20])
xlim([1 366])

% no3
figure;
plot(1:366, yp_w_no3_04,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_no3_af,'--','color','r','linew',2)
plot(1:366, yp_w_no3,'--','color','k','linew',2)
legend('1st','2nd','whole');
xlabel('time(days)','fontsize',16)
ylabel('NO3 (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial no3.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0.4 2.2])
xlim([1 366])

% nh4
figure;
plot(1:366, yp_w_nh4_04,'--','color','b','linew',2)
hold on
plot(1:366, yp_w_nh4_af,'--','color','r','linew',2)
plot(1:366, yp_w_nh4,'--','color','k','linew',2)
legend('1st','2nd','whole');
xlabel('time(days)','fontsize',16)
ylabel('NH4 (mg/L)','fontsize',16)
title('sumjin(songjung)-polynomial nh4.','fontsize',16)
grid on
set(gca,'fontsize',16)
ylim([0 0.5])
xlim([1 366])

txt

w_tn = raw_tn;
% t_tick_all = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
plot(w_tn,'.');
 xlabel('time(year)','fontsize',16)
ylabel('TN (mg/L)','fontsize',16)
grid on; set(gca,'fontsize',16)
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',t_tick_all(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_tn_04=w_tn(1:t_tick(end));
idx1 = find(isnan(w_tn_04) == 0);
w_tn_042=w_tn_04;
upper_bc(sig) = mean(w_tn_04(idx1)) + sig*std(w_tn_04(idx1));
lower_bc(sig) = mean(w_tn_04(idx1)) - sig*std(w_tn_04(idx1));
mean(w_tn_04(idx1))
w_tn_042(find(w_tn_04 > mean(w_tn_04(idx1)) + sig*std(w_tn_04(idx1))))=NaN;
w_tn_042(find(w_tn_04 < mean(w_tn_04(idx1)) - sig*std(w_tn_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_tn_042(~isnan(w_tn_042));
idx = find(isnan(w_tn_042) == 0);
nanaxisx=find(isnan(w_tn_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_tn_04)].*b1(sig) + b0(sig);
% plot(yCalc2{sig},'-.','color','k','linew',2)
% plot([1:1:length(w_tn_04)], upper_bc(sig).*ones(length(w_tn_04)),'color', color_spec(sig),'linew',2)
% plot([1:1:length(w_tn_04)], lower_bc(sig).*ones(length(w_tn_04)),'color', color_spec(sig),'linew',2)
line([1 length(w_tn_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([1 length(w_tn_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
% plot([1:1:length(w_tn_04)],mean(w_tn_04(idx1)).*ones(length(w_tn_04)),'-.','color',color_spec(sig-1),'linew',2);
line([1 length(w_tn_04)],[mean(w_tn_04(idx1)) mean(w_tn_04(idx1))],'color',color_spec(sig-1),'linew',2);
%     for i = 1:length(yCalc3)
%             w_tn_recon(i) = yCalc3(i)  + yp_w_tn(indx_366{i});
%     end
% plot([1:length(w_tn)], w_tn_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_tn_04=w_tn(t_tick(end)+1:end);
idx1 = find(isnan(w_tn_04) == 0);
w_tn_042=w_tn_04;
upper_bc(sig) = mean(w_tn_04(idx1)) + sig*std(w_tn_04(idx1));
lower_bc(sig) = mean(w_tn_04(idx1)) - sig*std(w_tn_04(idx1));
mean(w_tn_04(idx1))
w_tn_042(find(w_tn_04 > mean(w_tn_04(idx1)) + sig*std(w_tn_04(idx1))))=NaN;
w_tn_042(find(w_tn_04 < mean(w_tn_04(idx1)) - sig*std(w_tn_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_tn_042(~isnan(w_tn_042));
idx = find(isnan(w_tn_042) == 0);
nanaxisx=find(isnan(w_tn_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_tn_04)].*b1(sig) + b0(sig);
% plot([t_tick(end)+1:1:t_tick(end)+length(w_tn_04)],yCalc2{sig},'-.','color','k','linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_tn_04)], upper_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_tn_04))),'color', color_spec(sig),'linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_tn_04)], lower_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_tn_04))),'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_tn_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_tn_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
line([t_tick(end)+1 t_tick(end)+length(w_tn_04)],[mean(w_tn_04(idx1)) mean(w_tn_04(idx1))],'color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_tn_recon(i) = yCalc3(i)  + yp_w_tn(indx_366{i});
%     end
% plot([1:length(w_tn)], w_tn_recon,'-.','color',color_spec(sig),'linew',1)
end
title(['sumjin(songjung) TN raw  ' num2str(sig) '-sigma'],'fontsize',16)
ylim([0 inf])
set(gca,'xtick',[t_tick_all(1:4:end)]);
set(gca,'xlim',[1 t_tick_all(end)]);
set(gca,'xticklabel',1989:4:2019);



w_do = raw_do;
% t_tick_all = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
plot(w_do,'.');
 xlabel('time(year)','fontsize',16)
ylabel('do (mg/L)','fontsize',16)
grid on; set(gca,'fontsize',16)
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',t_tick_all(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_do_04=w_do(1:t_tick(end));
idx1 = find(isnan(w_do_04) == 0);
w_do_042=w_do_04;
upper_bc(sig) = mean(w_do_04(idx1)) + sig*std(w_do_04(idx1));
lower_bc(sig) = mean(w_do_04(idx1)) - sig*std(w_do_04(idx1));
mean(w_do_04(idx1))
w_do_042(find(w_do_04 > mean(w_do_04(idx1)) + sig*std(w_do_04(idx1))))=NaN;
w_do_042(find(w_do_04 < mean(w_do_04(idx1)) - sig*std(w_do_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_do_042(~isnan(w_do_042));
idx = find(isnan(w_do_042) == 0);
nanaxisx=find(isnan(w_do_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_do_04)].*b1(sig) + b0(sig);
% plot(yCalc2{sig},'-.','color','k','linew',2)
% plot([1:1:length(w_do_04)], upper_bc(sig).*ones(length(w_do_04)),'color', color_spec(sig),'linew',2)
% plot([1:1:length(w_do_04)], lower_bc(sig).*ones(length(w_do_04)),'color', color_spec(sig),'linew',2)
line([1 length(w_do_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([1 length(w_do_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
% plot([1:1:length(w_do_04)],mean(w_do_04(idx1)).*ones(length(w_do_04)),'-.','color',color_spec(sig-1),'linew',2);
line([1 length(w_do_04)],[mean(w_do_04(idx1)) mean(w_do_04(idx1))],'color',color_spec(sig-1),'linew',2);
%     for i = 1:length(yCalc3)
%             w_do_recon(i) = yCalc3(i)  + yp_w_do(indx_366{i});
%     end
% plot([1:length(w_do)], w_do_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_do_04=w_do(t_tick(end)+1:end);
idx1 = find(isnan(w_do_04) == 0);
w_do_042=w_do_04;
upper_bc(sig) = mean(w_do_04(idx1)) + sig*std(w_do_04(idx1));
lower_bc(sig) = mean(w_do_04(idx1)) - sig*std(w_do_04(idx1));
mean(w_do_04(idx1))
w_do_042(find(w_do_04 > mean(w_do_04(idx1)) + sig*std(w_do_04(idx1))))=NaN;
w_do_042(find(w_do_04 < mean(w_do_04(idx1)) - sig*std(w_do_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_do_042(~isnan(w_do_042));
idx = find(isnan(w_do_042) == 0);
nanaxisx=find(isnan(w_do_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_do_04)].*b1(sig) + b0(sig);
% plot([t_tick(end)+1:1:t_tick(end)+length(w_do_04)],yCalc2{sig},'-.','color','k','linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_do_04)], upper_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_do_04))),'color', color_spec(sig),'linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_do_04)], lower_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_do_04))),'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_do_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_do_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
line([t_tick(end)+1 t_tick(end)+length(w_do_04)],[mean(w_do_04(idx1)) mean(w_do_04(idx1))],'color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_do_recon(i) = yCalc3(i)  + yp_w_do(indx_366{i});
%     end
% plot([1:length(w_do)], w_do_recon,'-.','color',color_spec(sig),'linew',1)
end
title(['sumjin(songjung) do raw  ' num2str(sig) '-sigma'],'fontsize',16)
ylim([0 inf])
set(gca,'xtick',[t_tick_all(1:4:end)]);
set(gca,'xlim',[1 t_tick_all(end)]);
set(gca,'xticklabel',1989:4:2019);


w_chl = raw_chl;
% t_tick_all = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
plot(w_chl,'.');
 xlabel('time(year)','fontsize',16)
ylabel('chl (mg/L)','fontsize',16)
grid on; set(gca,'fontsize',16)
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',t_tick_all(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_chl_04=w_chl(1:t_tick(end));
idx1 = find(isnan(w_chl_04) == 0);
w_chl_042=w_chl_04;
upper_bc(sig) = mean(w_chl_04(idx1)) + sig*std(w_chl_04(idx1));
lower_bc(sig) = mean(w_chl_04(idx1)) - sig*std(w_chl_04(idx1));
mean(w_chl_04(idx1))
w_chl_042(find(w_chl_04 > mean(w_chl_04(idx1)) + sig*std(w_chl_04(idx1))))=NaN;
w_chl_042(find(w_chl_04 < mean(w_chl_04(idx1)) - sig*std(w_chl_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_chl_042(~isnan(w_chl_042));
idx = find(isnan(w_chl_042) == 0);
nanaxisx=find(isnan(w_chl_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_chl_04)].*b1(sig) + b0(sig);
% plot(yCalc2{sig},'-.','color','k','linew',2)
% plot([1:1:length(w_chl_04)], upper_bc(sig).*ones(length(w_chl_04)),'color', color_spec(sig),'linew',2)
% plot([1:1:length(w_chl_04)], lower_bc(sig).*ones(length(w_chl_04)),'color', color_spec(sig),'linew',2)
line([1 length(w_chl_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([1 length(w_chl_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
% plot([1:1:length(w_chl_04)],mean(w_chl_04(idx1)).*ones(length(w_chl_04)),'-.','color',color_spec(sig-1),'linew',2);
line([1 length(w_chl_04)],[mean(w_chl_04(idx1)) mean(w_chl_04(idx1))],'color',color_spec(sig-1),'linew',2);
%     for i = 1:length(yCalc3)
%             w_chl_recon(i) = yCalc3(i)  + yp_w_chl(indx_366{i});
%     end
% plot([1:length(w_chl)], w_chl_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_chl_04=w_chl(t_tick(end)+1:end);
idx1 = find(isnan(w_chl_04) == 0);
w_chl_042=w_chl_04;
upper_bc(sig) = mean(w_chl_04(idx1)) + sig*std(w_chl_04(idx1));
lower_bc(sig) = mean(w_chl_04(idx1)) - sig*std(w_chl_04(idx1));
mean(w_chl_04(idx1))
w_chl_042(find(w_chl_04 > mean(w_chl_04(idx1)) + sig*std(w_chl_04(idx1))))=NaN;
w_chl_042(find(w_chl_04 < mean(w_chl_04(idx1)) - sig*std(w_chl_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_chl_042(~isnan(w_chl_042));
idx = find(isnan(w_chl_042) == 0);
nanaxisx=find(isnan(w_chl_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_chl_04)].*b1(sig) + b0(sig);
% plot([t_tick(end)+1:1:t_tick(end)+length(w_chl_04)],yCalc2{sig},'-.','color','k','linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_chl_04)], upper_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_chl_04))),'color', color_spec(sig),'linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_chl_04)], lower_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_chl_04))),'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_chl_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_chl_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
line([t_tick(end)+1 t_tick(end)+length(w_chl_04)],[mean(w_chl_04(idx1)) mean(w_chl_04(idx1))],'color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_chl_recon(i) = yCalc3(i)  + yp_w_chl(indx_366{i});
%     end
% plot([1:length(w_chl)], w_chl_recon,'-.','color',color_spec(sig),'linew',1)
end
title(['sumjin(songjung) chl raw  ' num2str(sig) '-sigma'],'fontsize',16)
ylim([0 inf])
set(gca,'xtick',[t_tick_all(1:3:end)]);
set(gca,'xticklabel',1989:3:2019);
set(gca,'xlim',[t_tick_all(13) t_tick_all(end)]);
ylim([0 50])


w_nh4 = raw_nh4;
% t_tick_all = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
plot(w_nh4,'.');
 xlabel('time(year)','fontsize',16)
ylabel('nh4 (mg/L)','fontsize',16)
grid on; set(gca,'fontsize',16)
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',t_tick_all(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_nh4_04=w_nh4(1:t_tick(end));
idx1 = find(isnan(w_nh4_04) == 0);
w_nh4_042=w_nh4_04;
upper_bc(sig) = mean(w_nh4_04(idx1)) + sig*std(w_nh4_04(idx1));
lower_bc(sig) = mean(w_nh4_04(idx1)) - sig*std(w_nh4_04(idx1));
mean(w_nh4_04(idx1))
w_nh4_042(find(w_nh4_04 > mean(w_nh4_04(idx1)) + sig*std(w_nh4_04(idx1))))=NaN;
w_nh4_042(find(w_nh4_04 < mean(w_nh4_04(idx1)) - sig*std(w_nh4_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_nh4_042(~isnan(w_nh4_042));
idx = find(isnan(w_nh4_042) == 0);
nanaxisx=find(isnan(w_nh4_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_nh4_04)].*b1(sig) + b0(sig);
% plot(yCalc2{sig},'-.','color','k','linew',2)
% plot([1:1:length(w_nh4_04)], upper_bc(sig).*ones(length(w_nh4_04)),'color', color_spec(sig),'linew',2)
% plot([1:1:length(w_nh4_04)], lower_bc(sig).*ones(length(w_nh4_04)),'color', color_spec(sig),'linew',2)
line([1 length(w_nh4_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([1 length(w_nh4_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
% plot([1:1:length(w_nh4_04)],mean(w_nh4_04(idx1)).*ones(length(w_nh4_04)),'-.','color',color_spec(sig-1),'linew',2);
line([1 length(w_nh4_04)],[mean(w_nh4_04(idx1)) mean(w_nh4_04(idx1))],'color',color_spec(sig-1),'linew',2);
%     for i = 1:length(yCalc3)
%             w_nh4_recon(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
%     end
% plot([1:length(w_nh4)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_nh4_04=w_nh4(t_tick(end)+1:end);
idx1 = find(isnan(w_nh4_04) == 0);
w_nh4_042=w_nh4_04;
upper_bc(sig) = mean(w_nh4_04(idx1)) + sig*std(w_nh4_04(idx1));
lower_bc(sig) = mean(w_nh4_04(idx1)) - sig*std(w_nh4_04(idx1));
mean(w_nh4_04(idx1))
w_nh4_042(find(w_nh4_04 > mean(w_nh4_04(idx1)) + sig*std(w_nh4_04(idx1))))=NaN;
w_nh4_042(find(w_nh4_04 < mean(w_nh4_04(idx1)) - sig*std(w_nh4_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_nh4_042(~isnan(w_nh4_042));
idx = find(isnan(w_nh4_042) == 0);
nanaxisx=find(isnan(w_nh4_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_nh4_04)].*b1(sig) + b0(sig);
% plot([t_tick(end)+1:1:t_tick(end)+length(w_nh4_04)],yCalc2{sig},'-.','color','k','linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_nh4_04)], upper_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_nh4_04))),'color', color_spec(sig),'linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_nh4_04)], lower_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_nh4_04))),'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_nh4_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_nh4_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
line([t_tick(end)+1 t_tick(end)+length(w_nh4_04)],[mean(w_nh4_04(idx1)) mean(w_nh4_04(idx1))],'color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_nh4_recon(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
%     end
% plot([1:length(w_nh4)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
end
title(['sumjin(songjung) nh4 raw  ' num2str(sig) '-sigma'],'fontsize',16)
ylim([0 inf])
set(gca,'xtick',[t_tick_all(1:4:end)]);
set(gca,'xlim',[1 t_tick_all(end)]);
set(gca,'xticklabel',1989:4:2019);


w_nh4 = raw_nh4;
% t_tick_all = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
plot(w_nh4,'.');
 xlabel('time(year)','fontsize',16)
ylabel('nh4 (mg/L)','fontsize',16)
grid on; set(gca,'fontsize',16)
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',t_tick_all(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_nh4_04=w_nh4(1:t_tick(end));
idx1 = find(isnan(w_nh4_04) == 0);
w_nh4_042=w_nh4_04;
upper_bc(sig) = mean(w_nh4_04(idx1)) + sig*std(w_nh4_04(idx1));
lower_bc(sig) = mean(w_nh4_04(idx1)) - sig*std(w_nh4_04(idx1));
mean(w_nh4_04(idx1))
w_nh4_042(find(w_nh4_04 > mean(w_nh4_04(idx1)) + sig*std(w_nh4_04(idx1))))=NaN;
w_nh4_042(find(w_nh4_04 < mean(w_nh4_04(idx1)) - sig*std(w_nh4_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_nh4_042(~isnan(w_nh4_042));
idx = find(isnan(w_nh4_042) == 0);
nanaxisx=find(isnan(w_nh4_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_nh4_04)].*b1(sig) + b0(sig);
% plot(yCalc2{sig},'-.','color','k','linew',2)
% plot([1:1:length(w_nh4_04)], upper_bc(sig).*ones(length(w_nh4_04)),'color', color_spec(sig),'linew',2)
% plot([1:1:length(w_nh4_04)], lower_bc(sig).*ones(length(w_nh4_04)),'color', color_spec(sig),'linew',2)
line([1 length(w_nh4_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([1 length(w_nh4_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
% plot([1:1:length(w_nh4_04)],mean(w_nh4_04(idx1)).*ones(length(w_nh4_04)),'-.','color',color_spec(sig-1),'linew',2);
line([1 length(w_nh4_04)],[mean(w_nh4_04(idx1)) mean(w_nh4_04(idx1))],'color',color_spec(sig-1),'linew',2);
%     for i = 1:length(yCalc3)
%             w_nh4_recon(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
%     end
% plot([1:length(w_nh4)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_nh4_04=w_nh4(t_tick(end)+1:end);
idx1 = find(isnan(w_nh4_04) == 0);
w_nh4_042=w_nh4_04;
upper_bc(sig) = mean(w_nh4_04(idx1)) + sig*std(w_nh4_04(idx1));
lower_bc(sig) = mean(w_nh4_04(idx1)) - sig*std(w_nh4_04(idx1));
mean(w_nh4_04(idx1))
w_nh4_042(find(w_nh4_04 > mean(w_nh4_04(idx1)) + sig*std(w_nh4_04(idx1))))=NaN;
w_nh4_042(find(w_nh4_04 < mean(w_nh4_04(idx1)) - sig*std(w_nh4_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_nh4_042(~isnan(w_nh4_042));
idx = find(isnan(w_nh4_042) == 0);
nanaxisx=find(isnan(w_nh4_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_nh4_04)].*b1(sig) + b0(sig);
% plot([t_tick(end)+1:1:t_tick(end)+length(w_nh4_04)],yCalc2{sig},'-.','color','k','linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_nh4_04)], upper_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_nh4_04))),'color', color_spec(sig),'linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_nh4_04)], lower_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_nh4_04))),'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_nh4_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_nh4_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
line([t_tick(end)+1 t_tick(end)+length(w_nh4_04)],[mean(w_nh4_04(idx1)) mean(w_nh4_04(idx1))],'color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_nh4_recon(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
%     end
% plot([1:length(w_nh4)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
end
title(['sumjin(songjung) nh4 raw  ' num2str(sig) '-sigma'],'fontsize',16)
ylim([0 inf])
set(gca,'xtick',[t_tick_all(1:3:end)]);
set(gca,'xticklabel',1989:3:2019);
set(gca,'xlim',[t_tick_all(7) t_tick_all(end)]);


w_no3 = raw_no3;
% t_tick_all = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
plot(w_no3,'.');
 xlabel('time(year)','fontsize',16)
ylabel('no3 (mg/L)','fontsize',16)
grid on; set(gca,'fontsize',16)
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',t_tick_all(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_no3_04=w_no3(1:t_tick(end));
idx1 = find(isnan(w_no3_04) == 0);
w_no3_042=w_no3_04;
upper_bc(sig) = mean(w_no3_04(idx1)) + sig*std(w_no3_04(idx1));
lower_bc(sig) = mean(w_no3_04(idx1)) - sig*std(w_no3_04(idx1));
mean(w_no3_04(idx1))
w_no3_042(find(w_no3_04 > mean(w_no3_04(idx1)) + sig*std(w_no3_04(idx1))))=NaN;
w_no3_042(find(w_no3_04 < mean(w_no3_04(idx1)) - sig*std(w_no3_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_no3_042(~isnan(w_no3_042));
idx = find(isnan(w_no3_042) == 0);
nanaxisx=find(isnan(w_no3_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_no3_04)].*b1(sig) + b0(sig);
% plot(yCalc2{sig},'-.','color','k','linew',2)
% plot([1:1:length(w_no3_04)], upper_bc(sig).*ones(length(w_no3_04)),'color', color_spec(sig),'linew',2)
% plot([1:1:length(w_no3_04)], lower_bc(sig).*ones(length(w_no3_04)),'color', color_spec(sig),'linew',2)
line([1 length(w_no3_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([1 length(w_no3_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
% plot([1:1:length(w_no3_04)],mean(w_no3_04(idx1)).*ones(length(w_no3_04)),'-.','color',color_spec(sig-1),'linew',2);
line([1 length(w_no3_04)],[mean(w_no3_04(idx1)) mean(w_no3_04(idx1))],'color',color_spec(sig-1),'linew',2);
%     for i = 1:length(yCalc3)
%             w_no3_recon(i) = yCalc3(i)  + yp_w_no3(indx_366{i});
%     end
% plot([1:length(w_no3)], w_no3_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_no3_04=w_no3(t_tick(end)+1:end);
idx1 = find(isnan(w_no3_04) == 0);
w_no3_042=w_no3_04;
upper_bc(sig) = mean(w_no3_04(idx1)) + sig*std(w_no3_04(idx1));
lower_bc(sig) = mean(w_no3_04(idx1)) - sig*std(w_no3_04(idx1));
mean(w_no3_04(idx1))
w_no3_042(find(w_no3_04 > mean(w_no3_04(idx1)) + sig*std(w_no3_04(idx1))))=NaN;
w_no3_042(find(w_no3_04 < mean(w_no3_04(idx1)) - sig*std(w_no3_04(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_no3_042(~isnan(w_no3_042));
idx = find(isnan(w_no3_042) == 0);
nanaxisx=find(isnan(w_no3_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_no3_04)].*b1(sig) + b0(sig);
% plot([t_tick(end)+1:1:t_tick(end)+length(w_no3_04)],yCalc2{sig},'-.','color','k','linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_no3_04)], upper_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_no3_04))),'color', color_spec(sig),'linew',2)
% plot([t_tick_all(25)+1:1:t_tick_all(25)+length(w_no3_04)], lower_bc(sig).*ones(length(t_tick_all(25)+1:1:t_tick_all(25)+length(w_no3_04))),'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_no3_04)], [upper_bc(sig) upper_bc(sig)],'color', color_spec(sig),'linew',2)
line([t_tick(end)+1 t_tick(end)+length(w_no3_04)], [lower_bc(sig) lower_bc(sig)],'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
line([t_tick(end)+1 t_tick(end)+length(w_no3_04)],[mean(w_no3_04(idx1)) mean(w_no3_04(idx1))],'color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_no3_recon(i) = yCalc3(i)  + yp_w_no3(indx_366{i});
%     end
% plot([1:length(w_no3)], w_no3_recon,'-.','color',color_spec(sig),'linew',1)
end
title(['sumjin(songjung) no3 raw  ' num2str(sig) '-sigma'],'fontsize',16)
ylim([0 inf])
set(gca,'xtick',[t_tick_all(1:3:end)]);
set(gca,'xticklabel',1989:3:2019);
set(gca,'xlim',[t_tick_all(7) t_tick_all(end)]);


