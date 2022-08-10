close all; clear; clc; 
cd E:\장기생태\Dynamic\06_river\환경과학원
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

return

for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
end


% make 366 mm-dd
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

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, r_date_txt_c)] == 1)     
end

%confirm how may data on each the days
for i = 1:366; size_c(i)=length(indx{i}); end
plot(size_c)

% make quasi_climate
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do(i) = NaN;
        w_chl(i) = NaN; 
        w_no3(i) = NaN; 
        w_nh4(i) = NaN; 
    else
        w_do(i) = nanmean(r_do(indx{i}));
        w_chl(i) = nanmean(r_chl(indx{i}));
        w_no3(i) = nanmean(r_no3(indx{i}));
        w_nh4(i) = nanmean(r_nh4(indx{i}));
    end
end

% make quasi_climate scatter
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do_c{i} = NaN;
        w_chl_c{i} = NaN; 
        w_no3_c{i} = NaN; 
        w_nh4_c{i} = NaN; 
    else
        w_do_c{i} = r_do(indx{i});
        w_chl_c{i} = r_chl(indx{i});
        w_no3_c{i} = r_no3(indx{i});
        w_nh4_c{i} = r_nh4(indx{i});
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
for i = 1:17
    if i ==1
         t_tick(i) = t_tick_pre(i);
    else
        t_tick(i) = sum(t_tick_pre(1:i));
    end
end

return

figure; hold on ;for i =1:366; plot(i,w_do_c{i},'bo'); end
figure; hold on ;for i =1:366; plot(i,w_chl_c{i},'bo'); end
figure; hold on ;for i =1:366; plot(i,w_nh4_c{i},'bo'); end

xp = 1:366;
xp_w_do = find(isnan(w_do)==0);
pf_w_do = polyfit(xp_w_do,w_do(xp_w_do),3);
yp_w_do = polyval(pf_w_do,xp);


xp = 1:366;
xp_w_chl = find(isnan(w_chl)==0);
pf_w_chl = polyfit(xp_w_chl,w_chl(xp_w_chl),3);
yp_w_chl = polyval(pf_w_chl,xp);


xp = 1:366;
xp_w_no3 = find(isnan(w_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,w_no3(xp_w_no3),3);
yp_w_no3 = polyval(pf_w_no3,xp);


xp = 1:366;
xp_w_nh4 = find(isnan(w_nh4)==0);
pf_w_nh4 = polyfit(xp_w_nh4,w_nh4(xp_w_nh4),3);
yp_w_nh4 = polyval(pf_w_nh4,xp);


w_do(i) = nanmean(r_do(indx{i}));
        w_chl(i) = nanmean(r_chl(indx{i}));
        w_no3(i) = nanmean(r_no3(indx{i}));
        w_nh4(i) = nanmean(r_nh4(indx{i}));

%% DO
figure;
scatter(1:366,w_do);
hold on
plot(1:366, yp_w_do,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('sachoen-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% chl
figure;
scatter(1:366,w_chl);
hold on
plot(1:366, yp_w_chl,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title('sachoen-polynomial chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])


%% no3
figure;
scatter(1:366,w_no3);
hold on
plot(1:366, yp_w_no3,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title('sachoen-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% nh4
figure;
scatter(1:366,w_nh4);
hold on
plot(1:366, yp_w_nh4,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title('sachoen-polynomial nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% ALL
figure;
hold on;
plot(1:366, yp_w_do,'--','color','b','linew',2)
plot(1:366, yp_w_chl./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3,'--','color','r','linew',2)
plot(1:366, yp_w_nh4,'--','color','c','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Concentration (mg/L)','fontsize',13)
title('sachoen-polynomial.','fontsize',13)
grid on;
legend('DO','Chla','no3','nh4');
set(gca,'fontsize',13)
xlim([1 366])
save('sumjin(songjung)_polynomial_climate_to2004.mat'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd E:\장기생태\Dynamic\06_river\환경과학원
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
end


% make 366 mm-dd
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

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, r_date_txt_c)] == 1)     
end

%confirm how may data on each the days
for i = 1:366; size_c(i)=length(indx{i}); end
plot(size_c)

% make quasi_climate
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do(i) = NaN;
        w_chl(i) = NaN; 
        w_no3(i) = NaN; 
        w_nh4(i) = NaN; 
    else
        w_do(i) = nanmean(r_do(indx{i}));
        w_chl(i) = nanmean(r_chl(indx{i}));
        w_no3(i) = nanmean(r_no3(indx{i}));
        w_nh4(i) = nanmean(r_nh4(indx{i}));
    end
end

% make quasi_climate scatter
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do_c{i} = NaN;
        w_chl_c{i} = NaN; 
        w_no3_c{i} = NaN; 
        w_nh4_c{i} = NaN; 
    else
        w_do_c{i} = r_do(indx{i});
        w_chl_c{i} = r_chl(indx{i});
        w_no3_c{i} = r_no3(indx{i});
        w_nh4_c{i} = r_nh4(indx{i});
    end
end

figure; hold on ;for i =1:366; plot(i,w_do_c{i},'bo'); end
figure; hold on ;for i =1:366; plot(i,w_chl_c{i},'bo'); end
figure; hold on ;for i =1:366; plot(i,w_nh4_c{i},'bo'); end

xp = 1:366;
xp_w_do = find(isnan(w_do)==0);
pf_w_do = polyfit(xp_w_do,w_do(xp_w_do),3);
yp_w_do = polyval(pf_w_do,xp);


xp = 1:366;
xp_w_chl = find(isnan(w_chl)==0);
pf_w_chl = polyfit(xp_w_chl,w_chl(xp_w_chl),3);
yp_w_chl = polyval(pf_w_chl,xp);


xp = 1:366;
xp_w_no3 = find(isnan(w_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,w_no3(xp_w_no3),3);
yp_w_no3 = polyval(pf_w_no3,xp);


xp = 1:366;
xp_w_nh4 = find(isnan(w_nh4)==0);
pf_w_nh4 = polyfit(xp_w_nh4,w_nh4(xp_w_nh4),3);
yp_w_nh4 = polyval(pf_w_nh4,xp);


w_do(i) = nanmean(r_do(indx{i}));
        w_chl(i) = nanmean(r_chl(indx{i}));
        w_no3(i) = nanmean(r_no3(indx{i}));
        w_nh4(i) = nanmean(r_nh4(indx{i}));

%% DO
figure;
scatter(1:366,w_do);
hold on
plot(1:366, yp_w_do,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('Namgang-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% chl
figure;
scatter(1:366,w_chl);
hold on
plot(1:366, yp_w_chl,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title('Namgang-polynomial chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])


%% no3
figure;
scatter(1:366,w_no3);
hold on
plot(1:366, yp_w_no3,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title('Namgang-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% nh4
figure;
scatter(1:366,w_nh4);
hold on
plot(1:366, yp_w_nh4,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title('Namgang-polynomial nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% ALL
figure;
hold on;
plot(1:366, yp_w_do,'--','color','b','linew',2)
plot(1:366, yp_w_chl./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3,'--','color','r','linew',2)
plot(1:366, yp_w_nh4,'--','color','c','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Concentration (mg/L)','fontsize',13)
title('Namgang-polynomial.','fontsize',13)
grid on;
legend('DO','Chla','no3','nh4');
set(gca,'fontsize',13)
xlim([1 366])
save('namgang_polynomial_climate.mat'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  gahwa river %%%%%%%%%%
close all; clear; clc; 
cd E:\장기생태\Dynamic\06_river\환경과학원
% [raw txt]=xlsread('수질측정지점_구례_fix.xls','검색결과','');
% [raw txt]=xlsread('수질측정지점_남강댐1.xls','검색결과','');
[raw txt]=xlsread('수질측정지점_하동.xls','검색결과','');
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
end


% make 366 mm-dd
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

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, r_date_txt_c)] == 1)     
end

%confirm how may data on each the days
for i = 1:366; size_c(i)=length(indx{i}); end
plot(size_c)

% make quasi_climate
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do(i) = NaN;
        w_chl(i) = NaN; 
        w_no3(i) = NaN; 
        w_nh4(i) = NaN; 
    else
        w_do(i) = nanmean(r_do(indx{i}));
        w_chl(i) = nanmean(r_chl(indx{i}));
        w_no3(i) = nanmean(r_no3(indx{i}));
        w_nh4(i) = nanmean(r_nh4(indx{i}));
    end
end

% make quasi_climate scatter
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do_c{i} = NaN;
        w_chl_c{i} = NaN; 
        w_no3_c{i} = NaN; 
        w_nh4_c{i} = NaN; 
    else
        w_do_c{i} = r_do(indx{i});
        w_chl_c{i} = r_chl(indx{i});
        w_no3_c{i} = r_no3(indx{i});
        w_nh4_c{i} = r_nh4(indx{i});
    end
end

figure; hold on ;for i =1:366; plot(i,w_do_c{i},'bo'); end
figure; hold on ;for i =1:366; plot(i,w_chl_c{i},'bo'); end
figure; hold on ;for i =1:366; plot(i,w_nh4_c{i},'bo'); end

xp = 1:366;
xp_w_do = find(isnan(w_do)==0);
pf_w_do = polyfit(xp_w_do,w_do(xp_w_do),3);
yp_w_do = polyval(pf_w_do,xp);


xp = 1:366;
xp_w_chl = find(isnan(w_chl)==0);
pf_w_chl = polyfit(xp_w_chl,w_chl(xp_w_chl),3);
yp_w_chl = polyval(pf_w_chl,xp);


xp = 1:366;
xp_w_no3 = find(isnan(w_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,w_no3(xp_w_no3),3);
yp_w_no3 = polyval(pf_w_no3,xp);


xp = 1:366;
xp_w_nh4 = find(isnan(w_nh4)==0);
pf_w_nh4 = polyfit(xp_w_nh4,w_nh4(xp_w_nh4),3);
yp_w_nh4 = polyval(pf_w_nh4,xp);


w_do(i) = nanmean(r_do(indx{i}));
        w_chl(i) = nanmean(r_chl(indx{i}));
        w_no3(i) = nanmean(r_no3(indx{i}));
        w_nh4(i) = nanmean(r_nh4(indx{i}));

%% DO
figure;
scatter(1:366,w_do);
hold on
plot(1:366, yp_w_do,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('hadong-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% chl
figure;
scatter(1:366,w_chl);
hold on
plot(1:366, yp_w_chl,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title('hadong-polynomial chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])


%% no3
figure;
scatter(1:366,w_no3);
hold on
plot(1:366, yp_w_no3,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title('hadong-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% nh4
figure;
scatter(1:366,w_nh4);
hold on
plot(1:366, yp_w_nh4,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title('hadong-polynomial nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% ALL
figure;
hold on;
plot(1:366, yp_w_do,'--','color','b','linew',2)
plot(1:366, yp_w_chl./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3,'--','color','r','linew',2)
plot(1:366, yp_w_nh4,'--','color','c','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Concentration (mg/L)','fontsize',13)
title('hadong-polynomial.','fontsize',13)
grid on;
legend('DO','Chla','no3','nh4');
set(gca,'fontsize',13)
xlim([1 366])
save('hadong_polynomial_climate.mat'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SUMJIN (songjung)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd E:\장기생태\Dynamic\06_river\환경과학원
[raw txt]=xlsread('수질측정지점_구례_fix.xls','검색결과','');
% [raw txt]=xlsread('수질측정지점_남강댐1.xls','검색결과','');

dash_c = '.';
r_txt_ud = flipud(txt);
% r_txt_ud = txt;
r_date_txt=[char(r_txt_ud(1:end-1,5))];

r_do_txt=[r_txt_ud(1:end-1,9)];
for i = 1:length(r_do_txt)
    if strcmp(r_do_txt{i,1},'') == 1
       r_do(i) = NaN;
    elseif strcmp(r_do_txt{i,1},'') == 0
       r_do(i) = str2num(char(r_do_txt{i,1}));
    end
end

r_chl_txt=[r_txt_ud(1:end-1,16)];
for i = 1:length(r_chl_txt)
    if strcmp(r_chl_txt{i,1},'') == 1
       r_chl(i) = NaN;
    elseif strcmp(r_chl_txt{i,1},'') == 0
       r_chl(i) = str2num(char(r_chl_txt{i,1}));
    end
end


r_no3_txt=[r_txt_ud(1:end-1,20)];
for i = 1:length(r_no3_txt)
    if strcmp(r_no3_txt{i,1},'') == 1
       r_no3(i) = NaN;
    elseif strcmp(r_no3_txt{i,1},'') == 0
       r_no3(i) = str2num(char(r_no3_txt{i,1}));
    end
end

r_nh4_txt=[r_txt_ud(1:end-1,21)];
for i = 1:length(r_nh4_txt)
    if strcmp(r_nh4_txt{i,1},'') == 1
       r_nh4(i) = NaN;
    elseif strcmp(r_nh4_txt{i,1},'') == 0
       r_nh4(i) = str2num(char(r_nh4_txt{i,1}));
    end
end


for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
end

% make 366 mm-dd
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

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, r_date_txt_c)] == 1)     
end

%confirm how may data on each the days
for i = 1:366; size_c(i)=length(indx{i}); end
plot(size_c)

% make quasi_climate
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do(i) = NaN;
        w_chl(i) = NaN; 
        w_no3(i) = NaN; 
        w_nh4(i) = NaN; 
    else
        w_do(i) = nanmean(r_do(indx{i}));
        w_chl(i) = nanmean(r_chl(indx{i}));
        w_no3(i) = nanmean(r_no3(indx{i}));
        w_nh4(i) = nanmean(r_nh4(indx{i}));
    end
end

% make quasi_climate scatter
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do_c{i} = NaN;
        w_chl_c{i} = NaN; 
        w_no3_c{i} = NaN; 
        w_nh4_c{i} = NaN; 
    else
        w_do_c{i} = r_do(indx{i});
        w_chl_c{i} = r_chl(indx{i});
        w_no3_c{i} = r_no3(indx{i});
        w_nh4_c{i} = r_nh4(indx{i});
    end
end

% figure; hold on ;for i =1:366; plot(i,w_do_c{i},'bo'); end
% figure; hold on ;for i =1:366; plot(i,w_chl_c{i},'bo'); end
% figure; hold on ;for i =1:366; plot(i,w_nh4_c{i},'bo'); end

xp = 1:366;
xp_w_do = find(isnan(w_do)==0);
pf_w_do = polyfit(xp_w_do,w_do(xp_w_do),3);
yp_w_do = polyval(pf_w_do,xp);


xp = 1:366;
xp_w_chl = find(isnan(w_chl)==0);
pf_w_chl = polyfit(xp_w_chl,w_chl(xp_w_chl),3);
yp_w_chl = polyval(pf_w_chl,xp);


xp = 1:366;
xp_w_no3 = find(isnan(w_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,w_no3(xp_w_no3),3);
yp_w_no3 = polyval(pf_w_no3,xp);


xp = 1:366;
xp_w_nh4 = find(isnan(w_nh4)==0);
pf_w_nh4 = polyfit(xp_w_nh4,w_nh4(xp_w_nh4),3);
yp_w_nh4 = polyval(pf_w_nh4,xp);


w_do(i) = nanmean(r_do(indx{i}));
        w_chl(i) = nanmean(r_chl(indx{i}));
        w_no3(i) = nanmean(r_no3(indx{i}));
        w_nh4(i) = nanmean(r_nh4(indx{i}));

%% DO
figure;
scatter(1:366,w_do);
hold on
plot(1:366, yp_w_do,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% chl
figure;
scatter(1:366,w_chl);
hold on
plot(1:366, yp_w_chl,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title('sumjin(songjung)-polynomial chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])


%% no3
figure;
scatter(1:366,w_no3);
hold on
plot(1:366, yp_w_no3,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% nh4
figure;
scatter(1:366,w_nh4);
hold on
plot(1:366, yp_w_nh4,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])

%% ALL
figure;
hold on;
plot(1:366, yp_w_do,'--','color','b','linew',2)
plot(1:366, yp_w_chl./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3,'--','color','r','linew',2)
plot(1:366, yp_w_nh4,'--','color','c','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Concentration (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial.','fontsize',13)
grid on;
legend('DO','Chla','no3','nh4');
set(gca,'fontsize',13)
xlim([1 366])

save('songjung_polynomial_climate.mat'); 


figure; hold on ;for i =1:366; plot(i,w_do_c{i},'bo'); end
plot(1:366, yp_w_do,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])


figure; hold on ;for i =1:366; plot(i,w_chl_c{i},'bo'); end
plot(1:366, yp_w_chl,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title('sumjin(songjung)-polynomial chl.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])


figure; hold on ;for i =1:366; plot(i,w_no3_c{i},'bo'); end
plot(1:366, yp_w_no3,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])


figure; hold on ;for i =1:366; plot(i,w_nh4_c{i},'bo'); end
plot(1:366, yp_w_nh4,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial nh4.','fontsize',13)
grid on
set(gca,'fontsize',13)
xlim([1 366])
