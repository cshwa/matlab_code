%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  gahwa river %%%%%%%%%%
close all; clear; clc; 
cd E:\장기생태\Dynamic\06_river\환경과학원
% [raw txt]=xlsread('수질측정지점_구례_fix.xls','검색결과','');
[raw txt]=xlsread('수질측정지점_사천천.xls','검색결과','');

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
    r_date_txt_c{i,1} = r_date_txt(i,:); % raw
    r_yymm_c{i,1} = r_date_txt(i,1:end-3); % remove days 
    r_yy_c{i,1} = r_date_txt(i,1:4); % remove days 
    r_mmdd_c{i,1} = r_date_txt(i,6:end); % remove days
end


% make 1989~present
k=0
for i = 1989:2019
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0; m=0;
for i = 1:31
    l=0
        ref_yy(i,:)=[num2str(i+1988)];
    for n = 1:12
        m = m+1;
        ref_yymm(m,:)=[num2str(i+1988) '.' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mth_d_txt(k,:)=[num2str(i+1988) '.' num2str(n,'%02d') '.'  num2str(j,'%02d')];
        mmdd(k,:)=[num2str(n,'%02d') '.'  num2str(j,'%02d')];
    end
    end
end

% make 366 mm-dd
for i = 1:12
  eom_80(i) = eomday(1980,i); % 1980 is leap-yr
end

k=0
for i = 1:12
    l=0
    for j = 1:eom_80(i)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mmdd_366(k,:)=[num2str(i,'%02d') '.'  num2str(l,'%02d')]
    end
end

% make it to cell-array
for i = 1:length(mth_d_txt)
    mth_d_txt_c{i,1} = mth_d_txt(i,:); % raw
end

for i = 1:length(ref_yymm)
        ref_yymm_c{i,1} = ref_yymm(i,:);
end


for i = 1:length(ref_yy)
        ref_yy_c{i,1} = ref_yy(i,:);
end

for i = 1:length(mmdd)
        mmdd_c{i,1} = mmdd(i,:);
end

for i = 1:length(mmdd_366)
        mmdd_366_c{i,1} = mmdd_366(i,:);
end

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, r_date_txt_c)] == 1);
end

for i = 1:length(ref_yymm_c)
           indx_mth{i} = find([strcmp(ref_yymm_c{i}, r_yymm_c)] == 1);
end

for i = 1:length(ref_yy_c)
           indx_yy{i} = find([strcmp(ref_yy_c{i}, r_yy_c)] == 1);
end

for i = 1:length(mmdd_c)
           indx_366{i} = find([strcmp(mmdd_c{i}, mmdd_366_c)] == 1);
end

% how many times obs. on daily matter. 
for i = 1:length(indx); size_in(i)=length(indx{i}) ; end

for i = 1:length(indx_mth); size_mth(i)=length(indx_mth{i}) ; end

for i = 1:length(indx_yy); size_yy(i)=length(indx_yy{i}) ; end

% plot(size_yy,'.')


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

for i = 1:length(ref_yymm_c)
    if size(indx_mth{i},1) == 0
        w_do_mth(i) = NaN;
        w_chl_mth(i) = NaN; 
        w_no3_mth(i) = NaN; 
        w_nh4_mth(i) = NaN; 
    else
        w_do_mth(i) = nanmean(r_do(indx_mth{i}));
        w_chl_mth(i) = nanmean(r_chl(indx_mth{i}));
        w_no3_mth(i) = nanmean(r_no3(indx_mth{i}));
        w_nh4_mth(i) = nanmean(r_nh4(indx_mth{i}));
    end
end

for i = 1:length(ref_yy_c)
    if size(indx_yy{i},1) == 0
        w_do_yy(i) = NaN;
        w_chl_yy(i) = NaN; 
        w_no3_yy(i) = NaN; 
        w_nh4_mth(i) = NaN; 
    else
        w_do_yy(i) = nanmean(r_do(indx_yy{i}));
        w_chl_yy(i) = nanmean(r_chl(indx_yy{i}));
        w_no3_yy(i) = nanmean(r_no3(indx_yy{i}));
        w_nh4_yy(i) = nanmean(r_nh4(indx_yy{i}));
    end
end

return

clearvars b1 yCalc1 X b yCalc2
% raw
 figure;
 plot(w_do,'-bo');
 xlabel('time(year)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title(['sumjin(songjung)-polynomial DO yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13)
%regression
%slope y = b1*x
b1 = [1:length(w_do)]'\w_do'; % x\y for getting slop
yCalc1 = b1*[1:length(w_do)];
% Slope & Intercept y = b0 + b1*x
X = [ones(length(length(w_do)),1) [1:length(w_do)]]; %b_0 b_1 
b = X\w_do;
yCalc2 = X*b;
plot([1:length(w_do)],yCalc2,'--','color','r')
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)

 
figure;
plot(w_chl,'-bo');
xlabel('time(year)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title(['sumjin(songjung)-polynomial Chla yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13)


figure;
plot(w_nh4,'-bo');
xlabel('time(year)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title(['sumjin(songjung)-polynomial NH4 yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13);

 
figure;
plot(w_no3,'-bo');
xlabel('time(year)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title(['sumjin(songjung)-polynomial NO3 yearly mean'],'fontsize',13)
grid on
set(gca,'fontsize',13)




 % month-timeseries
for i = 1:12 
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_do_mth_f = w_do_mth(i:12:end); 
figure; hold on;
plot(w_do_mth_f,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title(['sachoen-regression DO ' num2str(i) '-month'],'fontsize',13)
grid on; 
ylim([0 22])
set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_do_mth_f(~isnan(w_do_mth_f));
idx = find(isnan(w_do_mth_f) == 0);
nanaxisx=find(isnan(w_do_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-regression DO ' num2str(i) '-month']);
end
 
for i = 1:12
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_chl_mth_f = w_chl_mth(i:12:end); 
 figure; hold on;
  plot(w_chl_mth_f,'-bo','linew',2)
  xlabel('time(year)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title(['sachoen-regression Chla ' num2str(i) '-month'],'fontsize',13)
grid on
set(gca,'fontsize',13)
  ylim([0 145])
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_chl_mth_f(~isnan(w_chl_mth_f));
idx = find(isnan(w_chl_mth_f) == 0);
nanaxisx=find(isnan(w_chl_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-regression Chla ' num2str(i) '-month']);
end
 
for i = 1:12
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_nh4_mth_f = w_nh4_mth(i:12:end);
     figure; hold on;
  plot(w_nh4_mth_f,'-bo','linew',2)
xlabel('time(year)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title(['sachoen-regression NH4 ' num2str(i) '-month'],'fontsize',13)
grid on
set(gca,'fontsize',13)
    ylim([0 2.5])
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
  %regression
%slope y = b1*x
nonan=w_nh4_mth_f(~isnan(w_nh4_mth_f));
idx = find(isnan(w_nh4_mth_f) == 0);
nanaxisx=find(isnan(w_nh4_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-regression NH4 ' num2str(i) '-month']);
end
 
for i = 1:12
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_no3_mth_f = w_no3_mth(i:12:end); 
figure; hold on;
  plot(w_no3_mth_f,'-bo','linew',2)
  xlabel('time(year)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title(['sachoen-regression NO3 ' num2str(i) '-month'],'fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
      ylim([0 5])
  %regression
%slope y = b1*x
nonan=w_no3_mth_f(~isnan(w_no3_mth_f));
idx = find(isnan(w_no3_mth_f) == 0);
nanaxisx=find(isnan(w_no3_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-regression NO3 ' num2str(i) '-month']);
end

close all;

% month-timeseries
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_do_mth_f = w_do_mth; 
figure; hold on;
plot(w_do_mth,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title(['sachoen-polynomial DO monthly'],'fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:36:12*31]);
set(gca,'xlim',[1 12*31]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_do_mth_f(~isnan(w_do_mth_f));
idx = find(isnan(w_do_mth_f) == 0);
nanaxisx=find(isnan(w_do_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-polynomial DO monthly']);

 
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_chl_mth_f = w_chl_mth; 
figure; hold on;
plot(w_chl_mth,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title(['sachoen-polynomial Chla monthly'],'fontsize',13)
grid on; set(gca,'fontsize',13);
set(gca,'xtick',[1:36:12*31]);
set(gca,'xlim',[1 12*31]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_chl_mth_f(~isnan(w_chl_mth_f));
idx = find(isnan(w_chl_mth_f) == 0);
nanaxisx=find(isnan(w_chl_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-regression Chla monthly']);



clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_nh4_mth_f = w_nh4_mth;
     figure; hold on;
plot(w_nh4_mth,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title(['sachoen-polynomial NH4 monthly'],'fontsize',13)
grid on; set(gca,'fontsize',13);
set(gca,'xtick',[1:36:12*31]);
set(gca,'xlim',[1 12*31]);
set(gca,'xticklabel',1989:3:2019);
  %regression
%slope y = b1*x
nonan=w_nh4_mth_f(~isnan(w_nh4_mth_f));
idx = find(isnan(w_nh4_mth_f) == 0);
nanaxisx=find(isnan(w_nh4_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-regression NH4 monthly']);

 
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_no3_mth_f = w_no3_mth; 
figure; hold on;
plot(w_no3_mth,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title(['sachoen-regression NO3 monthly'],'fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:36:12*31]);
set(gca,'xlim',[1 12*31]);
set(gca,'xticklabel',1989:3:2019);
 %regression
%slope y = b1*x
nonan=w_no3_mth_f(~isnan(w_no3_mth_f));
idx = find(isnan(w_no3_mth_f) == 0);
nanaxisx=find(isnan(w_no3_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-regression NO3 monthly']);

 % year-timeseries
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_do_mth_f = w_do_yy; 
figure; hold on;
plot(w_do_yy,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title(['sachoen-polynomial DO yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_do_mth_f(~isnan(w_do_mth_f));
idx = find(isnan(w_do_mth_f) == 0);
nanaxisx=find(isnan(w_do_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-polynomial DO yearly mean']);
 
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_chl_mth_f = w_chl_yy; 
figure; hold on;
plot(w_chl_yy,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title(['sachoen-polynomial Chla yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_chl_mth_f(~isnan(w_chl_mth_f));
idx = find(isnan(w_chl_mth_f) == 0);
nanaxisx=find(isnan(w_chl_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-polynomial Chla yearly mean']);

% w_nh4_yy(1:3)=NaN;
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_nh4_mth_f = w_nh4_yy;
     figure; hold on;
plot(w_nh4_yy,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title(['sachoen-polynomial NH4 yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13);
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
  %regression
%slope y = b1*x
nonan=w_nh4_mth_f(~isnan(w_nh4_mth_f));
idx = find(isnan(w_nh4_mth_f) == 0);
nanaxisx=find(isnan(w_nh4_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-polynomial NH4 yearly mean']);

 
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_no3_mth_f = w_no3_yy; 
figure; hold on;
plot(w_no3_yy,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title(['sachoen-polynomial NO3 yearly mean'],'fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_no3_mth_f(~isnan(w_no3_mth_f));
idx = find(isnan(w_no3_mth_f) == 0);
nanaxisx=find(isnan(w_no3_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sachoen-polynomial NO3 yearly mean']);

save('sachoen_regression_biovariable_climate.mat'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SUMJIN (songjung)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  gahwa river %%%%%%%%%%
close all; clear; clc; 
cd E:\장기생태\Dynamic\06_river\환경과학원
[raw txt]=xlsread('수질측정지점_하동.xls','검색결과','');
% [raw txt]=xlsread('수질측정지점_남강댐1.xls','검색결과','');

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
    r_date_txt_c{i,1} = r_date_txt(i,:); % raw
    r_yymm_c{i,1} = r_date_txt(i,1:end-3); % remove days 
    r_yy_c{i,1} = r_date_txt(i,1:4); % remove days 
end


% make 1989~present
k=0
for i = 1989:2019
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0; m=0;
for i = 1:31
    l=0
        ref_yy(i,:)=[num2str(i+1988)];
    for n = 1:12
        m = m+1;
        ref_yymm(m,:)=[num2str(i+1988) '.' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mth_d_txt(k,:)=[num2str(i+1988) '.' num2str(n,'%02d') '.'  num2str(j,'%02d')];
    end
    end
end


% make it to cell-array
for i = 1:length(mth_d_txt)
    mth_d_txt_c{i,1} = mth_d_txt(i,:); % raw
end

for i = 1:length(ref_yymm)
        ref_yymm_c{i,1} = ref_yymm(i,:);
end


for i = 1:length(ref_yy)
        ref_yy_c{i,1} = ref_yy(i,:);
end

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, r_date_txt_c)] == 1);
end

for i = 1:length(ref_yymm_c)
           indx_mth{i} = find([strcmp(ref_yymm_c{i}, r_yymm_c)] == 1);
end

for i = 1:length(ref_yy_c)
           indx_yy{i} = find([strcmp(ref_yy_c{i}, r_yy_c)] == 1);
end

% how many times obs. on daily matter. 
for i = 1:length(indx); size_in(i)=length(indx{i}) ; end

for i = 1:length(indx_mth); size_mth(i)=length(indx_mth{i}) ; end

for i = 1:length(indx_yy); size_yy(i)=length(indx_yy{i}) ; end

% plot(size_yy,'.')


return

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

for i = 1:length(ref_yymm_c)
    if size(indx_mth{i},1) == 0
        w_do_mth(i) = NaN;
        w_chl_mth(i) = NaN; 
        w_no3_mth(i) = NaN; 
        w_nh4_mth(i) = NaN; 
    else
        w_do_mth(i) = nanmean(r_do(indx_mth{i}));
        w_chl_mth(i) = nanmean(r_chl(indx_mth{i}));
        w_no3_mth(i) = nanmean(r_no3(indx_mth{i}));
        w_nh4_mth(i) = nanmean(r_nh4(indx_mth{i}));
    end
end

for i = 1:length(ref_yy_c)
    if size(indx_yy{i},1) == 0
        w_do_yy(i) = NaN;
        w_chl_yy(i) = NaN; 
        w_no3_yy(i) = NaN; 
        w_nh4_mth(i) = NaN; 
    else
        w_do_yy(i) = nanmean(r_do(indx_yy{i}));
        w_chl_yy(i) = nanmean(r_chl(indx_yy{i}));
        w_no3_yy(i) = nanmean(r_no3(indx_yy{i}));
        w_nh4_yy(i) = nanmean(r_nh4(indx_yy{i}));
    end
end

clearvars b1 yCalc1 X b yCalc2
% raw
 figure;
 plot(w_do,'-bo');
 xlabel('time(year)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title(['sumjin(hadong)-polynomial DO yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13)
%regression
%slope y = b1*x
b1 = [1:length(w_do)]'\w_do'; % x\y for getting slop
yCalc1 = b1*[1:length(w_do)];
% Slope & Intercept y = b0 + b1*x
X = [ones(length(length(w_do)),1) [1:length(w_do)]]; %b_0 b_1 
b = X\w_do;
yCalc2 = X*b;
plot([1:length(w_do)],yCalc2,'--','color','r')
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)

 
figure;
plot(w_chl,'-bo');
xlabel('time(year)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title(['sumjin(hadong)-polynomial Chla yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13)


figure;
plot(w_nh4,'-bo');
xlabel('time(year)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title(['sumjin(hadong)-polynomial NH4 yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13);

 
figure;
plot(w_no3,'-bo');
xlabel('time(year)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title(['sumjin(hadong)-polynomial NO3 yearly mean'],'fontsize',13)
grid on
set(gca,'fontsize',13)


 % month-timeseries
for i = 1:12 
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_do_mth_f = w_do_mth(i:12:end); 
figure; hold on;
plot(w_do_mth_f,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title(['sumjin(hadong)-regression DO ' num2str(i) '-month'],'fontsize',13)
grid on
ylim([0 22])
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_do_mth_f(~isnan(w_do_mth_f));
idx = find(isnan(w_do_mth_f) == 0);
nanaxisx=find(isnan(w_do_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
if length(nanaxisx) == 0
    plot(1:length(idx),yCalc2,'--','color','r','linew',2)
end
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-regression DO ' num2str(i) '-month']);
end
 
for i = 1:12
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_chl_mth_f = w_chl_mth(i:12:end); 
 figure; hold on;
  plot(w_chl_mth_f,'-bo','linew',2);
  xlabel('time(year)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title(['sumjin(hadong)-regression Chla ' num2str(i) '-month'],'fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
  ylim([0 145])
%regression
%slope y = b1*x
nonan=w_chl_mth_f(~isnan(w_chl_mth_f));
idx = find(isnan(w_chl_mth_f) == 0);
nanaxisx=find(isnan(w_chl_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
if length(nanaxisx) == 0
    plot(1:length(idx),yCalc2,'--','color','r','linew',2)
elseif length(nanaxisx) > 0
    plot(idx,yCalc2,'--','color','r','linew',2)
end
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-regression Chla ' num2str(i) '-month']);
end
 
for i = 1:12
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_nh4_mth_f = w_nh4_mth(i:12:end);
     figure; hold on;
  plot(w_nh4_mth_f,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title(['sumjin(hadong)-regression NH4 ' num2str(i) '-month'],'fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
  ylim([0 2.5])
  %regression
%slope y = b1*x
nonan=w_nh4_mth_f(~isnan(w_nh4_mth_f));
idx = find(isnan(w_nh4_mth_f) == 0);
nanaxisx=find(isnan(w_nh4_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
if length(nanaxisx) == 0
    plot(1:length(idx),yCalc2,'--','color','r','linew',2)
elseif length(nanaxisx) > 0
    plot(idx,yCalc2,'--','color','r','linew',2)
end
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-regression NH4 ' num2str(i) '-month']);
end
 
for i = 1:12
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_no3_mth_f = w_no3_mth(i:12:end); 
figure; hold on;
  plot(w_no3_mth_f,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title(['sumjin(hadong)-regression NO3 ' num2str(i) '-month'],'fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
  ylim([0 11])
  %regression
%slope y = b1*x
nonan=w_no3_mth_f(~isnan(w_no3_mth_f));
idx = find(isnan(w_no3_mth_f) == 0);
nanaxisx=find(isnan(w_no3_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-regression NO3 ' num2str(i) '-month']);
end

close all;

% month-timeseries
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_do_mth_f = w_do_mth; 
figure; hold on;
plot(w_do_mth,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title(['sumjin(hadong)-polynomial DO monthly'],'fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:36:12*31]);
set(gca,'xlim',[1 12*31]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_do_mth_f(~isnan(w_do_mth_f));
idx = find(isnan(w_do_mth_f) == 0);
nanaxisx=find(isnan(w_do_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
if length(nanaxisx) == 0
    plot(1:length(idx),yCalc2,'--','color','r','linew',2)
elseif length(nanaxisx) > 0
    plot(idx,yCalc2,'--','color','r','linew',2)
end
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-polynomial DO monthly']);

 
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_chl_mth_f = w_chl_mth; 
figure; hold on;
plot(w_chl_mth,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title(['sumjin(hadong)-polynomial Chla monthly'],'fontsize',13)
grid on; set(gca,'fontsize',13);
set(gca,'xtick',[1:36:12*31]);
set(gca,'xlim',[1 12*31]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_chl_mth_f(~isnan(w_chl_mth_f));
idx = find(isnan(w_chl_mth_f) == 0);
nanaxisx=find(isnan(w_chl_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-regression Chla monthly']);



clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_nh4_mth_f = w_nh4_mth;
     figure; hold on;
plot(w_nh4_mth,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title(['sumjin(hadong)-polynomial NH4 monthly'],'fontsize',13)
grid on; set(gca,'fontsize',13);
set(gca,'xtick',[1:36:12*31]);
set(gca,'xlim',[1 12*31]);
set(gca,'xticklabel',1989:3:2019);
  %regression
%slope y = b1*x
nonan=w_nh4_mth_f(~isnan(w_nh4_mth_f));
idx = find(isnan(w_nh4_mth_f) == 0);
nanaxisx=find(isnan(w_nh4_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-regression NH4 monthly']);

 
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_no3_mth_f = w_no3_mth; 
figure; hold on;
plot(w_no3_mth,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title(['sumjin(hadong)-regression NO3 monthly'],'fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:36:12*31]);
set(gca,'xlim',[1 12*31]);
set(gca,'xticklabel',1989:3:2019);
 %regression
%slope y = b1*x
nonan=w_no3_mth_f(~isnan(w_no3_mth_f));
idx = find(isnan(w_no3_mth_f) == 0);
nanaxisx=find(isnan(w_no3_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-regression NO3 monthly']);

 % year-timeseries
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_do_mth_f = w_do_yy; 
figure; hold on;
plot(w_do_yy,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title(['sumjin(hadong)-polynomial DO yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_do_mth_f(~isnan(w_do_mth_f));
idx = find(isnan(w_do_mth_f) == 0);
nanaxisx=find(isnan(w_do_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
if length(nanaxisx) == 0
    plot(1:length(idx),yCalc2,'--','color','r','linew',2)
elseif length(nanaxisx) > 0
    plot(idx,yCalc2,'--','color','r','linew',2)
end
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-polynomial DO yearly mean']);


clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_chl_mth_f = w_chl_yy; 
figure; hold on;
plot(w_chl_yy,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title(['sumjin(hadong)-polynomial Chla yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_chl_mth_f(~isnan(w_chl_mth_f));
idx = find(isnan(w_chl_mth_f) == 0);
nanaxisx=find(isnan(w_chl_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-polynomial Chla yearly mean']);


clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_nh4_mth_f = w_nh4_yy;
     figure; hold on;
plot(w_nh4_yy,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title(['sumjin(hadong)-polynomial NH4 yearly mean'],'fontsize',13)
grid on; set(gca,'fontsize',13);
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
  %regression
%slope y = b1*x
nonan=w_nh4_mth_f(~isnan(w_nh4_mth_f));
idx = find(isnan(w_nh4_mth_f) == 0);
nanaxisx=find(isnan(w_nh4_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-polynomial NH4 yearly mean']);

 
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx
w_no3_mth_f = w_no3_yy; 
figure; hold on;
plot(w_no3_yy,'-bo','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title(['sumjin(hadong)-polynomial NO3 yearly mean'],'fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:3:32]);
set(gca,'xlim',[1 32]);
set(gca,'xticklabel',1989:3:2019);
%regression
%slope y = b1*x
nonan=w_no3_mth_f(~isnan(w_no3_mth_f));
idx = find(isnan(w_no3_mth_f) == 0);
nanaxisx=find(isnan(w_no3_mth_f) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = X*b;
plot(idx,yCalc2,'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
print('-dpng', ['sumjin(hadong)-polynomial NO3 yearly mean']);

save('sumjin(hadong)_regression_biovariable_climate.mat');     