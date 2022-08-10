%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  gahwa river %%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\환경과학원
[raw txt]=xlsread('수질측정지점_구례_fix.xls','검색결과','');
% [raw txt]=xlsread('수질측정지점_사천천.xls','검색결과','');
load('songjung_polynomial_climate.mat','yp_*');

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
    r_mm_c{i,1} = r_date_txt(i,6:7); % remove days
end


% make 1989~present
k=0
for i = 1980:2019
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
for i = 1:40
    if i ==1
         t_tick(i) = t_tick_pre(i);
    else
        t_tick(i) = sum(t_tick_pre(1:i));
    end
end

k=0; m=0;
for i = 1:40
    l=0
        ref_yy(i,:)=[num2str(i+1979)];
    for n = 1:12
        m = m+1;
        ref_yymm(m,:)=[num2str(i+1979) '.' num2str(n,'%02d')];        
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mth_d_txt(k,:)=[num2str(i+1979) '.' num2str(n,'%02d') '.'  num2str(j,'%02d')];
        mmdd(k,:)=[num2str(n,'%02d') '.'  num2str(j,'%02d')];
        order_mm(k,:)=[num2str(n,'%02d')];  % 1989~2019 only month
    end
    end
end

for i = 1:12
ref_mm(i,:)=[num2str(i,'%02d')]; % 1~12mth
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

for i = 1:length(ref_mm)
        ref_mm_c{i,1} = ref_mm(i,:);
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

for i = 1:length(order_mm)
        order_mm_c{i,1} = order_mm(i,:); % 1989~2019 only month
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

% 1989~2019 already oredered matrix indx
for i = 1:length(ref_mm_c)
           indx_12mth{i} = find([strcmp(ref_mm_c{i}, order_mm_c)] == 1);
end


% how many times obs. on daily matter. 
for i = 1:length(indx); size_in(i)=length(indx{i}) ; end

for i = 1:length(indx_mth); size_mth(i)=length(indx_mth{i}) ; end

for i = 1:length(indx_yy); size_yy(i)=length(indx_yy{i}) ; end

% plot(size_yy,'.')


%% make quasi_climate
% raw
% make raw - climate (daily matters)
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do_mi(i) = NaN;
        w_chl_mi(i) = NaN; 
        w_no3_mi(i) = NaN; 
        w_nh4_mi(i) = NaN; 
    else
        w_do_mi(i) = nanmean(r_do(indx{i})) - yp_w_do(indx_366{i});
        w_chl_mi(i) = nanmean(r_chl(indx{i})) - yp_w_chl(indx_366{i});
        w_no3_mi(i) = nanmean(r_no3(indx{i})) - yp_w_no3(indx_366{i});
        w_nh4_mi(i) = nanmean(r_nh4(indx{i})) - yp_w_nh4(indx_366{i});        
    end
end


% make raw - climate (monthly matters)
for i = 1:length(mth_d_txt_c)
        do_day_clim(i) = yp_w_do(indx_366{i});
        chl_day_clim(i) = yp_w_chl(indx_366{i});
        no3_day_clim(i) = yp_w_no3(indx_366{i});
        nh4_day_clim(i) = yp_w_nh4(indx_366{i});
end

for i = 1:length(t_indx)
    if i == 1
    do_mth_clim(i) = mean(do_day_clim(1:t_indx(i)));
    chl_mth_clim(i) = mean(chl_day_clim(1:t_indx(i)));
    no3_mth_clim(i) =  mean(no3_day_clim(1:t_indx(i)));
    nh4_mth_clim(i) = mean(nh4_day_clim(1:t_indx(i)));
    else
    do_mth_clim(i) = mean(do_day_clim(t_indx(i-1)+1:t_indx(i)));
    chl_mth_clim(i) = mean(chl_day_clim(t_indx(i-1)+1:t_indx(i)));
    no3_mth_clim(i) =  mean(no3_day_clim(t_indx(i-1)+1:t_indx(i)));
    nh4_mth_clim(i) = mean(nh4_day_clim(t_indx(i-1)+1:t_indx(i)));
    end
end

%  
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

tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_do,'.');
 xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_do_04=w_do(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_do_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_do_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_do_04)],mean(w_do_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

    for i = 1:length(yCalc3)
            w_do_recon(i) = yCalc3(i)  + yp_w_do(indx_366{i});
    end
plot([1:length(w_do)], w_do_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_do_04=w_do(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)],mean(w_do_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_do_recon(i) = yCalc3(i)  + yp_w_do(indx_366{i});
%     end
% plot([1:length(w_do)], w_do_recon,'-.','color',color_spec(sig),'linew',1)
end

title(['sumjin(songjung) do raw  ' num2str(sig) '-sigma'],'fontsize',13)


tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_chl,'.');
 xlabel('time(year)','fontsize',13)
ylabel('chl (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_chl_04=w_chl(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_chl_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_chl_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_chl_04)],mean(w_chl_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_chl_recon(i) = yCalc3(i)  + yp_w_chl(indx_366{i});
%     end
% plot([1:length(w_chl)], w_chl_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_chl_04=w_chl(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)],mean(w_chl_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_chl_recon(i) = yCalc3(i)  + yp_w_chl(indx_366{i});
%     end
% plot([1:length(w_chl)], w_chl_recon,'-.','color',color_spec(sig),'linew',1)
end
ylim([0 50])
title(['sumjin(songjung) chl raw  ' num2str(sig) '-sigma'],'fontsize',13)



tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_nh4,'.');
 xlabel('time(year)','fontsize',13)
ylabel('nh4 (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_nh4_04=w_nh4(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_nh4_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_nh4_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_nh4_04)],mean(w_nh4_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_nh4_recon(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
%     end
% plot([1:length(w_nh4)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_nh4_04=w_nh4(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)],mean(w_nh4_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_nh4_recon(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
%     end
% plot([1:length(w_nh4)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
end
ylim([0 inf])
title(['sumjin(songjung) nh4 raw  ' num2str(sig) '-sigma'],'fontsize',13)


tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_no3,'.');
 xlabel('time(year)','fontsize',13)
ylabel('no3 (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_no3_04=w_no3(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_no3_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_no3_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_no3_04)],mean(w_no3_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_no3_recon(i) = yCalc3(i)  + yp_w_no3(indx_366{i});
%     end
% plot([1:length(w_no3)], w_no3_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_no3_04=w_no3(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)],mean(w_no3_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_no3_recon(i) = yCalc3(i)  + yp_w_no3(indx_366{i});
%     end
% plot([1:length(w_no3)], w_no3_recon,'-.','color',color_spec(sig),'linew',1)
end
ylim([0 inf])
title(['sumjin(songjung) no3 raw  ' num2str(sig) '-sigma'],'fontsize',13)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% %%% %%% %%% %%% %%%  %%% Namgang river %%% %%% %%% %%% %%% %%% %%% %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\환경과학원
% [raw txt]=xlsread('수질측정지점_하동.xls','검색결과','');
[raw txt]=xlsread('수질측정지점_남강댐1.xls','검색결과','');
load('namgang_polynomial_climate.mat','yp_*');

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
    r_mm_c{i,1} = r_date_txt(i,6:7); % remove days
end


% make 1989~present
k=0
for i = 1980:2019
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
for i = 1:40
    if i ==1
         t_tick(i) = t_tick_pre(i);
    else
        t_tick(i) = sum(t_tick_pre(1:i));
    end
end

k=0; m=0;
for i = 1:40
    l=0
        ref_yy(i,:)=[num2str(i+1979)];
    for n = 1:12
        m = m+1;
        ref_yymm(m,:)=[num2str(i+1979) '.' num2str(n,'%02d')];        
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mth_d_txt(k,:)=[num2str(i+1979) '.' num2str(n,'%02d') '.'  num2str(j,'%02d')];
        mmdd(k,:)=[num2str(n,'%02d') '.'  num2str(j,'%02d')];
        order_mm(k,:)=[num2str(n,'%02d')];  % 1989~2019 only month
    end
    end
end

for i = 1:12
ref_mm(i,:)=[num2str(i,'%02d')]; % 1~12mth
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

for i = 1:length(ref_mm)
        ref_mm_c{i,1} = ref_mm(i,:);
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

for i = 1:length(order_mm)
        order_mm_c{i,1} = order_mm(i,:); % 1989~2019 only month
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

% 1989~2019 already oredered matrix indx
for i = 1:length(ref_mm_c)
           indx_12mth{i} = find([strcmp(ref_mm_c{i}, order_mm_c)] == 1);
end


% how many times obs. on daily matter. 
for i = 1:length(indx); size_in(i)=length(indx{i}) ; end

for i = 1:length(indx_mth); size_mth(i)=length(indx_mth{i}) ; end

for i = 1:length(indx_yy); size_yy(i)=length(indx_yy{i}) ; end

% plot(size_yy,'.')



%% make quasi_climate

% make raw - climate (daily matters)
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do_mi(i) = NaN;
        w_chl_mi(i) = NaN; 
        w_no3_mi(i) = NaN; 
        w_nh4_mi(i) = NaN; 
    else
        w_do_mi(i) = nanmean(r_do(indx{i})) - yp_w_do(indx_366{i});
        w_chl_mi(i) = nanmean(r_chl(indx{i})) - yp_w_chl(indx_366{i});
        w_no3_mi(i) = nanmean(r_no3(indx{i})) - yp_w_no3(indx_366{i});
        w_nh4_mi(i) = nanmean(r_nh4(indx{i})) - yp_w_nh4(indx_366{i});
        
    end
end


% make raw - climate (monthly matters)
for i = 1:length(mth_d_txt_c)
        do_day_clim(i) = yp_w_do(indx_366{i});
        chl_day_clim(i) = yp_w_chl(indx_366{i});
        no3_day_clim(i) = yp_w_no3(indx_366{i});
        nh4_day_clim(i) = yp_w_nh4(indx_366{i});
end

for i = 1:length(t_indx)
    if i == 1
    do_mth_clim(i) = mean(do_day_clim(1:t_indx(i)));
    chl_mth_clim(i) = mean(chl_day_clim(1:t_indx(i)));
    no3_mth_clim(i) =  mean(no3_day_clim(1:t_indx(i)));
    nh4_mth_clim(i) = mean(nh4_day_clim(1:t_indx(i)));
    else
    do_mth_clim(i) = mean(do_day_clim(t_indx(i-1)+1:t_indx(i)));
    chl_mth_clim(i) = mean(chl_day_clim(t_indx(i-1)+1:t_indx(i)));
    no3_mth_clim(i) =  mean(no3_day_clim(t_indx(i-1)+1:t_indx(i)));
    nh4_mth_clim(i) = mean(nh4_day_clim(t_indx(i-1)+1:t_indx(i)));
    end
end

%  

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

tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_do,'.');
 xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_do_04=w_do(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_do_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_do_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_do_04)],mean(w_do_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_do_recon(i) = yCalc3(i)  + yp_w_do(indx_366{i});
%     end
% plot([1:length(w_do)], w_do_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_do_04=w_do(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)],mean(w_do_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_do_recon(i) = yCalc3(i)  + yp_w_do(indx_366{i});
%     end
% plot([1:length(w_do)], w_do_recon,'-.','color',color_spec(sig),'linew',1)
end
title(['Namgang do raw  ' num2str(sig) '-sigma'],'fontsize',13)


tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_chl,'.');
 xlabel('time(year)','fontsize',13)
ylabel('chl (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_chl_04=w_chl(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_chl_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_chl_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_chl_04)],mean(w_chl_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_chl_recon(i) = yCalc3(i)  + yp_w_chl(indx_366{i});
%     end
% plot([1:length(w_chl)], w_chl_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_chl_04=w_chl(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)],mean(w_chl_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_chl_recon(i) = yCalc3(i)  + yp_w_chl(indx_366{i});
%     end
% plot([1:length(w_chl)], w_chl_recon,'-.','color',color_spec(sig),'linew',1)
end
ylim([0 50])
title(['Namgang chl raw  ' num2str(sig) '-sigma'],'fontsize',13)



tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_nh4,'.');
 xlabel('time(year)','fontsize',13)
ylabel('nh4 (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_nh4_04=w_nh4(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_nh4_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_nh4_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_nh4_04)],mean(w_nh4_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_nh4_recon(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
%     end
% plot([1:length(w_nh4)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_nh4_04=w_nh4(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)],mean(w_nh4_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_nh4_recon(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
%     end
% plot([1:length(w_nh4)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
end
ylim([0 inf])
title(['Namgang nh4 raw  ' num2str(sig) '-sigma'],'fontsize',13)


tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_no3,'.');
 xlabel('time(year)','fontsize',13)
ylabel('no3 (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_no3_04=w_no3(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_no3_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_no3_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_no3_04)],mean(w_no3_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_no3_recon(i) = yCalc3(i)  + yp_w_no3(indx_366{i});
%     end
% plot([1:length(w_no3)], w_no3_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_no3_04=w_no3(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)],mean(w_no3_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_no3_recon(i) = yCalc3(i)  + yp_w_no3(indx_366{i});
%     end
% plot([1:length(w_no3)], w_no3_recon,'-.','color',color_spec(sig),'linew',1)
end
ylim([0 inf])
title(['Namgang no3 raw  ' num2str(sig) '-sigma'],'fontsize',13)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% jinwal river %%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\환경과학원
[raw txt]=xlsread('수질측정지점_진월_수정.xls','검색결과','');
% [raw txt]=xlsread('수질측정지점_사천천.xls','검색결과','');
load('songjung_polynomial_climate.mat','yp_*');

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
    r_mm_c{i,1} = r_date_txt(i,6:7); % remove days
end


% make 1989~present
k=0
for i = 1980:2019
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
for i = 1:40
    if i ==1
         t_tick(i) = t_tick_pre(i);
    else
        t_tick(i) = sum(t_tick_pre(1:i));
    end
end

k=0; m=0;
for i = 1:40
    l=0
        ref_yy(i,:)=[num2str(i+1979)];
    for n = 1:12
        m = m+1;
        ref_yymm(m,:)=[num2str(i+1979) '.' num2str(n,'%02d')];        
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mth_d_txt(k,:)=[num2str(i+1979) '.' num2str(n,'%02d') '.'  num2str(j,'%02d')];
        mmdd(k,:)=[num2str(n,'%02d') '.'  num2str(j,'%02d')];
        order_mm(k,:)=[num2str(n,'%02d')];  % 1989~2019 only month
    end
    end
end

for i = 1:12
ref_mm(i,:)=[num2str(i,'%02d')]; % 1~12mth
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

for i = 1:length(ref_mm)
        ref_mm_c{i,1} = ref_mm(i,:);
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

for i = 1:length(order_mm)
        order_mm_c{i,1} = order_mm(i,:); % 1989~2019 only month
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

% 1989~2019 already oredered matrix indx
for i = 1:length(ref_mm_c)
           indx_12mth{i} = find([strcmp(ref_mm_c{i}, order_mm_c)] == 1);
end


% how many times obs. on daily matter. 
for i = 1:length(indx); size_in(i)=length(indx{i}) ; end

for i = 1:length(indx_mth); size_mth(i)=length(indx_mth{i}) ; end

for i = 1:length(indx_yy); size_yy(i)=length(indx_yy{i}) ; end

% plot(size_yy,'.')


%% make quasi_climate
% raw
% make raw - climate (daily matters)
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_do_mi(i) = NaN;
        w_chl_mi(i) = NaN; 
        w_no3_mi(i) = NaN; 
        w_nh4_mi(i) = NaN; 
    else
        w_do_mi(i) = nanmean(r_do(indx{i})) - yp_w_do(indx_366{i});
        w_chl_mi(i) = nanmean(r_chl(indx{i})) - yp_w_chl(indx_366{i});
        w_no3_mi(i) = nanmean(r_no3(indx{i})) - yp_w_no3(indx_366{i});
        w_nh4_mi(i) = nanmean(r_nh4(indx{i})) - yp_w_nh4(indx_366{i});        
    end
end


% make raw - climate (monthly matters)
for i = 1:length(mth_d_txt_c)
        do_day_clim(i) = yp_w_do(indx_366{i});
        chl_day_clim(i) = yp_w_chl(indx_366{i});
        no3_day_clim(i) = yp_w_no3(indx_366{i});
        nh4_day_clim(i) = yp_w_nh4(indx_366{i});
end

for i = 1:length(t_indx)
    if i == 1
    do_mth_clim(i) = mean(do_day_clim(1:t_indx(i)));
    chl_mth_clim(i) = mean(chl_day_clim(1:t_indx(i)));
    no3_mth_clim(i) =  mean(no3_day_clim(1:t_indx(i)));
    nh4_mth_clim(i) = mean(nh4_day_clim(1:t_indx(i)));
    else
    do_mth_clim(i) = mean(do_day_clim(t_indx(i-1)+1:t_indx(i)));
    chl_mth_clim(i) = mean(chl_day_clim(t_indx(i-1)+1:t_indx(i)));
    no3_mth_clim(i) =  mean(no3_day_clim(t_indx(i-1)+1:t_indx(i)));
    nh4_mth_clim(i) = mean(nh4_day_clim(t_indx(i-1)+1:t_indx(i)));
    end
end

%  
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

tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_do,'.');
 xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_do_04=w_do(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_do_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_do_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_do_04)],mean(w_do_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_do_recon(i) = yCalc3(i)  + yp_w_do(indx_366{i});
%     end
% plot([1:length(w_do)], w_do_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_do_04=w_do(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_do_04)],mean(w_do_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_do_recon(i) = yCalc3(i)  + yp_w_do(indx_366{i});
%     end
% plot([1:length(w_do)], w_do_recon,'-.','color',color_spec(sig),'linew',1)
end

title(['sumjin(songjung) do raw  ' num2str(sig) '-sigma'],'fontsize',13)


tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_chl,'.');
 xlabel('time(year)','fontsize',13)
ylabel('chl (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_chl_04=w_chl(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_chl_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_chl_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_chl_04)],mean(w_chl_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_chl_recon(i) = yCalc3(i)  + yp_w_chl(indx_366{i});
%     end
% plot([1:length(w_chl)], w_chl_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_chl_04=w_chl(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_chl_04)],mean(w_chl_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_chl_recon(i) = yCalc3(i)  + yp_w_chl(indx_366{i});
%     end
% plot([1:length(w_chl)], w_chl_recon,'-.','color',color_spec(sig),'linew',1)
end
ylim([0 50])
title(['sumjin(songjung) chl raw  ' num2str(sig) '-sigma'],'fontsize',13)



tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_nh4,'.');
 xlabel('time(year)','fontsize',13)
ylabel('nh4 (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_nh4_04=w_nh4(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_nh4_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_nh4_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_nh4_04)],mean(w_nh4_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_nh4_recon(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
%     end
% plot([1:length(w_nh4)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_nh4_04=w_nh4(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_nh4_04)],mean(w_nh4_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_nh4_recon(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
%     end
% plot([1:length(w_nh4)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
end
ylim([0 inf])
title(['sumjin(songjung) nh4 raw  ' num2str(sig) '-sigma'],'fontsize',13)


tx_tick = [1 t_tick(1:end)];
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_no3,'.');
 xlabel('time(year)','fontsize',13)
ylabel('no3 (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:4:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_no3_04=w_no3(1:tx_tick(25));
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
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_no3_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_no3_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_no3_04)],mean(w_no3_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_no3_recon(i) = yCalc3(i)  + yp_w_no3(indx_366{i});
%     end
% plot([1:length(w_no3)], w_no3_recon,'-.','color',color_spec(sig),'linew',1)
end

for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_no3_04=w_no3(tx_tick(25)+1:end);
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
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)],yCalc2{sig},'-.','color','k','linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([tx_tick(25)+1:1:tx_tick(25)+length(w_no3_04)],mean(w_no3_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_no3_recon(i) = yCalc3(i)  + yp_w_no3(indx_366{i});
%     end
% plot([1:length(w_no3)], w_no3_recon,'-.','color',color_spec(sig),'linew',1)
end
ylim([0 inf])
title(['sumjin(songjung) no3 raw  ' num2str(sig) '-sigma'],'fontsize',13)








