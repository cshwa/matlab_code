%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SUMJIN (songjung)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  gahwa river %%%%%%%%%%
close all; clear; clc; 
cd E:\장기생태\Dynamic\06_river\환경과학원
% [raw txt]=xlsread('수질측정지점_하동.xls','검색결과','');
[raw txt]=xlsread('수질측정지점_남강댐1.xls','검색결과','');
load('Namgang_polynomial_climate.mat','yp_*');
load('ycf_nh4.mat')

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

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_do,'.');
 xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
for sig=1:6  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
idx1 = find(isnan(w_do_mi) == 0);
w_do_mi2=w_do_mi;
upper_bc(sig) = mean(w_do_mi(idx1)) + sig*std(w_do_mi(idx1));
lower_bc(sig) = mean(w_do_mi(idx1)) - sig*std(w_do_mi(idx1));
mean(w_do_mi(idx1))
w_do_mi2(find(w_do_mi > mean(w_do_mi(idx1)) + sig*std(w_do_mi(idx1))))=NaN;
w_do_mi2(find(w_do_mi < mean(w_do_mi(idx1)) - sig*std(w_do_mi(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_do_mi2(~isnan(w_do_mi2));
idx = find(isnan(w_do_mi2) == 0);
nanaxisx=find(isnan(w_do_mi2) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_do_mi)].*b1(sig) + b0(sig);
% plot(yCalc2{sig},'-.','color',color_spec(sig),'linew',1)
% line([1:1:length(w_do_mi)], upper_bc(sig),'color', color_spec(sig),'linew',2)
% line([1:1:length(w_do_mi)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
    for i = 1:length(yCalc3)
            w_do_recon(i) = yCalc3(i)  + yp_w_do(indx_366{i});
    end
plot([1:length(w_do_mi)], w_do_recon,'-.','color',color_spec(sig),'linew',1)
end
title(['Namgang-2-polynomial do raw - climate recon ' num2str(sig) '-sigma'],'fontsize',13)
print('-dpng', ['Namgang-2_do_raw-climate_recon ' num2str(sig) '-sigma']);
set(gca,'xtick',[1 t_tick(1:1:end)]);
set(gca,'xlim',[1 t_tick(10)]);
set(gca,'xticklabel',1980:1:2020);    
print('-dpng', ['Namgang-2_do_raw-climate-climate_recon_80s ' num2str(sig) '-sigma']);

%%%%%%%%%%%%%%%%%%%%%
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_chl,'.');
 xlabel('time(year)','fontsize',13)
ylabel('chl (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
for sig=1:6  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
idx1 = find(isnan(w_chl_mi) == 0);
w_chl_mi2=w_chl_mi;
upper_bc(sig) = mean(w_chl_mi(idx1)) + sig*std(w_chl_mi(idx1));
lower_bc(sig) = mean(w_chl_mi(idx1)) - sig*std(w_chl_mi(idx1));
mean(w_chl_mi(idx1))
w_chl_mi2(find(w_chl_mi > mean(w_chl_mi(idx1)) + sig*std(w_chl_mi(idx1))))=NaN;
w_chl_mi2(find(w_chl_mi < mean(w_chl_mi(idx1)) - sig*std(w_chl_mi(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_chl_mi2(~isnan(w_chl_mi2));
idx = find(isnan(w_chl_mi2) == 0);
nanaxisx=find(isnan(w_chl_mi2) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_chl_mi)].*b1(sig) + b0(sig);
% plot(yCalc2{sig},'-.','color',color_spec(sig),'linew',1)
% line([1:1:length(w_chl_mi)], upper_bc(sig),'color', color_spec(sig),'linew',2)
% line([1:1:length(w_chl_mi)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
    for i = 1:length(yCalc3)
            w_chl_recon(i) = yCalc3(i)  + yp_w_chl(indx_366{i});
    end
plot([1:length(w_chl_mi)], w_chl_recon,'-.','color',color_spec(sig),'linew',1)
end
ylim([-inf 30])
title(['Namgang-2-polynomial chl raw - climate recon ' num2str(sig) '-sigma'],'fontsize',13)
print('-dpng', ['Namgang-2_chl_raw-climate_recon ' num2str(sig) '-sigma']);
set(gca,'xtick',[1 t_tick(1:1:end)]);
set(gca,'xlim',[1 t_tick(10)]);
set(gca,'xticklabel',1980:1:2020);    
print('-dpng', ['Namgang-2_chl_raw-climateclimate_recon_80s ' num2str(sig) '-sigma']);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc yCalc3 yCalc4 *_recon_1 *_recon_2
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_nh4,'.');
 xlabel('time(year)','fontsize',13)
ylabel('nh4 (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
for i = 1:length(ycf_nh4_3sig)
            w_nh4_recon_cf(i) = ycf_nh4_3sig(i)  + yp_w_nh4(indx_366{i});
end
plot([1:length(w_nh4_mi)], w_nh4_recon_cf,'-.','color',color_spec(sig),'linew',1)
title(['Namgang-2-polynomial nh4 raw - climate recon ' num2str(sig) '-sigma'],'fontsize',13)
print('-dpng', ['Namgang-2_nh4_raw-climate_recon ' num2str(sig) '-sigma']);
set(gca,'xtick',[1 t_tick(1:1:end)]);
set(gca,'xlim',[1 t_tick(10)]);
set(gca,'xticklabel',1980:1:2020);    
print('-dpng', ['Namgang-2-regression_nh4_raw-climate_month(daily)_recon_80s' num2str(sig) '-sigma']);
    
 

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc yCalc3 yCalc4 *_recon_1 *_recon_2
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_nh4,'.');
 xlabel('time(year)','fontsize',13)
ylabel('nh4 (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
for sig=1:6  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3 *_dx
% raw
idx1 = find(isnan(w_nh4_mi) == 0);
w_nh4_mi2=w_nh4_mi;
upper_bc(sig) = mean(w_nh4_mi(idx1)) + sig*std(w_nh4_mi(idx1));
lower_bc(sig) = mean(w_nh4_mi(idx1)) - sig*std(w_nh4_mi(idx1));
mean(w_nh4_mi(idx1))
w_nh4_mi2(find(w_nh4_mi > mean(w_nh4_mi(idx1)) + sig*std(w_nh4_mi(idx1))))=NaN;
w_nh4_mi2(find(w_nh4_mi < mean(w_nh4_mi(idx1)) - sig*std(w_nh4_mi(idx1))))=NaN;
 
%regression
%slope y = b1*x
w_nh4_sig{sig}=w_nh4_mi2;
nonan=w_nh4_mi2(~isnan(w_nh4_mi2));
idx = find(isnan(w_nh4_mi2) == 0);
p1_dx = find(idx <= 9559);
p2_dx = find(idx >= 9560);
nanaxisx=find(isnan(w_nh4_mi2) == 1);
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx(p1_dx)'),1) idx(p1_dx)']; %b_0 b_1 
b = X\nonan(p1_dx)';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_nh4_mi)].*b1(sig) + b0(sig);
X = [ones(length(idx(p2_dx)'),1) idx(p2_dx)']; %b_0 b_1 
b = X\nonan(p2_dx)';
b0_1(sig) = b(1);  b1_1(sig) = b(2);
yCalc1{sig} = [1:length(w_nh4_mi)].*b1_1(sig) + b0_1(sig);
% plot(yCalc2{sig},'-.','color',color_spec(sig),'linew',1)
% line([1:1:length(w_nh4_mi)], upper_bc(sig),'color', color_spec(sig),'linew',2)
% line([1:1:length(w_nh4_mi)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
yCalc4 = yCalc1{sig};
    for i = 1:length(yCalc3)
            w_nh4_recon_2(i) = yCalc3(i)  + yp_w_nh4(indx_366{i});
            w_nh4_recon_1(i) = yCalc4(i)  + yp_w_nh4(indx_366{i});
    end
% plot([1:length(w_nh4_mi)], w_nh4_recon,'-.','color',color_spec(sig),'linew',1)
plot([9560:length(w_nh4_mi2)], w_nh4_recon_1(9560:end),'-.','color',color_spec(sig),'linew',1)
plot([1:9559], w_nh4_recon_2(1:9559),'-.','color',color_spec(sig),'linew',1)
end
title(['Namgang-2-polynomial nh4 raw - climate recon ' num2str(sig) '-sigma'],'fontsize',13)
print('-dpng', ['Namgang-2_nh4_raw-climate_recon ' num2str(sig) '-sigma']);
set(gca,'xtick',[1 t_tick(1:1:end)]);
set(gca,'xlim',[1 t_tick(10)]);
set(gca,'xticklabel',1980:1:2020);    
print('-dpng', ['Namgang-2-regression_nh4_raw-climate_month(daily)_recon_80s' num2str(sig) '-sigma']);
    
 
%%%

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
 plot(w_no3,'.');
 xlabel('time(year)','fontsize',13)
ylabel('no3 (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
for sig=1:6  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
idx1 = find(isnan(w_no3_mi) == 0);
w_no3_mi2=w_no3_mi;
upper_bc(sig) = mean(w_no3_mi(idx1)) + sig*std(w_no3_mi(idx1));
lower_bc(sig) = mean(w_no3_mi(idx1)) - sig*std(w_no3_mi(idx1));
mean(w_no3_mi(idx1))
w_no3_mi2(find(w_no3_mi > mean(w_no3_mi(idx1)) + sig*std(w_no3_mi(idx1))))=NaN;
w_no3_mi2(find(w_no3_mi < mean(w_no3_mi(idx1)) - sig*std(w_no3_mi(idx1))))=NaN;
 
%regression
%slope y = b1*x
nonan=w_no3_mi2(~isnan(w_no3_mi2));
idx = find(isnan(w_no3_mi2) == 0);
nanaxisx=find(isnan(w_no3_mi2) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_no3_mi)].*b1(sig) + b0(sig);
% plot(yCalc2{sig},'-.','color',color_spec(sig),'linew',1)
% line([1:1:length(w_no3_mi)], upper_bc(sig),'color', color_spec(sig),'linew',2)
% line([1:1:length(w_no3_mi)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
    for i = 1:length(yCalc3)
            w_no3_recon(i) = yCalc3(i)  + yp_w_no3(indx_366{i});
    end
plot([1:length(w_no3_mi)], w_no3_recon,'-.','color',color_spec(sig),'linew',1)
end
title(['Namgang-2-polynomial no3 raw - climate recon ' num2str(sig) '-sigma'],'fontsize',13)
print('-dpng', ['Namgang-2_no3_raw-climate_recon ' num2str(sig) '-sigma']);
  
%%% 

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
% raw
 figure; hold on;
 plot(w_no3,'.');
 xlabel('time(year)','fontsize',13)
ylabel('no3 (mg/L)','fontsize',13)
title(['Namgang-2-polynomial no3 raw - climate'],'fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
%regression
%slope y = b1*x
nonan=w_no3_mi(~isnan(w_no3_mi));
idx = find(isnan(w_no3_mi) == 0);
nanaxisx=find(isnan(w_no3_mi) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = [1:length(w_no3_mi)].*b(2) + b(1);
for i = 1:length(yCalc2)
        w_no3_recon(i) = yCalc2(i)  + yp_w_no3(indx_366{i});
        yp_w_no3_raw(i) = yp_w_no3(indx_366{i});
end
% plot([1:length(w_chl_mi)], yCalc2,'--','color','c','linew',1)
plot([1:length(w_no3_mi)], w_no3_recon,'--','color','c','linew',1)
% plot([1:length(w_chl_mi)],yp_w_chl_raw,'--','color','g','linew',1)
% print('-dpng', ['Namgang-2_no3_raw&recon']);
%regression
%slope y = b1*x
nonan=w_no3_mi2(~isnan(w_no3_mi2));
idx = find(isnan(w_no3_mi2) == 0);
nanaxisx=find(isnan(w_no3_mi2) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = [1:length(w_no3_mi2)].*b(2) + b(1);
for i = 1:length(yCalc2)
        w_no3_recon(i) = yCalc2(i)  + yp_w_no3(indx_366{i});
        yp_w_no3_raw(i) = yp_w_no3(indx_366{i});
end
% plot([1:length(w_no3_mi2)], yCalc2,'--','color','c','linew',1)
plot([1:length(w_no3_mi2)], w_no3_recon,'--','color','r','linew',1)


%% month-timeseries raw

do_recon_mth = NaN(6,length(w_do));
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 w_do_mi_p sig_txt *_mi_p *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
plot(w_do,'.','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
grid on; 
ylim([0 inf])
% set(gca,'ytick',[nanmin(w_do_mth_f):linspace(nanmin(w_do_mth_f),nanmax(w_do_mth_f),20):nanmax(w_do_mth_f)]);
set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
for i = 1:12
%     i = 1;
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 w_do_mi_p sig_txt *_mi_p *_bc
    for sig = 1:6 
        clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 w_do_mi_p sig_txt *_mi_p recon_dx
        w_do_mi_p=w_do(indx_12mth{i});
        idx1 = find(isnan(w_do_mi_p) == 0);
        sig_txt = sig*std(w_do_mi_p(idx1))
        upper_bc(sig) = mean(w_do_mi_p(idx1)) + sig*std(w_do_mi_p(idx1));
        lower_bc(sig) = mean(w_do_mi_p(idx1)) - sig*std(w_do_mi_p(idx1));
        w_do_mi_p(find(w_do_mi_p >  mean(w_do_mi_p(idx1)) + sig*std(w_do_mi_p(idx1))))=NaN;
        w_do_mi_p(find(w_do_mi_p < mean(w_do_mi_p(idx1)) - sig*std(w_do_mi_p(idx1))))=NaN;
        w_do_mth_f = w_do_mi_p;

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
        b0(sig) = b(1);  b1(sig) = b(2);
        yCalc2{i,sig} = [1:length(w_do_mi_p)] .*b1(sig) + b0(sig);
%         yCalc2{i,sig} = X*b;
%         plot(idx,yCalc2{sig},'--','color',color_spec(sig),'linew',2)
        recon_dx = indx_12mth{i};
        do_recon_mth(sig,recon_dx) = yCalc2{i,sig};
%         plot([1:length(w_do_mi_p)],y_eqa,'--','color','c','linew',2)
        % gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
        % gtext([ num2str(sig) '-sigma = ' num2str(sig_txt)],'Color','k','FontSize',16)
%         line([1:1:length(w_do_mi_p)], upper_bc(sig),'color', color_spec(sig),'linew',2)
%         line([1:1:length(w_do_mi_p)], lower_bc(sig),'color', color_spec(sig),'linew',2)
    if i== 12
    plot(do_recon_mth(sig,:),'--','color', color_spec(sig),'linew',2)
    end
    end
end
ylim([4 inf])
    title(['Namgang-2-regression DO raw month(daily) recon ' num2str(sig) '-sigma'],'fontsize',13)
    print('-dpng', ['Namgang-2-regression_DO_raw_month(daily)_recon_' num2str(sig) '-sigma']);
set(gca,'xtick',[1 t_tick(1:1:end)]);
set(gca,'xlim',[1 t_tick(10)]);
set(gca,'xticklabel',1980:1:2020);    
    print('-dpng', ['Namgang-2-regression_DO_raw_month(daily)_recon_80s' num2str(sig) '-sigma']);
    
%%%%

chl_recon_mth = NaN(6,length(w_chl));
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 w_chl_mi_p sig_txt *_mi_p *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
plot(w_chl,'.','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('chl (mg/L)','fontsize',13)
grid on; 
ylim([0 inf])
% set(gca,'ytick',[nanmin(w_chl_mth_f):linspace(nanmin(w_chl_mth_f),nanmax(w_chl_mth_f),20):nanmax(w_chl_mth_f)]);
set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
for i = 1:12
%     i = 1;
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 w_chl_mi_p sig_txt *_mi_p *_bc
    for sig = 1:6 
        clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 w_chl_mi_p sig_txt *_mi_p recon_dx
        w_chl_mi_p=w_chl(indx_12mth{i});
        idx1 = find(isnan(w_chl_mi_p) == 0);
        sig_txt = sig*std(w_chl_mi_p(idx1))
        upper_bc(sig) = mean(w_chl_mi_p(idx1)) + sig*std(w_chl_mi_p(idx1));
        lower_bc(sig) = mean(w_chl_mi_p(idx1)) - sig*std(w_chl_mi_p(idx1));
        w_chl_mi_p(find(w_chl_mi_p >  mean(w_chl_mi_p(idx1)) + sig*std(w_chl_mi_p(idx1))))=NaN;
        w_chl_mi_p(find(w_chl_mi_p < mean(w_chl_mi_p(idx1)) - sig*std(w_chl_mi_p(idx1))))=NaN;
        w_chl_mth_f = w_chl_mi_p;

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
        b0(sig) = b(1);  b1(sig) = b(2);
        yCalc2{i,sig} = [1:length(w_chl_mi_p)] .*b1(sig) + b0(sig);
%         yCalc2{i,sig} = X*b;
%         plot(idx,yCalc2{sig},'--','color',color_spec(sig),'linew',2)
        recon_dx = indx_12mth{i};
        chl_recon_mth(sig,recon_dx) = yCalc2{i,sig};
%         plot([1:length(w_chl_mi_p)],y_eqa,'--','color','c','linew',2)
        % gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
        % gtext([ num2str(sig) '-sigma = ' num2str(sig_txt)],'Color','k','FontSize',16)
%         line([1:1:length(w_chl_mi_p)], upper_bc(sig),'color', color_spec(sig),'linew',2)
%         line([1:1:length(w_chl_mi_p)], lower_bc(sig),'color', color_spec(sig),'linew',2)
    if i== 12
    plot(chl_recon_mth(sig,:),'--','color', color_spec(sig),'linew',2)
%         plot(chl_recon_mth(sig,:),'--','color', color_spec(sig),'linew',1)
    end
    end
end
ylim([-inf 45])
    title(['Namgang-2-regression chl raw month(daily) recon ' num2str(sig) '-sigma'],'fontsize',13)
    print('-dpng', ['Namgang-2-regression_chl_raw_month(daily)_recon_' num2str(sig) '-sigma']);
set(gca,'xtick',[1 t_tick(1:1:end)]);
set(gca,'xlim',[1 t_tick(10)]);
set(gca,'xticklabel',1980:1:2020);    
    print('-dpng', ['Namgang-2-regression_chl_raw_month(daily)_recon_80s' num2str(sig) '-sigma']);

    
nh4_recon_mth = NaN(6,length(w_nh4));
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 w_nh4_mi_p sig_txt *_mi_p *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
plot(w_nh4,'.','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('nh4 (mg/L)','fontsize',13)
grid on; 
% ylim([0 inf])
% set(gca,'ytick',[nanmin(w_nh4_mth_f):linspace(nanmin(w_nh4_mth_f),nanmax(w_nh4_mth_f),20):nanmax(w_nh4_mth_f)]);
set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
for i = 1:12
%     i = 1;
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 w_nh4_mi_p sig_txt *_mi_p *_bc
    for sig = 1:6 
        clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 w_nh4_mi_p sig_txt *_mi_p recon_dx
        w_nh4_mi_p=w_nh4(indx_12mth{i});
        idx1 = find(isnan(w_nh4_mi_p) == 0);
        sig_txt = sig*std(w_nh4_mi_p(idx1))
        upper_bc(sig) = mean(w_nh4_mi_p(idx1)) + sig*std(w_nh4_mi_p(idx1));
        lower_bc(sig) = mean(w_nh4_mi_p(idx1)) - sig*std(w_nh4_mi_p(idx1));
        w_nh4_mi_p(find(w_nh4_mi_p >  mean(w_nh4_mi_p(idx1)) + sig*std(w_nh4_mi_p(idx1))))=NaN;
        w_nh4_mi_p(find(w_nh4_mi_p < mean(w_nh4_mi_p(idx1)) - sig*std(w_nh4_mi_p(idx1))))=NaN;
        w_nh4_mth_f = w_nh4_mi_p;

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
        b0(sig) = b(1);  b1(sig) = b(2);
        yCalc2{i,sig} = [1:length(w_nh4_mi_p)] .*b1(sig) + b0(sig);
%         yCalc2{i,sig} = X*b;
%         plot(idx,yCalc2{sig},'--','color',color_spec(sig),'linew',2)
        recon_dx = indx_12mth{i};
        nh4_recon_mth(sig,recon_dx) = yCalc2{i,sig};
%         plot([1:length(w_nh4_mi_p)],y_eqa,'--','color','c','linew',2)
        % gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
        % gtext([ num2str(sig) '-sigma = ' num2str(sig_txt)],'Color','k','FontSize',16)
%         line([1:1:length(w_nh4_mi_p)], upper_bc(sig),'color', color_spec(sig),'linew',2)
%         line([1:1:length(w_nh4_mi_p)], lower_bc(sig),'color', color_spec(sig),'linew',2)
    if i== 12
    plot(nh4_recon_mth(sig,:),'--','color', color_spec(sig),'linew',2)
%         plot(nh4_recon_mth(sig,:),'--','color', color_spec(sig),'linew',1)
    end
    end
end
ylim([-0.1 0.5])
    title(['Namgang-2-regression nh4 raw month(daily) recon ' num2str(sig) '-sigma'],'fontsize',13)
    print('-dpng', ['Namgang-2-regression_nh4_raw_month(daily)_recon_' num2str(sig) '-sigma']);
set(gca,'xtick',[1 t_tick(1:1:end)]);
set(gca,'xlim',[1 t_tick(10)]);
set(gca,'xticklabel',1980:1:2020);    
    print('-dpng', ['Namgang-2-regression_nh4_raw_month(daily)_recon_80s' num2str(sig) '-sigma']);    

no3_recon_mth = NaN(6,length(w_no3));
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 w_no3_mi_p sig_txt *_mi_p *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
figure; hold on;
plot(w_no3,'.','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('no3 (mg/L)','fontsize',13)
grid on; 
% ylim([0 inf])
% set(gca,'ytick',[nanmin(w_no3_mth_f):linspace(nanmin(w_no3_mth_f),nanmax(w_no3_mth_f),20):nanmax(w_no3_mth_f)]);
set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
for i = 1:12
%     i = 1;
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 w_no3_mi_p sig_txt *_mi_p *_bc
    for sig = 1:6 
        clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 w_no3_mi_p sig_txt *_mi_p recon_dx
        w_no3_mi_p=w_no3(indx_12mth{i});
        idx1 = find(isnan(w_no3_mi_p) == 0);
        sig_txt = sig*std(w_no3_mi_p(idx1))
        upper_bc(sig) = mean(w_no3_mi_p(idx1)) + sig*std(w_no3_mi_p(idx1));
        lower_bc(sig) = mean(w_no3_mi_p(idx1)) - sig*std(w_no3_mi_p(idx1));
        w_no3_mi_p(find(w_no3_mi_p >  mean(w_no3_mi_p(idx1)) + sig*std(w_no3_mi_p(idx1))))=NaN;
        w_no3_mi_p(find(w_no3_mi_p < mean(w_no3_mi_p(idx1)) - sig*std(w_no3_mi_p(idx1))))=NaN;
        w_no3_mth_f = w_no3_mi_p;

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
        b0(sig) = b(1);  b1(sig) = b(2);
        yCalc2{i,sig} = [1:length(w_no3_mi_p)] .*b1(sig) + b0(sig);
%         yCalc2{i,sig} = X*b;
%         plot(idx,yCalc2{sig},'--','color',color_spec(sig),'linew',2)
        recon_dx = indx_12mth{i};
        no3_recon_mth(sig,recon_dx) = yCalc2{i,sig};
%         plot([1:length(w_no3_mi_p)],y_eqa,'--','color','c','linew',2)
        % gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
        % gtext([ num2str(sig) '-sigma = ' num2str(sig_txt)],'Color','k','FontSize',16)
%         line([1:1:length(w_no3_mi_p)], upper_bc(sig),'color', color_spec(sig),'linew',2)
%         line([1:1:length(w_no3_mi_p)], lower_bc(sig),'color', color_spec(sig),'linew',2)
    if i== 12
    plot(no3_recon_mth(sig,:),'--','color', color_spec(sig),'linew',2)
%         plot(no3_recon_mth(sig,:),'--','color', color_spec(sig),'linew',1)
    end
    end
end
% ylim([-0.1 0.5])
    title(['Namgang-2-regression no3 raw month(daily) recon ' num2str(sig) '-sigma'],'fontsize',13)
    print('-dpng', ['Namgang-2-regression_no3_raw_month(daily)_recon_' num2str(sig) '-sigma']);
set(gca,'xtick',[1 t_tick(1:1:end)]);
set(gca,'xlim',[1 t_tick(10)]);
set(gca,'xticklabel',1980:1:2020);    
    print('-dpng', ['Namgang-2-regression_no3_raw_month(daily)_recon_80s' num2str(sig) '-sigma']);
    
%%%%%%%%%%%%%%%%%%
    do_recon_mth = NaN(1,length(w_do));
for i = 1:12
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 w_do_mi_p sig_txt recon_dx
w_do_mi_p=w_do_mi(indx_12mth{i});
idx1 = find(isnan(w_do_mi_p) == 0);
sig_txt = 3*std(w_do_mi_p(idx1))
w_do_mi_p(find(w_do_mi_p > sig*std(w_do_mi_p(idx1))))=NaN;
w_do_mi_p(find(w_do_mi_p < -sig*std(w_do_mi_p(idx1))))=NaN;
w_do_mth_f = w_do_mi_p;
figure; hold on;
plot(w_do_mth_f,'.','linew',2);
xlabel('time(year)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title(['sumjin(hadong)-regression DO raw-climate ' num2str(i) '-month(daily) ' num2str(sig) '-sigma'],'fontsize',13)
grid on; 
ylim([-10 5])
set(gca,'ytick',[-10:2:5]);
set(gca,'fontsize',13)
set(gca,'xtick',[1:31*3:length(w_do_mth_f)]);
set(gca,'xlim',[1 length(w_do_mth_f)]);
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
yCalc2{i} = X*b;
plot(idx,yCalc2{i},'--','color','r','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
gtext(['3-sigma = ' num2str(sig_txt)],'Color','k','FontSize',16)
y_eqa = b(2).*[1:length(w_do_mi_p)] + b(1); 
% print('-dpng', ['sumjin(hadong)-regression DO_raw-climate ' num2str(i) '-month(daily) ' num2str(sig) '-sigma']);
recon_dx = indx_12mth{i};
do_recon_mth(recon_dx) = y_eqa;
plot([1:length(w_do_mi_p)],y_eqa,'--','color','c','linew',2)
end

chl_recon_mth = NaN(1,length(w_chl));
for i = 1:12
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt w_chl_mi_p recon_dx
% raw
w_chl_mi_p=w_chl_mi(indx_12mth{i});
idx1 = find(isnan(w_chl_mi_p) == 0);
sig_txt = 3*std(w_chl_mi_p(idx1));
w_chl_mi_p(find(w_chl_mi_p > sig_txt))=NaN;
w_chl_mi_p(find(w_chl_mi_p < -sig_txt))=NaN;
w_chl_mth_f = w_chl_mi_p;
 figure; hold on;
  plot(w_chl_mth_f,'.','linew',2)
  xlabel('time(year)','fontsize',13)
ylabel('Chla (mg/m^3)','fontsize',13)
title(['sumjin(hadong)-regression Chla raw-climate ' num2str(i) '-month(daily) '  num2str(sig) '-sigma'],'fontsize',13)
grid on
set(gca,'fontsize',13)
  ylim([-10 30])
set(gca,'xtick',[1:31*3:length(w_chl_mth_f)]);
set(gca,'xlim',[1 length(w_chl_mth_f)]);
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
gtext(['3-sigma = ' num2str(sig_txt)],'Color','k','FontSize',16)
% print('-dpng', ['sumjin(hadong)-regression Chla_raw-climate ' num2str(i) '-month(daily) '  num2str(sig) '-sigma']);
y_eqa = b(2).*[1:length(w_chl_mi_p)] + b(1); 
recon_dx = indx_12mth{i};
chl_recon_mth(recon_dx) = y_eqa;
plot([1:length(w_chl_mi_p)],y_eqa,'--','color','c','linew',2)
end
 
for i = 1:12
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_nh4_mth_f = w_nh4_mi(indx_12mth{i});
     figure; hold on;
  plot(w_nh4_mth_f,'.','linew',2)
xlabel('time(year)','fontsize',13)
ylabel('NH4 (mg/L)','fontsize',13)
title(['Namgang-2-regression NH4 raw-climate ' num2str(i) '-month(daily)'],'fontsize',13)
grid on
set(gca,'fontsize',13)
    ylim([-0.2 2])
set(gca,'xtick',[1:31*3:length(w_nh4_mth_f)]);
set(gca,'xlim',[1 length(w_nh4_mth_f)]);
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
print('-dpng', ['Namgang-2-regression NH4_raw-climate ' num2str(i) '-month(daily)']);
end
 
for i = 1:12
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f
w_no3_mth_f = w_no3_mi(indx_12mth{i}); 
figure; hold on;
  plot(w_no3_mth_f,'.','linew',2)
  xlabel('time(year)','fontsize',13)
ylabel('NO3 (mg/L)','fontsize',13)
title(['Namgang-2-regression NO3 raw-climate ' num2str(i) '-month(daily)'],'fontsize',13)
grid on
set(gca,'fontsize',13)
set(gca,'xtick',[1:31*3:length(w_no3_mth_f)]);
set(gca,'xlim',[1 length(w_no3_mth_f)]);
set(gca,'xticklabel',1989:3:2019);
      ylim([-2 4])
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
print('-dpng', ['Namgang-2-regression NO3_raw-climate ' num2str(i) '-month(daily)']);
end
