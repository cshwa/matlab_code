close all; clear; clc;
% -------------------------------------------------------------------------
t_year = 1989:2018
% _________________________________________________________________________
% cd D:\장기생태\Dynamic\06_river
% load jinwall_raw_data_daily.mat %there has to be song vs. song_transp
%jinwall vs. jinwall_transp
% load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','raw_*')
cd D:\장기생태\Dynamic\06_river\환경과학원
% load songjung_yymm_koem_data_89to19_3sig.mat  % polynomial_fitting_on_ecological_variable_to2004_new_v6_3sigma.m
load songjung_yymm_monthly_data_89to19_3sig.mat  % polynomial_fitting_on_ecological_variable_to2004_new_v6_3sigma.m
% addP=load('songjung_bio_v3_day_koem_yr.mat'); % 2004_regime_v3.m

cd D:\장기생태\Dynamic\06_river\data\sj_1980to1996
%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
load songjung_discharge_1980to2018.mat  % pre_merg_dis is discharge
% load excel 2019 trans
[raw_sj txt_sj]=xlsread('D:\장기생태\Dynamic\06_river\환경과학원\수리·수문·기상_유량_구례_20190101-20191231.xls','유량');


% cd D:\장기생태\Dynamic\KOEM
% % yj=load(['koem_timeseires_koem_gy_only_1to3sig_v2_each_st.mat']); % full (9point) data (extract_sigma_for_each_st)
% yj=load('koem_timeseires_monthly_gy_only_1to3sig_v2_each_st.mat'); % full (9point) data (extract_sigma_for_each_st)
cd D:\장기생태\Dynamic\06_river
yj=load('koem_timeseires_monthly_gy_only_1to3sig_v2_each_st_16p.mat'); % 16 points 
% *_sur(bot)_clim = raw data koem  (no sigma extract)

cd D:\장기생태\Dynamic\06_river
% approximation : Sumjin(songjung) & Namgang have same discharge during 1980~1989.
load('sumjin_recons_water_temp_present.mat','merg_recon_w_c');
sumjin_re_w_c = merg_recon_w_c;

sj=load('songjung_plt.mat'); % from plot_obs_trans_vs_bio_2_songjung_monthly.m

clearvars temperature sj_trans_out
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
% tempo_temp=sumjin_re_w_c{t_year(order_i)-1989};
tempo_trans = dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
%     temperature=tempo_temp;
    sj_trans_out=tempo_trans;
else
%     temperature=cat(1,temperature,tempo_temp);
    sj_trans_out=cat(1,sj_trans_out,tempo_trans);
end
end
sj_trans_out(end+1:end+365) = raw_sj(:,4);

%% KOEM (yoonja) make it to 1989~2018 form (origin 1997~2019)
si_MW = 28.085530; 
gap_yr = length(1989:2019)*12 - size(yj.regime_temp_3s,2);


koem_temp_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_temp_b_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_salt_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_salt_b_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_do_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_do_b_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_chl_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_chl_b_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_din_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_din_b_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_nh4_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_nh4_b_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_no3_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_no3_b_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_po4_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_po4_b_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_secchi_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_ss_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_ss_b_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_si_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);
koem_si_b_raw=NaN(size(yj.regime_po4_3s,1),length(1989:2019)*12);

% sp mean
koem_temp=NaN(1,length(1989:2019)*12);
koem_temp_b=NaN(1,length(1989:2019)*12);
koem_salt=NaN(1,length(1989:2019)*12);
koem_salt_b=NaN(1,length(1989:2019)*12);
koem_din=NaN(1,length(1989:2019)*12);
koem_din_b=NaN(1,length(1989:2019)*12);
koem_chl=NaN(1,length(1989:2019)*12);
koem_chl_b=NaN(1,length(1989:2019)*12);
koem_nh4=NaN(1,length(1989:2019)*12);
koem_nh4_b=NaN(1,length(1989:2019)*12);
koem_no3=NaN(1,length(1989:2019)*12);
koem_no3_b=NaN(1,length(1989:2019)*12);
koem_po4=NaN(1,length(1989:2019)*12);
koem_po4_b=NaN(1,length(1989:2019)*12);
koem_secchi=NaN(1,length(1989:2019)*12);
koem_ss=NaN(1,length(1989:2019)*12);
koem_ss_b=NaN(1,length(1989:2019)*12);
koem_si=NaN(1,length(1989:2019)*12);
koem_si_b=NaN(1,length(1989:2019)*12);

% 3 sigma extract data
koem_temp_raw(:,gap_yr+1:end)=yj.regime_temp_3s(:,1:end);
koem_temp_b_raw(:,gap_yr+1:end)=yj.regime_temp_b_3s(:,1:end);
koem_salt_raw(:,gap_yr+1:end)=yj.regime_salt_3s(:,1:end);
koem_salt_b_raw(:,gap_yr+1:end)=yj.regime_salt_b_3s(:,1:end);
koem_chl_raw(:,gap_yr+1:end)=yj.regime_chl_3s(:,1:end);
koem_chl_b_raw(:,gap_yr+1:end)=yj.regime_chl_b_3s(:,1:end);
koem_do_raw(:,gap_yr+1:end)=yj.regime_do_3s(:,1:end);
koem_do_b_raw(:,gap_yr+1:end)=yj.regime_do_b_3s(:,1:end);

koem_nh4_raw(:,gap_yr+1:end)=yj.regime_nh4_3s(:,1:end)./ yj.N_MW;
koem_nh4_b_raw(:,gap_yr+1:end)=yj.regime_nh4_b_3s(:,1:end)./ yj.N_MW;
koem_din_raw(:,gap_yr+1:end)=yj.regime_din_3s(:,1:end)./ yj.N_MW;
koem_din_b_raw(:,gap_yr+1:end)=yj.regime_din_b_3s(:,1:end)./ yj.N_MW;
koem_no3_raw(:,gap_yr+1:end)=yj.regime_no3_3s(:,1:end)./ yj.N_MW;
koem_no3_b_raw(:,gap_yr+1:end)=yj.regime_no3_b_3s(:,1:end)./ yj.N_MW;
koem_po4_raw(:,gap_yr+1:end)=yj.regime_po4_3s(:,1:end)./ yj.P_MW;
koem_po4_b_raw(:,gap_yr+1:end)=yj.regime_po4_b_3s(:,1:end)./ yj.P_MW;

koem_secchi_raw(:,gap_yr+1:end)=yj.regime_secchi_3s(:,1:end);
koem_ss_raw(:,gap_yr+1:end)=yj.regime_ss_3s(:,1:end);
koem_ss_b_raw(:,gap_yr+1:end)=yj.regime_ss_b_3s(:,1:end);
koem_si_raw(:,gap_yr+1:end)=yj.regime_si_3s(:,1:end)./ si_MW;
koem_si_b_raw(:,gap_yr+1:end)=yj.regime_si_b_3s(:,1:end)./ si_MW;

koem_temp=nanmean(koem_temp_raw,1);
koem_temp_b=nanmean(koem_temp_b_raw,1);
koem_salt=nanmean(koem_salt_raw,1);
koem_salt_b=nanmean(koem_salt_b_raw,1);
koem_din=nanmean(koem_din_raw,1);
koem_din_b=nanmean(koem_din_b_raw,1);
koem_chl=nanmean(koem_chl_raw,1);
koem_chl_b=nanmean(koem_chl_b_raw,1);
koem_nh4=nanmean(koem_nh4_raw,1);
koem_nh4_b=nanmean(koem_nh4_b_raw,1);
koem_no3=nanmean(koem_no3_raw,1);
koem_no3_b=nanmean(koem_no3_b_raw,1);
koem_po4=nanmean(koem_po4_raw,1);
koem_po4_b=nanmean(koem_po4_b_raw,1);
koem_secchi=nanmean(koem_secchi_raw,1);
koem_ss=nanmean(koem_ss_raw,1);
koem_ss_b=nanmean(koem_ss_b_raw,1);
koem_si=nanmean(koem_si_raw,1);
koem_si_b=nanmean(koem_si_b_raw,1);

%% regime shif detection data
% % % spmean 
% clearvars koem_chl koem_chl_b
% koem_chl=nanmean(yj.regime_chl_3s(sp_9p_st,:),1)';
% koem_chl_b=nanmean(yj.regime_chl_b_3s(sp_9p_st,:),1)';
% nandx=find(isnan(koem_chl)==1);
% nandx_b=find(isnan(koem_chl_b)==1);
% in_date = yj.ref_date;
% in_date_b = yj.ref_date;
% koem_chl(nandx)=[];
% koem_chl_b(nandx_b)=[];
% in_date(nandx)=[];
% in_date_b(nandx_b)=[];
% 
% % % st
% clearvars koem_chl koem_chl_b
% i=9
% koem_chl=yj.regime_chl_3s(sp_9p_st(i),:)';
% koem_chl_b=yj.regime_chl_b_3s(sp_9p_st(i),:)';
% nandx=find(isnan(koem_chl)==1);
% nandx_b=find(isnan(koem_chl_b)==1);
% in_date = yj.ref_date;
% in_date_b = yj.ref_date;
% koem_chl(nandx)=[];
% koem_chl_b(nandx_b)=[];
% in_date(nandx)=[];
% in_date_b(nandx_b)=[];

%%

% no sigma extraction
% koem_temp_raw(:,gap_yr+1:end)=yj.temp_sur_clim(:,1:end);
% koem_temp_b_raw(:,gap_yr+1:end)=yj.temp_bot_clim(:,1:end);
% koem_chl_raw(:,gap_yr+1:end)=yj.chl_sur_clim(:,1:end);
% koem_chl_b_raw(:,gap_yr+1:end)=yj.chl_bot_clim(:,1:end);
% koem_nh4_raw(:,gap_yr+1:end)=yj.nh4_sur_clim(:,1:end)./ yj.N_MW;
% koem_nh4_b_raw(:,gap_yr+1:end)=yj.nh4_bot_clim(:,1:end)./ yj.N_MW;
% koem_no3_raw(:,gap_yr+1:end)=yj.no3_sur_clim(:,1:end)./ yj.N_MW;
% koem_no3_b_raw(:,gap_yr+1:end)=yj.no3_bot_clim(:,1:end)./ yj.N_MW;
% koem_din_raw(:,gap_yr+1:end)=yj.din_sur_clim(:,1:end)./ yj.N_MW;
% koem_din_b_raw(:,gap_yr+1:end)=yj.din_bot_clim(:,1:end)./ yj.N_MW;
% koem_po4_raw(:,gap_yr+1:end)=yj.po4_sur_clim(:,1:end)./ yj.P_MW;
% koem_po4_b_raw(:,gap_yr+1:end)=yj.po4_bot_clim(:,1:end)./ yj.P_MW;
% 
% koem_secchi_raw(:,gap_yr+1:end)=yj.secchi_sur_clim(:,1:end);
% koem_ss_raw(:,gap_yr+1:end)=yj.ss_sur_clim(:,1:end);
% koem_ss_b_raw(:,gap_yr+1:end)=yj.ss_bot_clim(:,1:end);
% koem_si_raw(:,gap_yr+1:end)=yj.si_sur_clim(:,1:end)./ si_MW;
% koem_si_b_raw(:,gap_yr+1:end)=yj.si_bot_clim(:,1:end)./ si_MW;
% 
% koem_temp=nanmean(koem_temp_raw,1);
% koem_temp_b=nanmean(koem_temp_b_raw,1);
% koem_chl=nanmean(koem_chl_raw,1);
% koem_chl_b=nanmean(koem_chl_b_raw,1);
% koem_nh4=nanmean(koem_nh4_raw,1);
% koem_nh4_b=nanmean(koem_nh4_b_raw,1);
% koem_no3=nanmean(koem_no3_raw,1);
% koem_no3_b=nanmean(koem_no3_b_raw,1);
% koem_po4=nanmean(koem_po4_raw,1);
% koem_po4_b=nanmean(koem_po4_b_raw,1);
% koem_secchi=nanmean(koem_secchi_raw,1);
% koem_ss=nanmean(koem_ss_raw,1);
% koem_ss_b=nanmean(koem_ss_b_raw,1);
% koem_si=nanmean(koem_si_raw,1);
% koem_si_b=nanmean(koem_si_b_raw,1);

% koem_chl=nanmean(koem_chl_raw,1);
% koem_chl_b=nanmean(koem_chl_b_raw,1);
% koem_nh4=nanmean(koem_nh4_raw,1);
% koem_nh4_b=nanmean(koem_nh4_b_raw,1);
% koem_no3=nanmean(koem_no3_raw,1);
% koem_no3_b=nanmean(koem_no3_b_raw,1);
% koem_po4=nanmean(koem_po4_raw,1);
% koem_po4_b=nanmean(koem_po4_b_raw,1);


% fill missing value
koem_in_chl = koem_chl;
koem_in_no3 = koem_no3;
koem_in_nh4 = koem_nh4;
koem_in_po4 = koem_po4;
koem_in_secchi = koem_secchi;
koem_in_ss = koem_ss;
koem_in_si = koem_si;

koem_in_chl_b = koem_chl_b;
koem_in_no3_b = koem_no3_b;
koem_in_nh4_b = koem_nh4_b;
koem_in_po4_b = koem_po4_b;
koem_in_ss_b = koem_ss_b;
koem_in_si_b = koem_si_b;

t=1:length(koem_chl);
koem_in_chl(isnan(koem_in_chl))=interp1(t(~isnan(koem_in_chl)),koem_in_chl(~isnan(koem_in_chl)),t(isnan(koem_in_chl)));
koem_in_no3(isnan(koem_in_no3))=interp1(t(~isnan(koem_in_no3)),koem_in_no3(~isnan(koem_in_no3)),t(isnan(koem_in_no3)));
koem_in_nh4(isnan(koem_in_nh4))=interp1(t(~isnan(koem_in_nh4)),koem_in_nh4(~isnan(koem_in_nh4)),t(isnan(koem_in_nh4)));
koem_in_po4(isnan(koem_in_po4))=interp1(t(~isnan(koem_in_po4)),koem_in_po4(~isnan(koem_in_po4)),t(isnan(koem_in_po4)));
koem_in_secchi(isnan(koem_in_secchi))=interp1(t(~isnan(koem_in_secchi)),koem_in_secchi(~isnan(koem_in_secchi)),t(isnan(koem_in_secchi)));
koem_in_ss(isnan(koem_in_ss))=interp1(t(~isnan(koem_in_ss)),koem_in_ss(~isnan(koem_in_ss)),t(isnan(koem_in_ss)));
koem_in_si(isnan(koem_in_si))=interp1(t(~isnan(koem_in_si)),koem_in_si(~isnan(koem_in_si)),t(isnan(koem_in_si)));

koem_in_chl_b(isnan(koem_in_chl_b))=interp1(t(~isnan(koem_in_chl_b)),koem_in_chl_b(~isnan(koem_in_chl_b)),t(isnan(koem_in_chl_b)));
koem_in_no3_b(isnan(koem_in_no3_b))=interp1(t(~isnan(koem_in_no3_b)),koem_in_no3_b(~isnan(koem_in_no3_b)),t(isnan(koem_in_no3_b)));
koem_in_nh4_b(isnan(koem_in_nh4_b))=interp1(t(~isnan(koem_in_nh4_b)),koem_in_nh4_b(~isnan(koem_in_nh4_b)),t(isnan(koem_in_nh4_b)));
koem_in_po4_b(isnan(koem_in_po4_b))=interp1(t(~isnan(koem_in_po4_b)),koem_in_po4_b(~isnan(koem_in_po4_b)),t(isnan(koem_in_po4_b)));
koem_in_ss_b(isnan(koem_in_ss_b))=interp1(t(~isnan(koem_in_ss_b)),koem_in_ss_b(~isnan(koem_in_ss_b)),t(isnan(koem_in_ss_b)));
koem_in_si_b(isnan(koem_in_si_b))=interp1(t(~isnan(koem_in_si_b)),koem_in_si_b(~isnan(koem_in_si_b)),t(isnan(koem_in_si_b)));

% 2019 has to be remove to fit the term for discharge (1989:2018)
% koem_chl(361:end) = [];
% koem_temp(361:end) = [];
% koem_no3(361:end) = [];
% koem_nh4(361:end) = [];
% koem_po4(361:end) = [];
% koem_secchi(361:end) = [];
% koem_ss(361:end) = [];
% koem_si(361:end) = [];
% 
% koem_chl_b(361:end) = [];
% koem_temp_b(361:end) = [];
% koem_no3_b(361:end) = [];
% koem_nh4_b(361:end) = [];
% koem_po4_b(361:end) = [];
% koem_ss_b(361:end) = [];
% koem_si_b(361:end) = [];
% 
% koem_in_chl(361:end) = [];
% koem_in_no3(361:end) = [];
% koem_in_nh4(361:end) = [];
% koem_in_po4(361:end) = [];
% koem_in_secchi(361:end) = [];
% koem_in_ss(361:end) = [];
% koem_in_si(361:end) = [];
% 
% koem_in_chl_b(361:end) = [];
% koem_in_no3_b(361:end) = [];
% koem_in_nh4_b(361:end) = [];
% koem_in_po4_b(361:end) = [];
% koem_in_ss_b(361:end) = [];
% koem_in_si_b(361:end) = [];

%% each st. YEARLY
clearvars name_tag
name_tag{1}={'광양항01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['광양만',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['여수연안',num2str(k,'%02d')]};
end

%% 9 point
sp_9p_st = [1:6, 14:16];

%% monthly mean (yy.mm) form after over 3sigma value extraction.

for i=1:length(t_indx)% time for num. of months (until 2018.12.31)
if i ==1
   monthly_trans(i) = nanmean(sj_trans_out(1:t_indx(i)));
else
   monthly_trans(i) = nanmean(sj_trans_out(t_indx(i-1)+1:t_indx(i)));
end
end

%% pick specific year of transp.
% 1989:2019
% pick_yr=[20,25,27,28,29];
% 2008,2013,2015,2016,2017
for i=1:length(pick_yr)% time for num. of months (until 2018.12.31)
        pick_yr_trans{i} = sj_trans_out(t_indx((pick_yr(i)-1)*12)+1:t_indx((pick_yr(i))*12));
end

figure; hold on;
for i=1:17*12 % time for num. of months (until 2018.12.31)
if i ==1
   plot(sj_trans_out(1:t_indx(i)),'color',[.5 .5 .5]);
else
   plot(sj_trans_out(t_indx(i-1)+1:t_indx(i)),'color',[.5 .5 .5]);
end
end
figure
hold on;
for i=17*12+1:27*12 % time for num. of months (until 2018.12.31)
if i ==1
   plot(sj_trans_out(1:t_indx(i)),'color','g');
else
   plot(sj_trans_out(t_indx(i-1)+1:t_indx(i)),'color','g');
end
end

figure
hold on;
for i=27*12+1:31*12 % time for num. of months (until 2018.12.31)
if i ==1
   plot(sj_trans_out(1:t_indx(i)),'color','r');
else
   plot(sj_trans_out(t_indx(i-1)+1:t_indx(i)),'color','r');
end
end

% t_indx(end) = 10957

figure;
plot(sj_trans_out); hold on
line(1:5844,repmat(mean(sj_trans_out(1:5844)),5844,1),'color','r') % ~2004
line(5845:9861,repmat(mean(sj_trans_out(5845:9861)),length(5845:9861),1),'color','g')% ~2015
line(9862:length(sj_trans_out),repmat(mean(sj_trans_out(9862:length(sj_trans_out))),length(9862:length(sj_trans_out)),1),'color','m')% % 2016 ~ 

figure; plot(monthly_trans)
% hold on; for i = 6:12:360; xline(i,'color','r'); end

for i = 1:31 %year
yr_trans(i) = nanmean(monthly_trans((i-1)*12+1:12*i));
yr_temp(i) = nanmean(koem_temp((i-1)*12+1:12*i));
yr_salt(i) = nanmean(koem_salt((i-1)*12+1:12*i));
yr_chl(i) = nanmean(koem_chl((i-1)*12+1:12*i));
yr_nh4(i) = nanmean(koem_nh4((i-1)*12+1:12*i));
yr_no3(i) = nanmean(koem_no3((i-1)*12+1:12*i));
yr_po4(i) = nanmean(koem_po4((i-1)*12+1:12*i));
yr_secchi(i) = nanmean(koem_secchi((i-1)*12+1:12*i));
yr_ss(i) = nanmean(koem_ss((i-1)*12+1:12*i));
yr_si(i) = nanmean(koem_si((i-1)*12+1:12*i));
yr_in_chl(i) = nanmean(koem_in_chl((i-1)*12+1:12*i));
yr_in_nh4(i) = nanmean(koem_in_nh4((i-1)*12+1:12*i));
yr_in_no3(i) = nanmean(koem_in_no3((i-1)*12+1:12*i));
yr_in_po4(i) = nanmean(koem_in_po4((i-1)*12+1:12*i));
yr_in_secchi(i) = nanmean(koem_in_secchi((i-1)*12+1:12*i));
yr_in_ss(i) = nanmean(koem_in_ss((i-1)*12+1:12*i));
yr_in_si(i) = nanmean(koem_in_si((i-1)*12+1:12*i));

yr_temp_b(i) = nanmean(koem_temp_b((i-1)*12+1:12*i));
yr_salt_b(i) = nanmean(koem_salt_b((i-1)*12+1:12*i));
yr_chl_b(i) = nanmean(koem_chl_b((i-1)*12+1:12*i));
yr_nh4_b(i) = nanmean(koem_nh4_b((i-1)*12+1:12*i));
yr_no3_b(i) = nanmean(koem_no3_b((i-1)*12+1:12*i));
yr_po4_b(i) = nanmean(koem_po4_b((i-1)*12+1:12*i));
yr_ss_b(i) = nanmean(koem_ss_b((i-1)*12+1:12*i));
yr_si_b(i) = nanmean(koem_si_b((i-1)*12+1:12*i));
yr_in_chl_b(i) = nanmean(koem_in_chl_b((i-1)*12+1:12*i));
yr_in_nh4_b(i) = nanmean(koem_in_nh4_b((i-1)*12+1:12*i));
yr_in_no3_b(i) = nanmean(koem_in_no3_b((i-1)*12+1:12*i));
yr_in_po4_b(i) = nanmean(koem_in_po4_b((i-1)*12+1:12*i));
yr_in_ss_b(i) = nanmean(koem_in_ss_b((i-1)*12+1:12*i));
yr_in_si_b(i) = nanmean(koem_in_si_b((i-1)*12+1:12*i));
end

figure;
plot(yr_trans); hold on
line(1:17,repmat(mean(yr_trans(1:17)),17,1),'color','r') % ~2004
line(18:27,repmat(mean(yr_trans(18:27)),length(18:27),1),'color','g')% ~2015
line(28:31,repmat(mean(yr_trans(28:31)),length(28:31),1),'color','m')% % 2016 ~ 

for i = 1:31 %year
yr_temp_9p(i) = nanmean(nanmean(koem_temp_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_do_9p(i) = nanmean(nanmean(koem_do_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_salt_9p(i) = nanmean(nanmean(koem_salt_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_chl_9p(i) = nanmean(nanmean(koem_chl_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_din_9p(i) = nanmean(nanmean(koem_din_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_nh4_9p(i) = nanmean(nanmean(koem_nh4_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_no3_9p(i) = nanmean(nanmean(koem_no3_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_po4_9p(i) = nanmean(nanmean(koem_po4_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_secchi_9p(i) = nanmean(nanmean(koem_secchi_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_ss_9p(i) = nanmean(nanmean(koem_ss_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_si_9p(i) = nanmean(nanmean(koem_si_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_chl_9p(i) = nanmean(nanmean(koem_in_chl_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_nh4_9p(i) = nanmean(nanmean(koem_in_nh4_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_no3_9p(i) = nanmean(nanmean(koem_in_no3_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_po4_9p(i) = nanmean(nanmean(koem_in_po4_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_secchi_9p(i) = nanmean(nanmean(koem_in_secchi_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_ss_9p(i) = nanmean(nanmean(koem_in_ss_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_si_9p(i) = nanmean(nanmean(koem_in_si_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);

yr_temp_b_9p(i) = nanmean(nanmean(koem_temp_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_do_b_9p(i) = nanmean(nanmean(koem_do_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_salt_b_9p(i) = nanmean(nanmean(koem_salt_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_chl_b_9p(i) = nanmean(nanmean(koem_chl_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_din_b_9p(i) = nanmean(nanmean(koem_din_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_nh4_b_9p(i) = nanmean(nanmean(koem_nh4_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_no3_b_9p(i) = nanmean(nanmean(koem_no3_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_po4_b_9p(i) = nanmean(nanmean(koem_po4_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_ss_b_9p(i) = nanmean(nanmean(koem_ss_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
yr_si_b_9p(i) = nanmean(nanmean(koem_si_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_chl_b_9p(i) = nanmean(nanmean(koem_in_chl_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_nh4_b_9p(i) = nanmean(nanmean(koem_in_nh4_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_no3_b_9p(i) = nanmean(nanmean(koem_in_no3_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_po4_b_9p(i) = nanmean(nanmean(koem_in_po4_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_ss_b_9p(i) = nanmean(nanmean(koem_in_ss_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
% yr_in_si_b_9p(i) = nanmean(nanmean(koem_in_si_b_raw(sp_9p_st,(i-1)*12+1:12*i),1),2);
end

for i = 1:31
   yr_salt_raw(:,i) = nanmean(koem_salt_raw(:,(i-1)*12+1:12*i),2);
   yr_salt_b_raw(:,i) = nanmean(koem_salt_b_raw(:,(i-1)*12+1:12*i),2);
   yr_secchi_raw(:,i) = nanmean(koem_secchi_raw(:,(i-1)*12+1:12*i),2);
   yr_chl_raw(:,i) = nanmean(koem_chl_raw(:,(i-1)*12+1:12*i),2);
   yr_temp_raw(:,i) = nanmean(koem_temp_raw(:,(i-1)*12+1:12*i),2);
   yr_temp_b_raw(:,i) = nanmean(koem_temp_b_raw(:,(i-1)*12+1:12*i),2);
   yr_nh4_raw(:,i) = nanmean(koem_nh4_raw(:,(i-1)*12+1:12*i),2); 
   yr_no3_raw(:,i) = nanmean(koem_no3_raw(:,(i-1)*12+1:12*i),2); 
   yr_nh4_b_raw(:,i) = nanmean(koem_nh4_b_raw(:,(i-1)*12+1:12*i),2); 
   yr_no3_b_raw(:,i) = nanmean(koem_no3_b_raw(:,(i-1)*12+1:12*i),2); 
   yr_din_raw(:,i) = nanmean(koem_din_raw(:,(i-1)*12+1:12*i),2); 
   yr_din_b_raw(:,i) = nanmean(koem_din_b_raw(:,(i-1)*12+1:12*i),2); 
   yr_po4_raw(:,i) = nanmean(koem_po4_raw(:,(i-1)*12+1:12*i),2); 
   yr_si_raw(:,i) = nanmean(koem_si_raw(:,(i-1)*12+1:12*i),2); 
   yr_ss_raw(:,i) = nanmean(koem_ss_raw(:,(i-1)*12+1:12*i),2); 
   yr_ss_b_raw(:,i) = nanmean(koem_ss_b_raw(:,(i-1)*12+1:12*i),2); 
   yr_si_b_raw(:,i) = nanmean(koem_si_b_raw(:,(i-1)*12+1:12*i),2); 
end

%% 6~8

for i = 1:31 %year
summer_trans(i) = nanmean(monthly_trans((i-1)*12+1+5:(i-1)*12+1+7));
summer_chl(i) = nanmean(koem_chl((i-1)*12+1+5:(i-1)*12+1+7));
summer_nh4(i) = nanmean(koem_nh4((i-1)*12+1+5:(i-1)*12+1+7));
summer_no3(i) = nanmean(koem_no3((i-1)*12+1+5:(i-1)*12+1+7));
summer_po4(i) = nanmean(koem_po4((i-1)*12+1+5:(i-1)*12+1+7));
summer_secchi(i) = nanmean(koem_secchi((i-1)*12+1+5:(i-1)*12+1+7));
summer_ss(i) = nanmean(koem_ss((i-1)*12+1+5:(i-1)*12+1+7));
summer_si(i) = nanmean(koem_si((i-1)*12+1+5:(i-1)*12+1+7));

summer_in_chl(i) = nanmean(koem_in_chl((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_nh4(i) = nanmean(koem_in_nh4((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_no3(i) = nanmean(koem_in_no3((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_po4(i) = nanmean(koem_in_po4((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_secchi(i) = nanmean(koem_in_secchi((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_ss(i) = nanmean(koem_in_ss((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_si(i) = nanmean(koem_in_si((i-1)*12+1+5:(i-1)*12+1+7));

summer_chl_b(i) = nanmean(koem_chl_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_nh4_b(i) = nanmean(koem_nh4_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_no3_b(i) = nanmean(koem_no3_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_po4_b(i) = nanmean(koem_po4_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_ss_b(i) = nanmean(koem_ss_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_si_b(i) = nanmean(koem_si_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_chl_b(i) = nanmean(koem_in_chl_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_nh4_b(i) = nanmean(koem_in_nh4_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_no3_b(i) = nanmean(koem_in_no3_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_po4_b(i) = nanmean(koem_in_po4_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_ss_b(i) = nanmean(koem_in_ss_b((i-1)*12+1+5:(i-1)*12+1+7));
summer_in_si_b(i) = nanmean(koem_in_si_b((i-1)*12+1+5:(i-1)*12+1+7));
end

for i = 1:31
   summer_secchi_raw(:,i) = nanmean(koem_secchi_raw(:,(i-1)*12+1+5:(i-1)*12+1+7),2);
   summer_ss_raw(:,i) = nanmean(koem_ss_raw(:,(i-1)*12+1+5:(i-1)*12+1+7),2); 
   summer_si_raw(:,i) = nanmean(koem_si_raw(:,(i-1)*12+1+5:(i-1)*12+1+7),2); 
   summer_ss_b_raw(:,i) = nanmean(koem_ss_b_raw(:,(i-1)*12+1+5:(i-1)*12+1+7),2); 
   summer_si_b_raw(:,i) = nanmean(koem_si_b_raw(:,(i-1)*12+1+5:(i-1)*12+1+7),2); 
end

% save('plot_yr_trans_vs_koem.mat'); % 9 points
save('plot_yr_trans_vs_koem_16points.mat'); % 16 points

% cd D:\장기생태\Dynamic\06_river\환경과학원
% load 2007_2019_time_axis.mat
% t_tick(2:end) = t_tick(2:end)+1;
t_tick = 1:12:372 ; %1989.01.01~2020.01.01




% 2008,2013,2015,2016,2017
fig = figure; hold on;
for i = 1:length(pick_yr_trans)
% plot(pick_yr_trans{i},'b*','linew',2);
plot(pick_yr_trans{i},'linew',2);
        title(['river discharge']);
        xlabel('time','fontsize',13)
        ylabel('discharge (m/s)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
%         set(gca,'xtick',1:31);
        legend("2008","2013","2015","2016","2017");
        xlim([1 150])
%         ylim([0 1.2])
%         set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
% yyaxis right
% plot(yr_chl,'g-','linew',2);
% ylim([0 10]);
% xtickangle(45)

end

%% YEARLY
%% Kd vs. CHL
plt_kd = 1.7 ./yr_secchi;
corrcoef(yr_chl(~isnan(yr_chl)),plt_kd(~isnan(yr_chl)))
fig = figure; hold on;
        plot(plt_kd,'b*','linew',2);
        plot(plt_kd,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. Kd']);
        xlabel('time','fontsize',13)
        ylabel('Kd (m^-1)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
ylim([0 10]);
xtickangle(45)

plt_kd = 1.7 ./yr_secchi_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. Kd']);
        xlabel('time','fontsize',13)
        ylabel('Kd (m^-1)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl_9p,'b*','linew',2);
plot(yr_chl_9p,'b-','linew',2);
plot(1:18,repmat(nanmean(yr_chl_9p(1:18)),18,1),':','color','b','linew',2);
plot(19:27,repmat(nanmean(yr_chl_9p(19:27)),length(19:27),1),':','color','b','linew',2);
plot(28:31,repmat(nanmean(yr_chl_9p(28:31)),length(28:31),1),':','color','b','linew',2);
ylabel('Chl (ug/L)');
ylim([0 7]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)


plt_kd = yr_secchi_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. secchi depth']);
        xlabel('time','fontsize',13)
        ylabel('secchi depth (m)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl_9p,'b*','linew',2);
plot(yr_chl_9p,'b-','linew',2);
plot(1:18,repmat(nanmean(yr_chl_9p(1:18)),18,1),':','color','b','linew',2);
plot(19:27,repmat(nanmean(yr_chl_9p(19:27)),length(19:27),1),':','color','b','linew',2);
plot(28:31,repmat(nanmean(yr_chl_9p(28:31)),length(28:31),1),':','color','b','linew',2);
ylabel('Chl (ug/L)');
ylim([0 7]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)
title(['KOEM yearly OBS secchi depth']);

plt_kd = yr_nh4_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. NH4-N']);
        xlabel('time','fontsize',13)
        ylabel('NH4-N (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl_9p,'b*','linew',2);
plot(yr_chl_9p,'b-','linew',2);
plot(1:18,repmat(nanmean(yr_chl_9p(1:18)),18,1),':','color','b','linew',2);
plot(19:27,repmat(nanmean(yr_chl_9p(19:27)),length(19:27),1),':','color','b','linew',2);
plot(28:31,repmat(nanmean(yr_chl_9p(28:31)),length(28:31),1),':','color','b','linew',2);
ylabel('Chl (ug/L)');
ylim([0 7]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)
title(['KOEM yearly OBS NH4']);


plt_kd = yr_no3_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. NO3-N']);
        xlabel('time','fontsize',13)
        ylabel('NO3-N (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl_9p,'b*','linew',2);
plot(yr_chl_9p,'b-','linew',2);
plot(1:18,repmat(nanmean(yr_chl_9p(1:18)),18,1),':','color','b','linew',2);
plot(19:27,repmat(nanmean(yr_chl_9p(19:27)),length(19:27),1),':','color','b','linew',2);
plot(28:31,repmat(nanmean(yr_chl_9p(28:31)),length(28:31),1),':','color','b','linew',2);
ylabel('Chl (ug/L)');
ylim([0 7]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)
title(['KOEM yearly OBS NO3']);

plt_kd = yr_din_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. DIN']);
        xlabel('time','fontsize',13)
        ylabel('DIN (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl_9p,'b*','linew',2);
plot(yr_chl_9p,'b-','linew',2);
plot(1:18,repmat(nanmean(yr_chl_9p(1:18)),18,1),':','color','b','linew',2);
plot(19:27,repmat(nanmean(yr_chl_9p(19:27)),length(19:27),1),':','color','b','linew',2);
plot(28:31,repmat(nanmean(yr_chl_9p(28:31)),length(28:31),1),':','color','b','linew',2);
ylabel('Chl (ug/L)');
ylim([0 7]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)
title(['KOEM yearly OBS DIN']);

plt_kd = yr_po4_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. PO4-P']);
        xlabel('time','fontsize',13)
        ylabel('PO4-P (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl_9p,'b*','linew',2);
plot(yr_chl_9p,'b-','linew',2);
plot(1:18,repmat(nanmean(yr_chl_9p(1:18)),18,1),':','color','b','linew',2);
plot(19:27,repmat(nanmean(yr_chl_9p(19:27)),length(19:27),1),':','color','b','linew',2);
plot(28:31,repmat(nanmean(yr_chl_9p(28:31)),length(28:31),1),':','color','b','linew',2);
ylabel('Chl (ug/L)');
ylim([0 7]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)


plt_kd = yr_temp_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. temp']);
        xlabel('time','fontsize',13)
        ylabel('temperature (^oC)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl_9p,'b*','linew',2);
plot(yr_chl_9p,'b-','linew',2);
plot(1:18,repmat(nanmean(yr_chl_9p(1:18)),18,1),':','color','b','linew',2);
plot(19:27,repmat(nanmean(yr_chl_9p(19:27)),length(19:27),1),':','color','b','linew',2);
plot(28:31,repmat(nanmean(yr_chl_9p(28:31)),length(28:31),1),':','color','b','linew',2);
ylabel('Chl (ug/L)');
ylim([0 7]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)
title(['KOEM yearly OBS temp']);

plt_kd = yr_ss_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl_9p,'b*','linew',2);
plot(yr_chl_9p,'b-','linew',2);
plot(1:18,repmat(nanmean(yr_chl_9p(1:18)),18,1),':','color','b','linew',2);
plot(19:27,repmat(nanmean(yr_chl_9p(19:27)),length(19:27),1),':','color','b','linew',2);
plot(28:31,repmat(nanmean(yr_chl_9p(28:31)),length(28:31),1),':','color','b','linew',2);
ylabel('Chl (ug/L)');
ylim([0 7]);
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)
title(['KOEM yearly OBS SS']);


plt_kd = yr_chl_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['KOEM yearly OBS Chl']);
        xlabel('time','fontsize',13)
        ylabel('Chl (ug/L)');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 7]);
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)


plt_kd = yr_do_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['KOEM yearly OBS DO']);
        xlabel('time','fontsize',13)
        ylabel('DO (mg/L)');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 7]);
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)


plt_kd = yr_salt_9p;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['KOEM yearly OBS salt']);
        xlabel('time','fontsize',13)
        ylabel('salt (PSU)');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 7]);
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
xtickangle(45)
xlim([19 31])
set(gca,'fontsize',16)

plt_kd = yr_trans;
corrcoef(yr_chl_9p(~isnan(yr_chl_9p)),plt_kd(~isnan(yr_chl_9p)))
fig = figure; hold on;
        plot(plt_kd,'k*','linew',2);
        plot(plt_kd,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_kd(1:18)),18,1),':','color',[.5 .5 .5],'linew',2);
        plot(19:27,repmat(nanmean(plt_kd(19:27)),length(19:27),1),':','color',[.5 .5 .5],'linew',2);
        plot(28:31,repmat(nanmean(plt_kd(28:31)),length(28:31),1),':','color',[.5 .5 .5],'linew',2);
        title(['KOEM yearly OBS songjung transp.']);
        xlabel('time','fontsize',13)
        ylabel('river transport (m^3/s)');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 7]);
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
xtickangle(45)
xlim([19 31])
set(gca,'fontsize',16)




%% Temp vs. CHL
corrcoef(yr_chl(~isnan(yr_chl)),yr_temp(~isnan(yr_chl)))
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_temp,'b*','linew',2);
        plot(yr_temp,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
ylim([0 9]);
xtickangle(45)


corrcoef(yr_chl(~isnan(yr_chl)),yr_temp(~isnan(yr_chl)))
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_temp_9p,'b*','linew',2);
        plot(yr_temp_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl_9p,'g-','linew',2);
ylim([0 9]);
xtickangle(45)

%% SS vs. CHL
corrcoef(yr_chl(~isnan(yr_chl)),yr_ss(~isnan(yr_chl)))
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_ss,'b*','linew',2);
        plot(yr_ss,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS Chl. vs. SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
ylim([0 7]);
xtickangle(45)

%% SI vs CHL
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_si .* yr_trans,'b*','linew',2);
        plot(yr_si .* yr_trans,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. Si']);
        xlabel('time','fontsize',13)
        ylabel('Si (mmol si / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
xtickangle(45)


%% NH4 vs CHL
plt_nh4=sj.yr_nh4 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 
plt_chl=yr_chl
r_coef=corrcoef(plt_chl(~isnan(plt_chl)),plt_nh4(~isnan(plt_chl)))

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sj.yr_nh4 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 ,'ko','linew',2);
        plot(sj.yr_nh4 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 ,'k','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS chl vs. NH4 load']);
        xlabel('time','fontsize',13)
        ylabel('NH4-N load (ton/yr)','fontsize',13)
        text(2,250, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);      
yyaxis right
plot(yr_chl ,'b-','linew',2);
xtickangle(45)
 print(fig,strcat('NH4_loads_vs_chl_sp'),'-dpng') 
 
% nh4 load vs. chl 
clearvars plt_nh4
plt_nh4=sj.yr_nh4 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 
for i = 1:16
clearvars r_coef plt_chl coe_*
    plt_chl = yr_chl_raw(i,:);
    coe_nh4=plt_nh4(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_nh4)),coe_nh4(~isnan(coe_nh4)))    
fig = figure; hold on;
        plot(plt_nh4,'ko','linew',2);
        plot(plt_nh4,'k','linew',2);
        title(['1989-2019 yearly KOEM Chl. vs. nh4 load st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('NH4-N load (ton/yr)','fontsize',13)
        text(2,300, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_nh4_load_koem_st_',char(name_tag{i})]),'-dpng')
end


% sj chl load vs. chl 
clearvars load_chl
load_chl=sj.yr_chl .* yr_trans .* 31536000 ./ 10^6 
for i = 1:16
clearvars r_coef plt_chl coe_*
    plt_chl = yr_chl_raw(i,:);
    coe_load_chl=load_chl(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_load_chl)),coe_load_chl(~isnan(coe_load_chl)))    
fig = figure; hold on;
        plot(load_chl,'ko','linew',2);
        plot(load_chl,'k','linew',2);
        title(['1989-2019 yearly KOEM Chl. vs. load chl load st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('chl load (ton/yr)','fontsize',13)
        text(2,300, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_load_chl_koem_st_',char(name_tag{i})]),'-dpng')
end

%% NO3 vs CHL
plt_no3=sj.yr_no3 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 
plt_chl=yr_chl
r_coef=corrcoef(plt_chl(~isnan(plt_chl)),plt_no3(~isnan(plt_chl)))

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sj.yr_no3 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 ,'ko','linew',2);
        plot(sj.yr_no3 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 ,'k','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS chl vs. NO3 load']);
        xlabel('time','fontsize',13)
        ylabel('NO3-N load (ton/yr)','fontsize',13)
        text(2,2000, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(yr_chl ,'b-','linew',2);
xtickangle(45)
print(fig,strcat('NO3_loads_vs_chl_sp'),'-dpng') 

% no3 load vs. chl 
clearvars plt_no3
plt_no3=sj.yr_no3 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 
for i = 1:16
clearvars r_coef plt_chl coe_*
    plt_chl = yr_chl_raw(i,:);
    coe_no3=plt_no3(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_no3)),coe_no3(~isnan(coe_no3)))    
fig = figure; hold on;
        plot(plt_no3,'ko','linew',2);
        plot(plt_no3,'k','linew',2);
        title(['1989-2019 yearly KOEM Chl. vs. no3 load st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('NO3-N load (ton/yr)','fontsize',13)
        text(3,2000, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_no3_load_koem_st_',char(name_tag{i})]),'-dpng')
end

% load nut. 
clearvars plt_nh4 plt_no3 plt_po4
plt_nh4=sj.yr_nh4 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 
plt_no3=sj.yr_no3 ./1000 .*14.006720.* yr_trans .* 31536000 ./ 10^6 
plt_po4=sj.yr_po4 ./1000 .*30.973762.* yr_trans .* 31536000 ./ 10^6
fig = figure; hold on;
        plot(plt_no3,'ro','linew',2);
        plot(plt_nh4,'ko','linew',2);
         plot(plt_no3,'r','linew',2);
        plot(plt_nh4,'k','linew',2);
        plot(1:18,repmat(nanmean(plt_no3(1:18)),18,1),':','color','r','linew',2);
        plot(19:27,repmat(nanmean(plt_no3(19:27)),length(19:27),1),':','color','r','linew',2);
        plot(28:31,repmat(nanmean(plt_no3(28:31)),length(28:31),1),':','color','r','linew',2);
        plot(1:18,repmat(nanmean(plt_nh4(1:18)),18,1),':','color','k','linew',2);
        plot(19:27,repmat(nanmean(plt_nh4(19:27)),length(19:27),1),':','color','k','linew',2);
        plot(28:31,repmat(nanmean(plt_nh4(28:31)),length(28:31),1),':','color','k','linew',2);
        title(['NO3-N, NH4-N, PO4-P load songjung']);
        xlabel('time','fontsize',13)
        ylabel('load (ton/yr)','fontsize',13)
        legend('NO3-N','NH4-N');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([7 31])
        ylim([0 inf])
%         ylim([0 5000])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_po4,'bo','linew',2);
plot(plt_po4,'b','linew',2);
plot(13:18,repmat(nanmean(plt_po4(13:18)),length(13:18),1),':','color','b','linew',2);
plot(19:27,repmat(nanmean(plt_po4(19:27)),length(19:27),1),':','color','b','linew',2);
plot(28:31,repmat(nanmean(plt_po4(28:31)),length(28:31),1),':','color','b','linew',2);
ylabel('PO4-P load (ton/yr)','fontsize',13)
ylim([0 inf])
xtickangle(45)
print(fig,strcat(['load_nut_river_yearly']),'-dpng')
ax = gca;
% ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
xtickangle(45)
% 2&3reg
xlim([19 31])
set(gca,'fontsize',16)



clearvars plt_nh4 plt_no3 plt_po4 plt_ss
plt_ss=sj.yr_ss .* yr_trans.* 31536000 ./ 10^6
fig = figure; hold on;
        plot(plt_ss,'ko','linew',2);
        plot(plt_ss,'k','linew',2);
        title(['SS load songjung']);
        xlabel('time','fontsize',13)
        ylabel('load (ton/yr)','fontsize',13)
%         legend('NO3-N','NH4-N');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
xtickangle(45)
print(fig,strcat(['load_ss_river_yearly']),'-dpng')


clearvars plt_chl
plt_chl=sj.yr_chl .* yr_trans.* 31536000 ./ 10^9  % chl is ug/L (= mg/m^3)
fig = figure; hold on;
        plot(plt_chl,'ko','linew',2);
        plot(plt_chl,'k','linew',2);
        title(['Chl load songjung']);
        xlabel('time','fontsize',13)
        ylabel('load (ton/yr)','fontsize',13)
%         legend('NO3-N','NH4-N');
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
xtickangle(45)
print(fig,strcat(['load_chl_river_yearly']),'-dpng')

%% PO4 vs CHL
plt_po4=sj.yr_po4 ./1000 .*30.973762.* yr_trans .* 31536000 ./ 10^6 
plt_chl=yr_chl
r_coef=corrcoef(plt_chl(~isnan(plt_po4)),plt_po4(~isnan(plt_po4)))

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sj.yr_po4 ./1000 .*30.973762.* yr_trans .* 31536000 ./ 10^6,'ko','linew',2);
        plot(sj.yr_po4 ./1000 .*30.973762.* yr_trans .* 31536000 ./ 10^6,'k','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS chl vs. PO4 load']);
        xlabel('time','fontsize',13)
        ylabel('PO4-P load (ton/yr)','fontsize',13)
        text(3,30, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(yr_chl,'b-','linew',2);
xtickangle(45)
print(fig,strcat('po4_loads_vs_chl_sp'),'-dpng') 

isplt=yr_chl.* yr_trans
ispltnan=find(isnan(isplt)==1);
plt=NaN(1,30); plt(13:end)=1;
plot(yr_po4 .* yr_trans , yr_chl.* yr_trans,'ro');

% PO4 vs CHL each
clearvars plt_po4
plt_po4=sj.yr_po4 ./1000 .*30.973762.* yr_trans .* 31536000 ./ 10^6 ; % ton/yr
for i = 1:16
clearvars  r_coef plt_chl coe_*
    plt_chl = yr_chl_raw(i,:);
    coe_po4=plt_po4(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_po4)),coe_po4(~isnan(coe_po4)))    
fig = figure; hold on;
        plot(plt_po4,'ko','linew',2);
        plot(plt_po4,'k','linew',2);
        title(['1989-2019 yearly KOEM Chl. vs. PO4 load st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('PO4-P load (ton/yr)','fontsize',13)
        text(5,30, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 130])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_PO4_load_koem_st_',char(name_tag{i})]),'-dpng')
end


%% ss vs CHL
plt_ss=sj.yr_ss .* yr_trans.* 31536000 ./ 10^6
plt_chl=yr_chl
r_coef=corrcoef(plt_chl(~isnan(plt_chl)),plt_ss(~isnan(plt_chl)))

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sj.yr_ss .* yr_trans .* 31536000 ./ 10^6,'ko','linew',2);
        plot(sj.yr_ss .* yr_trans .* 31536000 ./ 10^6,'k','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS chl vs. SS load']);
        xlabel('time','fontsize',13)
        ylabel('SS load (ton/yr)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl ,'b-','linew',2);
% ylim([2 9])
xtickangle(45)
print(fig,strcat('ss_loads_vs_chl_sp'),'-dpng') 

% ss vs chl each st
clearvars plt_ss
plt_ss=sj.yr_ss .* yr_trans.* 31536000 ./ 10^6
for i = 1:16
clearvars r_coef plt_chl coe_*
    plt_chl = yr_chl_raw(i,:);
    coe_ss=plt_ss(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_ss)),coe_ss(~isnan(coe_ss)))    
fig = figure; hold on;
        plot(plt_ss,'ko','linew',2);
        plot(plt_ss,'k','linew',2);
        title(['1989-2019 yearly KOEM Chl. vs. SS load st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('ss (ton/yr)','fontsize',13)
        text(5,25000, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_ss_load_koem_st_',char(name_tag{i})]),'-dpng')
end


%% SS
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_ss,'b*','linew',2);
        plot(yr_ss,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)

%% Si
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_si,'b*','linew',2);
        plot(yr_si,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. Si']);
        xlabel('time','fontsize',13)
        ylabel('Si (mmol si / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(yr_nh4./yr_si,'b*','linew',2);
%         plot(yr_nh4./yr_si,'r','linew',2);
%         plot(yr_po4./yr_si,'b*','linew',2);
%         plot(yr_po4./yr_si,'r','linew',2);
%         plot((yr_no3 + yr_nh4)./yr_si,'b*','linew',2);
%         plot((yr_no3+yr_nh4)./yr_si,'r','linew',2);
        
        plot((yr_no3 + yr_nh4)./yr_po4,'b*','linew',2);
        plot((yr_no3 + yr_nh4)./yr_po4,'r','linew',2);
        
        title(['GY 1989-2019 yearly OBS transp vs. Si']);
        xlabel('time','fontsize',13)
        ylabel('Si (mmol si / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(yr_nh4./yr_si,'b*','linew',2);
%         plot(yr_nh4./yr_si,'r','linew',2);
%         plot(yr_po4./yr_si,'b*','linew',2);
%         plot(yr_po4./yr_si,'r','linew',2);
%         plot((yr_no3 + yr_nh4)./yr_si,'b*','linew',2);
%         plot((yr_no3+yr_nh4)./yr_si,'r','linew',2);    
%         plot((sj.yr_no3 + sj.yr_nh4)./sj.yr_po4,'b*','linew',2);
%         plot((sj.yr_no3 + sj.yr_nh4)./sj.yr_po4,'r','linew',2);        
%         plot((sj.yr_no3)./sj.yr_po4,'b*','linew',2);
%         plot((sj.yr_no3)./sj.yr_po4,'r','linew',2);
%         plot((sj.yr_nh4)./sj.yr_po4,'b*','linew',2);
%         plot((sj.yr_nh4)./sj.yr_po4,'r','linew',2);
        plot((sj.yr_no3 + sj.yr_nh4),'b*','linew',2);
        plot((sj.yr_no3 + sj.yr_nh4),'r','linew',2);   
        
        title(['GY 1989-2019 yearly OBS transp vs. Si']);
        xlabel('time','fontsize',13)
        ylabel('Si (mmol si / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_chl,'g-','linew',2);
xtickangle(45)

%% secchi
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_secchi,'b*','linew',2);
        plot(yr_secchi,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. secchi']);
        xlabel('time','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(1.7 ./yr_secchi,'b*','linew',2);
        plot(1.7 ./yr_secchi,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. Kd']);
        xlabel('time','fontsize',13)
        ylabel('Kd (m^-1)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)

%% chl
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_chl,'b*','linew',2);
        plot(yr_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end

%% NH4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_nh4,'b*','linew',2);
        plot(yr_nh4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. nh4']);
        xlabel('time','fontsize',13)
        ylabel('NH4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_no3,'b*','linew',2);
        plot(yr_no3,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. no3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_po4,'b*','linew',2);
        plot(yr_po4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. PO4']);
        xlabel('time','fontsize',13)
        ylabel('PO4 (mmol P / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)

%% bot

%% Si
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_si_b,'b*','linew',2);
        plot(yr_si_b,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. Si bot']);
        xlabel('time','fontsize',13)
        ylabel('Si (mmol Si / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)

%% SS
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_ss_b,'b*','linew',2);
        plot(yr_ss_b,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. SS bot']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end

%% CHL
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_chl_b,'b*','linew',2);
        plot(yr_chl_b,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. chl bot']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end

%% NH4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_nh4_b,'b*','linew',2);
        plot(yr_nh4_b,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. nh4']);
        xlabel('time','fontsize',13)
        ylabel('NH4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_no3_b,'b*','linew',2);
        plot(yr_no3_b,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. no3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_po4,'b*','linew',2);
        plot(yr_po4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. PO4']);
        xlabel('time','fontsize',13)
        ylabel('PO4 (mmol P / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)

%% summer
%% chl
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(summer_chl,'b*','linew',2);
        plot(summer_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 summer OBS transp vs. chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(summer_trans,'g-','linew',2);
xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end

%% NH4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(summer_nh4,'b*','linew',2);
        plot(summer_nh4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 summer OBS transp vs. nh4']);
        xlabel('time','fontsize',13)
        ylabel('NH4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(summer_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(summer_no3,'b*','linew',2);
        plot(summer_no3,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 summer OBS transp vs. no3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(summer_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(summer_po4,'b*','linew',2);
        plot(summer_po4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 summer OBS transp vs. PO4']);
        xlabel('time','fontsize',13)
        ylabel('PO4 (mmol P / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(summer_trans,'g-','linew',2); xtickangle(45)

%% bot
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(summer_chl_b,'b*','linew',2);
        plot(summer_chl_b,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 summer OBS transp vs. chl bot']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(summer_trans,'g-','linew',2);
xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end

%% NH4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(summer_nh4_b,'b*','linew',2);
        plot(summer_nh4_b,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 summer OBS transp vs. nh4']);
        xlabel('time','fontsize',13)
        ylabel('NH4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(summer_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(summer_no3_b,'b*','linew',2);
        plot(summer_no3_b,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 summer OBS transp vs. no3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(summer_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(summer_po4,'b*','linew',2);
        plot(summer_po4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 summer OBS transp vs. PO4']);
        xlabel('time','fontsize',13)
        ylabel('PO4 (mmol P / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(summer_trans,'g-','linew',2); xtickangle(45)


%% MONTHLY
%% chl
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(koem_chl,'b*','linew',2);
        plot(koem_in_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 monthly OBS transp vs. chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(koem_chl)])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(monthly_trans,'g-','linew',2);
% for i = 6:12:360; xline(i,'color','m'); end

%% NH4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(koem_nh4,'b*','linew',2);
        plot(koem_in_nh4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 monthly OBS transp vs. nh4']);
        xlabel('time','fontsize',13)
        ylabel('NH4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(koem_chl)])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(monthly_trans,'g-','linew',2);
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(koem_no3,'b*','linew',2);
        plot(koem_in_no3,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 monthly OBS transp vs. no3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(koem_chl)])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(monthly_trans,'g-','linew',2);
for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(koem_po4,'b*','linew',2);
        plot(koem_in_po4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 monthly OBS transp vs. po4']);
        xlabel('time','fontsize',13)
        ylabel('po4 (mmol P / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(koem_chl)])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(monthly_trans,'g-','linew',2);
for i = 6:12:360; xline(i,'color','m'); end

return

corrcoef(monthly_trans(~isnan(koem_in_chl)),koem_in_chl(~isnan(koem_in_chl)))
% corrcoef(sj_trans_out,interp_chl,'Rows','complete')
corrcoef(monthly_trans(~isnan(koem_chl)),koem_chl(~isnan(koem_chl)))

figure;
plot(monthly_trans(~isnan(koem_in_chl)),koem_in_chl(~isnan(koem_in_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('Chla (ug/L)','fontsize',13)
grid on;
title(['GY 1989-2019 monthly interped chl vs. transport']);

figure;
plot(monthly_trans(~isnan(koem_chl)),koem_chl(~isnan(koem_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('Chla (ug/L)','fontsize',13)
grid on;
title(['GY 1989-2019 monthly chl vs. transport']);


clearvars r lags
[r,lags] = xcorr(sj_trans_out_x(~isnan(interp_chl)),interp_chl_x(~isnan(interp_chl)),'normalized') 
% [r,lags] = xcorr(sj_trans_out_x(~isnan(raw_chl)),raw_chl_x(~isnan(raw_chl)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily chl vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL chl vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

sj_tran_in = sj_trans_out(~isnan(raw_chl));
raw_chl_in = raw_chl(~isnan(raw_chl));
clearvars r lags
[r,lags] = xcorr(sj_tran_in - nanmean(sj_tran_in),raw_chl_in - nanmean(raw_chl_in), 'normalized') 
% [r,lags] = xcorr(sj_trans_out(~isnan(interp_chl)),interp_chl(~isnan(interp_chl)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily chl vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');



%% nh4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(raw_nh4,'b*','linew',2);
        plot(interp_nh4,'r','linew',2);
%         plot(1:length(obm_gy_nh4),obm_gy_nh4 + obs_std_gy_nh4,'m-','linew',2);
%         plot(1:length(obm_gy_nh4),obm_gy_nh4 - obs_std_gy_nh4,'m-','linew',2);
        title(['GY 2001 daily OBS transp vs. nh4']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(raw_chl)])
        set(gca,'xticklabel',2007:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_nh4'),'-dpng') 
yyaxis right
plot(sj_trans_out,'g-','linew',2);
xlim([1 length(raw_nh4)])

return

corrcoef(sj_trans_out(~isnan(interp_nh4)),interp_nh4(~isnan(interp_nh4)))
% corrcoef(sj_trans_out,interp_nh4,'Rows','complete')
corrcoef(sj_trans_out(~isnan(raw_nh4)),raw_nh4(~isnan(raw_nh4)))

figure;
plot(sj_trans_out(~isnan(interp_nh4)),interp_nh4(~isnan(interp_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['GY 2001 daily interped nh4 vs. transport']);

figure;
plot(sj_trans_out(~isnan(raw_nh4)),raw_nh4(~isnan(raw_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['GY 2001 daily nh4 vs. transport']);



clearvars r lags
% [r,lags] = xcorr(sj_trans_out_x(~isnan(interp_nh4)),interp_nh4_x(~isnan(interp_nh4)),'normalized') 
[r,lags] = xcorr(sj_trans_out_x(~isnan(raw_nh4)),raw_nh4_x(~isnan(raw_nh4)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily nh4 vs. transport']);
n_nan_nh4 = find(isnan(raw_nh4)==0)



text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL nh4 vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');


%% no3
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(raw_no3,'b*','linew',2);
        plot(interp_no3,'r','linew',2);
%         plot(1:length(obm_gy_no3),obm_gy_no3 + obs_std_gy_no3,'m-','linew',2);
%         plot(1:length(obm_gy_no3),obm_gy_no3 - obs_std_gy_no3,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL no3 + transport']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(raw_chl)])
        set(gca,'xticklabel',2007:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3'),'-dpng') 
yyaxis right
plot(sj_trans_out,'g-','linew',2);
xlim([1 length(raw_no3)])

return

corrcoef(sj_trans_out(~isnan(interp_no3)),interp_no3(~isnan(interp_no3)))
% corrcoef(sj_trans_out,interp_no3,'Rows','complete')
corrcoef(sj_trans_out(~isnan(raw_no3)),raw_no3(~isnan(raw_no3)))

figure;
plot(sj_trans_out(~isnan(interp_no3)),interp_no3(~isnan(interp_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['GY 2001 daily interped no3 vs. transport']);

figure;
plot(sj_trans_out(~isnan(raw_no3)),raw_no3(~isnan(raw_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['GY 2001 daily no3 vs. transport']);



clearvars r lags
% [r,lags] = xcorr(sj_trans_out_x(~isnan(interp_no3)),interp_no3_x(~isnan(interp_no3)),'normalized') 
[r,lags] = xcorr(sj_trans_out_x(~isnan(raw_no3)),raw_no3_x(~isnan(raw_no3)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily no3 vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL no3 vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');



%% do
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(raw_do,'b*','linew',2);
        plot(interp_do,'r','linew',2);
%         plot(1:length(obm_gy_do),obm_gy_do + obs_std_gy_do,'m-','linew',2);
%         plot(1:length(obm_gy_do),obm_gy_do - obs_std_gy_do,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL do + transport']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('do (mmol / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[t_tick(1:end)]);
        xlim([1 length(raw_chl)])
        set(gca,'xticklabel',2007:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_do'),'-dpng') 
yyaxis right
plot(sj_trans_out,'g-','linew',2);
xlim([1 length(raw_do)])

return

corrcoef(sj_trans_out(~isnan(interp_do)),interp_do(~isnan(interp_do)))
% corrcoef(sj_trans_out,interp_do,'Rows','complete')
corrcoef(sj_trans_out(~isnan(raw_do)),raw_do(~isnan(raw_do)))

figure;
plot(sj_trans_out(~isnan(interp_do)),interp_do(~isnan(interp_do)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('do (mmol/ m^3)','fontsize',13)
grid on;
title(['GY 2001 daily interped do vs. transport']);

figure;
plot(sj_trans_out(~isnan(raw_do)),raw_do(~isnan(raw_do)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('do (mmol/ m^3)','fontsize',13)
grid on;
title(['GY 2001 daily do vs. transport']);


clearvars r lags
% [r,lags] = xcorr(sj_trans_out_x(~isnan(interp_do)),interp_do_x(~isnan(interp_do)),'normalized') 
[r,lags] = xcorr(sj_trans_out_x(~isnan(raw_do)),raw_do_x(~isnan(raw_do)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily do vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL do vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');



%% po4
%% chl

corrcoef(monthly_trans(~isnan(koem_in_chl)),koem_in_chl(~isnan(koem_in_chl)))
% corrcoef(monthly_trans,koem_in_chl,'Rows','complete')
corrcoef(monthly_trans(~isnan(koem_chl)),koem_chl(~isnan(koem_chl)))
corrcoef(yr_trans(~isnan(yr_chl)),yr_chl(~isnan(yr_chl)))
corrcoef(summer_trans(~isnan(summer_chl)),summer_chl(~isnan(summer_chl)))

figure;
plot(monthly_trans(~isnan(koem_in_chl)),koem_in_chl(~isnan(koem_in_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
grid on;
title(['monthly interped chl vs. transport']);

figure;
plot(monthly_trans(~isnan(koem_chl)),koem_chl(~isnan(koem_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
grid on;
title(['monthly chl vs. transport']);

figure;
plot(yr_trans(~isnan(yr_chl)),yr_chl(~isnan(yr_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
grid on;
title(['yearly chl vs. transport']);
ylim([0 11]);  xlim([0 140])

figure;
plot(summer_trans(~isnan(summer_chl)),summer_chl(~isnan(summer_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
grid on;
title(['yearly chl vs. transport']);
ylim([0 11]);  xlim([0 140])

%% nh4

corrcoef(monthly_trans(~isnan(koem_in_nh4)),koem_in_nh4(~isnan(koem_in_nh4)))
% corrcoef(monthly_trans,koem_in_nh4,'Rows','complete')
corrcoef(monthly_trans(~isnan(koem_nh4)),koem_nh4(~isnan(koem_nh4)))
corrcoef(yr_trans(~isnan(yr_nh4)),yr_nh4(~isnan(yr_nh4)))
corrcoef(summer_trans(~isnan(summer_nh4)),summer_nh4(~isnan(summer_nh4)))

figure;
plot(monthly_trans(~isnan(koem_in_nh4)),koem_in_nh4(~isnan(koem_in_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly interped nh4 vs. transport']);

figure;
plot(monthly_trans(~isnan(koem_nh4)),koem_nh4(~isnan(koem_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly nh4 vs. transport']);


figure;
plot(yr_trans(~isnan(yr_nh4)),yr_nh4(~isnan(yr_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly nh4 vs. transport']);
xlim([0 140])


figure;
plot(summer_trans(~isnan(summer_nh4)),summer_nh4(~isnan(summer_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly nh4 vs. transport']);
xlim([0 140])


%% no3

corrcoef(monthly_trans(~isnan(koem_in_no3)),koem_in_no3(~isnan(koem_in_no3)))
% corrcoef(monthly_trans,koem_in_no3,'Rows','complete')
corrcoef(monthly_trans(~isnan(koem_no3)),koem_no3(~isnan(koem_no3)))
corrcoef(yr_trans(~isnan(yr_no3)),yr_no3(~isnan(yr_no3)))
corrcoef(summer_trans(~isnan(summer_no3)),summer_no3(~isnan(summer_no3)))
figure;
plot(monthly_trans(~isnan(koem_in_no3)),koem_in_no3(~isnan(koem_in_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly interped no3 vs. transport']);

figure;
plot(monthly_trans(~isnan(koem_no3)),koem_no3(~isnan(koem_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly no3 vs. transport']);


figure;
plot(yr_trans(~isnan(yr_no3)),yr_no3(~isnan(yr_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly no3 vs. transport']);
xlim([0 140])


figure;
plot(summer_trans(~isnan(summer_no3)),summer_no3(~isnan(summer_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly no3 vs. transport']);
xlim([0 140])


%% po4

corrcoef(monthly_trans(~isnan(koem_in_po4)),koem_in_po4(~isnan(koem_in_po4)))
% corrcoef(monthly_trans,koem_in_po4,'Rows','complete')
corrcoef(monthly_trans(~isnan(koem_po4)),koem_po4(~isnan(koem_po4)))
corrcoef(yr_trans(~isnan(yr_po4)),yr_po4(~isnan(yr_po4)))
corrcoef(summer_trans(~isnan(summer_po4)),summer_po4(~isnan(summer_po4)))

figure;
plot(monthly_trans(~isnan(koem_in_po4)),koem_in_po4(~isnan(koem_in_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly interped po4 vs. transport']);

figure;
plot(monthly_trans(~isnan(koem_po4)),koem_po4(~isnan(koem_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly po4 vs. transport']);


figure;
plot(yr_trans(~isnan(yr_po4)),yr_po4(~isnan(yr_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly po4 vs. transport']);
xlim([0 140])

figure;
plot(summer_trans(~isnan(summer_po4)),summer_po4(~isnan(summer_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly po4 vs. transport']);
xlim([0 140])


corrcoef(monthly_trans(~isnan(koem_in_po4)),koem_in_po4(~isnan(koem_in_po4)))
% corrcoef(monthly_trans,koem_in_po4,'Rows','complete')
corrcoef(monthly_trans(~isnan(koem_po4)),koem_po4(~isnan(koem_po4)))
corrcoef(yr_trans(~isnan(yr_po4)),yr_po4(~isnan(yr_po4)))

figure;
plot(monthly_trans(~isnan(koem_in_po4)),koem_in_po4(~isnan(koem_in_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly interped po4 vs. transport']);

figure;
plot(monthly_trans(~isnan(koem_po4)),koem_po4(~isnan(koem_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['monthly po4 vs. transport']);


figure;
plot(yr_trans(~isnan(yr_po4)),yr_po4(~isnan(yr_po4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('po4 (mmol N / m^3)','fontsize',13)
grid on;
title(['yearly po4 vs. transport']);



clearvars r lags
% [r,lags] = xcorr(monthly_trans_x(~isnan(koem_in_po4)),koem_in_po4_x(~isnan(koem_in_po4)),'normalized') 
[r,lags] = xcorr(sj_trans_out_x(~isnan(koem_po4)),koem_po4_x(~isnan(koem_po4)),'coeff')
figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily po4 vs. transport']);
n_nan_po4 = find(isnan(koem_po4)==0)



text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');

figure;hold on;
plot(lags,r); 
xline(lags(find(max(r) == r)),'color','r')
xlabel('time lag (days)','fontsize',13)
ylabel('corrcoef','fontsize',13)
grid on;
title(['GY 2001 daily MODEL po4 vs. transport']);
text(-360,0.7,['corrcoef. = ',num2str(max(r),'%.1f'),' ( lag = ', num2str(lags(find(max(r) == r))), 'day )'],'color','r');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% June concentration vs. transp. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHL
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(june_chl,'b*','linew',2);
        plot(june_in_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 june OBS transp vs. chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[1:length(1989:2018)]);
        xlim([1 length(1989:2018)])
        set(gca,'xticklabel',1989:2018,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(june_trans,'g-','linew',2);
% for i = 6:12:360; xline(i,'color','m'); end

return

corrcoef(june_trans(~isnan(june_in_chl)),june_in_chl(~isnan(june_in_chl)))
% corrcoef(monthly_trans,interp_chl,'Rows','complete')
corrcoef(june_trans(~isnan(june_chl)),june_chl(~isnan(june_chl)))

figure;
plot(june_trans(~isnan(june_in_chl)),june_in_chl(~isnan(june_in_chl)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('Chla (ug/L)','fontsize',13)
grid on;
title(['GY 1989-2019 june interped chl vs. transport']);

% NH4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(june_nh4,'b*','linew',2);
        plot(june_in_nh4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 june OBS transp vs. nh4']);
        xlabel('time','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[1:length(1989:2018)]);
        xlim([1 length(1989:2018)])
        set(gca,'xticklabel',1989:2018,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(june_trans,'g-','linew',2);

corrcoef(june_trans(~isnan(june_in_nh4)),june_in_nh4(~isnan(june_in_nh4)))
% corrcoef(sj_trans_out,interp_nh4,'Rows','complete')
corrcoef(june_trans(~isnan(june_nh4)),june_nh4(~isnan(june_nh4)))

figure;
plot(june_trans(~isnan(june_in_nh4)),june_in_nh4(~isnan(june_in_nh4)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('nh4 (mmol N /m^3)','fontsize',13)
grid on;
title(['GY 1989-2019 june interped nh4 vs. transport']);

% NO3
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(june_no3,'b*','linew',2);
        plot(june_in_no3,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 june OBS transp vs. no3']);
        xlabel('time','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[1:length(1989:2018)]);
        xlim([1 length(1989:2018)])
        set(gca,'xticklabel',1989:2018,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(june_trans,'g-','linew',2);

corrcoef(june_trans(~isnan(june_in_no3)),june_in_no3(~isnan(june_in_no3)))
% corrcoef(sj_trans_out,interp_no3,'Rows','complete')
corrcoef(june_trans(~isnan(june_no3)),june_no3(~isnan(june_no3)))

figure;
plot(june_trans(~isnan(june_in_no3)),june_in_no3(~isnan(june_in_no3)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('no3 (mmol N /m^3)','fontsize',13)
grid on;
title(['GY 1989-2019 june interped no3 vs. transport']);

% DO
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(june_do,'b*','linew',2);
        plot(june_in_do,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 june OBS transp vs. do']);
        xlabel('time','fontsize',13)
        ylabel('do (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[1:length(1989:2018)]);
        xlim([1 length(1989:2018)])
        set(gca,'xticklabel',1989:2018,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(june_trans,'g-','linew',2);

corrcoef(june_trans(~isnan(june_in_do)),june_in_do(~isnan(june_in_do)))
% corrcoef(sj_trans_out,interp_do,'Rows','complete')
corrcoef(june_trans(~isnan(june_do)),june_do(~isnan(june_do)))

figure;
plot(june_trans(~isnan(june_in_do)),june_in_do(~isnan(june_in_do)),'*')
xlabel('transport (m^3/s)','fontsize',13)
ylabel('do (mmol/m^3)','fontsize',13)
grid on;
title(['GY 1989-2019 june interped do vs. transport']);




%% total mass through the river

figure;
plot((monthly_trans) .* koem_in_chl)  % m^3/s * mg/m^3 = mg/s
hold; 
yyaxis right
plot(monthly_trans,'r');


figure;
plot((monthly_trans) .* (koem_in_nh4 ./1000 .* 14 ))  % m^3/s * g/m^3 = g/s
hold; 
yyaxis right
plot(monthly_trans,'r');


figure;
plot((monthly_trans) .* (koem_in_no3 ./1000 .* 14 ))  % m^3/s * g/m^3 = g/s
hold; 
yyaxis right
plot(monthly_trans,'r');

% figure;
% plot((monthly_trans) .* (koem_in_no3 ./1000 .* 14 ))  % m^3/s * g/m^3 = g/s
% hold; 
% yyaxis right
% plot(monthly_trans,'r');

cd D:\장기생태\Dynamic\06_river
load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','yp_w_*_04');   % polynomial_fitting_on_ecological_variable_to2004_new_v5_3sigma.m 
% ~2003 : advanced 전체기간
% 2004~ : 04~ data's climate
cd D:\장기생태\Dynamic\06_river\환경과학원
load('sumjin(songjung)_polynomial_climate_to2004_3sig.mat','yp_w_*_af'); % ~2004 : 1st regime, 2004~2018 : 2nd regime
% sumjin_do = yp_w_do_04'*0.7*44.661;
% sumjin_chl = yp_w_chl_04';
% sumjin_no3 = yp_w_no3_04'*1000/14;
% sumjin_nh4 = yp_w_nh4_04'*1000/14;

clearvars all_recon_*
for i = 1:30
   yy = t_year(i);
   clearvars tem_*
   if i >= 15  %1989~2003
   tem_do = yp_w_do_04;
   tem_chl = yp_w_chl_04;
   tem_no3 = yp_w_no3_04;
   tem_nh4 = yp_w_nh4_04;
   else %2004~
   tem_do = yp_w_do_af;
   tem_chl = yp_w_chl_af;
   tem_no3 = yp_w_no3_af;
   tem_nh4 = yp_w_nh4_af; 
   end
if leapyear(yy) == 0
   tem_do(60) = [];
   tem_chl(60) = [];
   tem_no3(60) = [];
   tem_nh4(60) = [];
end
if i == 1
    all_recon_do = tem_do;
    all_recon_chl = tem_chl;
    all_recon_no3 = tem_no3;
    all_recon_nh4 = tem_nh4;
else
   all_recon_do=cat(2,all_recon_do,tem_do);
   all_recon_chl=cat(2,all_recon_chl,tem_chl);
   all_recon_no3=cat(2,all_recon_no3,tem_no3);
   all_recon_nh4=cat(2,all_recon_nh4,tem_nh4);
end

end

% daily transport of nutrients too many data for plot it
figure;hold; 
plot((sj_trans_out) .* all_recon_chl,'b')  % m^3/s * mg/m^3 = mg/s
hold on;
yyaxis right
plot(sj_trans_out,'r');
hold off

figure;hold; 
plot(((sj_trans_out) .* (all_recon_nh4))./1000,'b')  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on;
yyaxis right
plot(sj_trans_out,'r');
hold off


figure;hold; 
plot(((sj_trans_out) .* (all_recon_no3))./1000,'b');  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
hold on;
% g/s /1000 = kg/s
yyaxis right
plot(sj_trans_out,'r'); hold off

% figure;
% plot((sj_trans_out) .* (koem_in_no3 ./1000 .* 14 ))  % m^3/s * mg/m^3 = mg/s
% hold; 
% yyaxis right
% plot(sj_trans_out,'r');

% save('transp_daily_and_nutrients.mat','all_recon_*','sj_trans_out');

% daily on 2001 only

figure;hold; 
plot((dis_pre_total{2001-1979}) .* yp_w_chl_04,'color','b','linew',2)  % m^3/s * mg/m^3 = mg/s
ylabel('Chla x transport (mg/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
hold off

figure;hold; 
plot((dis_pre_total{2001-1979}) .* yp_w_nh4_04,'color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
ylabel('nh4 x transport (g/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
hold off


figure;hold; 
plot((dis_pre_total{2001-1979}) .* yp_w_no3_04,'color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
ylabel('no3 x transport (g/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
hold off


tt_tick = [1 t_indx(1:11)+1];
figure;hold; 
for i = 8:30
plot(double(dis_pre_total{t_year(i)-1979}),'b'); hold on;  % m^3/s * mg/m^3 = mg/s
end
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title('1996~2018 songjung transport')
% set(gca,'xtick',[tt_tick(1:end)]);
% set(gca,'xticklabel',,'fontsize',10);


%% find when 6 month's transp. is larger than summer

for i = 8:30
figure; 
plot(double(dis_pre_total{t_year(i)-1979}),'b');  % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_year(i)) ' songjung transport'])
saveas(gcf,[num2str(t_year(i)) '-songjung-transport.png']);
close all;
end

for i = 8:30
figure; 
plot(double(dis_pre_total{t_year(i)-1979}),'b');  % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_year(i)) ' songjung transport'])
ylim([0 5000])
saveas(gcf,[num2str(t_year(i)) '-songjung-transport.png']);
close all;
end

cd D:\장기생태\Dynamic\06_river\환경과학원
clearvars raw_*
load('sumjin(jinwall)_polynomial_climate_advanced(v4)_3sig.mat','raw_*'); 

r_do_jinwal = raw_do *0.7*44.661;
r_chl_jinwal = raw_chl;
r_no3_jinwal = raw_no3 .* 1000 ./14;
r_nh4_jinwal = raw_nh4 .* 1000 ./14;

t_07=load('2007_2019_time_axis.mat','t_tick');
t_07_tick = t_07.t_tick;

cd D:\장기생태\Dynamic\06_river
clearvars raw_*
load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat'); 
r_do_song = raw_do *0.7*44.661;
r_chl_song = raw_chl;
r_no3_song = raw_no3 .* 1000 ./14;
r_nh4_song = raw_nh4 .* 1000 ./14;

t_indx_year = t_indx(12:12:end); %1989~2019

% when 6 month's transp. larger
% 1996 : 8th year
% 1998 : 10th year
% 2001 : 13th year
% 2008 : 20th year 
nut_big_y_indx = [8, 10, 13, 20];
t_big_year = [1996, 1998, 2001, 2008];
for i = 1:4
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS nh4'])
ylim([0 5000])
yyaxis right
plot(r_nh4_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('NH4 (mmol N /m^3)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-nh4.png']);
close all;
end

for i = 1:4
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS no3'])
ylim([0 5000])
yyaxis right
plot(r_no3_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('NO3 (mmol N /m^3)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-no3.png']);
close all;
end

for i = 1:4
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS Chla'])
ylim([0 5000])
yyaxis right
plot(r_chl_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('Chla (ug/L)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-Chla.png']);
close all;
end

%%
figure;hold;
clearvars temp_val
temp_val = yp_w_chl_04;
temp_val(60)=[];
plot((dis_pre_total{2001-1979}) .* temp_val','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s
hold on; 
plot((dis_pre_total{2001-1979}) .* (r_chl_song(t_indx_year(nut_big_y_indx(3)-1)+1:t_indx_year(nut_big_y_indx(3))) )','r*','linew',2)
ylabel('Chla x transport (mg/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off

figure;hold; 
clearvars temp_val
temp_val = yp_w_nh4_04;
temp_val(60)=[];
plot((dis_pre_total{2001-1979}) .* temp_val','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2001-1979}) .* (r_nh4_song(t_indx_year(nut_big_y_indx(3)-1)+1:t_indx_year(nut_big_y_indx(3))) ./ 1000 .*14)','r*','linew',2)
ylabel('nh4 x transport (g/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off


figure;hold; 
clearvars temp_val
temp_val = yp_w_no3_04;
temp_val(60)=[];
plot((dis_pre_total{2001-1979}) .* temp_val','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2001-1979}) .* (r_no3_song(t_indx_year(nut_big_y_indx(3)-1)+1:t_indx_year(nut_big_y_indx(3))) ./ 1000 .*14)','r*','linew',2)
ylabel('no3 x transport (g/s)')
xlabel('days on 2001'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2001-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off

%% 2008
figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_chl_af','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_chl_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) )','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_chl_jinwal(t_07_tick(2)+1:t_07_tick(3)) )','c*','linew',2)
ylabel('Chla x transport (mg/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off

figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_nh4_af','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_nh4_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) ./ 1000 .*14)','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_nh4_jinwal(t_07_tick(2)+1:t_07_tick(3)) ./ 1000 .*14)','c*','linew',2)
ylabel('nh4 x transport (g/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off


figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_no3_af','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_no3_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) ./ 1000 .*14)','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_no3_jinwal(t_07_tick(2)+1:t_07_tick(3)) ./ 1000 .*14)','c*','linew',2)
ylabel('no3 x transport (g/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off



% set(gca,'xtick',[tt_tick(1:end)]);
% set(gca,'xticklabel',,'fontsize',10);

t_indx_year = t_indx(12:12:end); %1989~2019
%nh4, no3 exist from 1996

nut_big_y_indx = 8:length(t_indx_year)-1;
t_big_year = 1996:2018;
for i = 1:23
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS nh4'])
ylim([0 5000])
yyaxis right
plot(r_nh4_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('NH4 (mmol N /m^3)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-nh4.png']);
close all;
end

for i = 1:23
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS no3'])
ylim([0 5000])
yyaxis right
plot(r_no3_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('NO3 (mmol N /m^3)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-no3.png']);
close all;
end

for i = 1:23
figure; 
plot(double(dis_pre_total{t_big_year(i)-1979})); % m^3/s * mg/m^3 = mg/s
ylabel('transport (m^3/s)','fontsize',13)
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS Chla'])
ylim([0 5000])
yyaxis right
plot(r_chl_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))),'*')
ylabel('Chla (ug/L)','fontsize',13)
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-Chla.png']);
close all;
end

%% + jinwall all plot
%% since 2007 : 19th year 
figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_chl_af','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_chl_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) )','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_chl_jinwal(t_07_tick(2)+1:t_07_tick(3)) )','c*','linew',2)
ylabel('Chla x transport (mg/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off

figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_nh4_af','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_nh4_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) ./ 1000 .*14)','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_nh4_jinwal(t_07_tick(2)+1:t_07_tick(3)) ./ 1000 .*14)','c*','linew',2)
ylabel('nh4 x transport (g/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off


figure;hold; 
plot((dis_pre_total{2008-1979}) .* yp_w_no3_af','color','b','linew',2)  % m^3/s * mg/L = m^3/s * g/m^3 = g / s
% g/s /1000 = kg/s
hold on; 
plot((dis_pre_total{2008-1979}) .* (r_no3_song(t_indx_year(nut_big_y_indx(4)-1)+1:t_indx_year(nut_big_y_indx(4))) ./ 1000 .*14)','r*','linew',2)
hold on;
plot((dis_pre_total{2008-1979}) .* (r_no3_jinwal(t_07_tick(2)+1:t_07_tick(3)) ./ 1000 .*14)','c*','linew',2)
ylabel('no3 x transport (g/s)')
xlabel('days on 2008'); xlim([1 365]); grid on;
hold on;
yyaxis right
plot(dis_pre_total{2008-1979},'linew',2);
ylim([0 5000])
% gca.YAxis(2).Color = 'r';
ylabel('tansport(m^3/s)');
xticks([tt_tick(1:end)]);
xticklabels(1:12);
hold off


t_indx_year = t_indx(12:12:end); %1989~2019
%nh4, no3 exist from 1996
nut_big_y_indx = 19:length(t_indx_year)-1;
t_big_year = 2007:2018;
for i = 1:12
figure; 
clearvars yp_w_*_af2
if leapyear(t_big_year(i))==1
     yp_w_chl_af2 = yp_w_chl_af;  
else
     yp_w_chl_af2 = yp_w_chl_af;
     yp_w_chl_af2(60)=[];
end
    plot((dis_pre_total{t_big_year(i)-1979}) .* yp_w_chl_af2','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s

hold on;
plot((dis_pre_total{t_big_year(i)-1979}) .* (r_chl_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))) )','r*','linew',2)
hold on;
if i == 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_chl_jinwal(t_07_tick(i):t_07_tick(i+1)))','c*','linew',2)
elseif i ~= 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_chl_jinwal(t_07_tick(i)+1:t_07_tick(i+1)))','c*','linew',2)
end
hold on; grid on;
ylabel('Chla x transport (mg/s)')
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS Chla'])
yyaxis right
plot(double(dis_pre_total{t_big_year(i)-1979}),'linew',2); % m^3/s * mg/m^3 = mg/s
ylim([0 5000])
ylabel('tansport(m^3/s)');
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-Chla-jin.png']);
close all;
end

t_indx_year = t_indx(12:12:end); %1989~2019
%nh4, no3 exist from 1996
nut_big_y_indx = 19:length(t_indx_year)-1;
t_big_year = 2007:2018;
for i = 1:12
figure; 
clearvars yp_w_*_af2
if leapyear(t_big_year(i))==1
     yp_w_nh4_af2 = yp_w_nh4_af;  
else
     yp_w_nh4_af2 = yp_w_nh4_af;
     yp_w_nh4_af2(60)=[];
end
    plot((dis_pre_total{t_big_year(i)-1979}) .* yp_w_nh4_af2','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s

hold on;
plot((dis_pre_total{t_big_year(i)-1979}) .* (r_nh4_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))) ./ 1000 .* 14 )','r*','linew',2)
hold on;
if i == 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_nh4_jinwal(t_07_tick(i):t_07_tick(i+1))./ 1000 .* 14)','c*','linew',2)
elseif i ~= 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_nh4_jinwal(t_07_tick(i)+1:t_07_tick(i+1))./ 1000 .* 14)','c*','linew',2)
end
hold on;grid on;
ylim([0 300])
ylabel('nh4 x transport (mg/s)')
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS nh4'])
yyaxis right
plot(double(dis_pre_total{t_big_year(i)-1979}),'linew',2); % m^3/s * mg/m^3 = mg/s
ylim([0 5000])
ylabel('tansport(m^3/s)');
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-nh4-jin.png']);
close all;
end


t_indx_year = t_indx(12:12:end); %1989~2019
%no3, no3 exist from 1996
nut_big_y_indx = 19:length(t_indx_year)-1;
t_big_year = 2007:2018;
for i = 1:12
figure; 
clearvars yp_w_*_af2
if leapyear(t_big_year(i))==1
     yp_w_no3_af2 = yp_w_no3_af;  
else
     yp_w_no3_af2 = yp_w_no3_af;
     yp_w_no3_af2(60)=[];
end
    plot((dis_pre_total{t_big_year(i)-1979}) .* yp_w_no3_af2','color','b','linew',2)  % m^3/s * mg/m^3 = mg/s

hold on;
plot((dis_pre_total{t_big_year(i)-1979}) .* (r_no3_song(t_indx_year(nut_big_y_indx(i)-1)+1:t_indx_year(nut_big_y_indx(i))) ./ 1000 .* 14 )','r*','linew',2)
hold on;
if i == 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_no3_jinwal(t_07_tick(i):t_07_tick(i+1))./ 1000 .* 14)','c*','linew',2)
elseif i ~= 1
    plot((dis_pre_total{t_big_year(i)-1979}) .* (r_no3_jinwal(t_07_tick(i)+1:t_07_tick(i+1))./ 1000 .* 14)','c*','linew',2)
end
hold on; grid on;
ylabel('no3 x transport (mg/s)')
xlabel('month on each year','fontsize',13); xlim([1 366]); grid on;
xticks([tt_tick(1:end)]);
xticklabels(1:12);
title([num2str(t_big_year(i)) ' songjung transport vs. OBS no3'])
yyaxis right
plot(double(dis_pre_total{t_big_year(i)-1979}),'linew',2); % m^3/s * mg/m^3 = mg/s
ylim([0 5000])
ylabel('tansport(m^3/s)');
saveas(gcf,[num2str(t_big_year(i)) '-songjung-trans-no3-jin.png']);
close all;
end



%% Kd vs. CHL
for i = 1:16
    clearvars Kd r_coef plt_chl
    Kd = 1.7 ./yr_secchi_raw(i,:);
    plt_chl = yr_chl_raw(i,:);
    r_coef= corrcoef(plt_chl(~isnan(plt_chl)),Kd(~isnan(plt_chl)))
fig = figure; hold on;
        plot(Kd,'ko','linew',2);
        plot(Kd,'k','linew',2);
        
        title(['1989-2019 yearly OBS Chl. vs. Kd st.', char(name_tag{i})]);
        text(5,2.5, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        xlabel('time','fontsize',13)
        ylabel('Kd (m^-1)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 2.6])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_Kd_koem_st_',char(name_tag{i})]),'-dpng')
end



    clearvars Kd r_coef plt_chl
    Kd = 1.7 ./nanmean(yr_secchi_raw([1:6,14:16],:),1);
    plt_chl = nanmean(yr_chl_raw([1:6,14:16],:),1);
    r_coef= corrcoef(plt_chl(~isnan(plt_chl)),Kd(~isnan(plt_chl)))
fig = figure; hold on;
        plot(Kd,'ko','linew',2);
        plot(Kd,'k','linew',2);
        
        title(['1989-2019 yearly OBS Chl. vs. Kd st. 9point']);
        text(5,2.5, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        xlabel('time','fontsize',13)
        ylabel('Kd (m^-1)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0.5 1.2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 8]);
xtickangle(45)
% print(fig,strcat(['yearly_chl_Kd_koem_st_',char(name_tag{i})]),'-dpng')


% for i = 1:9
%     clearvars Kd r_coef plt_chl
%     Kd = 1.79 ./(yr_secchi_raw(i,:).^0.978) ;
%     plt_chl = yr_chl_raw(i,:);
%     r_coef= corrcoef(plt_chl(~isnan(plt_chl)),Kd(~isnan(plt_chl)))
% fig = figure; hold on;
%         plot(Kd,'b*','linew',2);
%         plot(Kd,'r','linew',2);
%         
%         title(['1989-2019 yearly OBS Chl. vs. Kd st.', nn(L(i),:)]);
%         text(5,2.5, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
%         xlabel('time','fontsize',13)
%         ylabel('Kd (m^-1)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         set(gca,'xtick',1:31);
%         xlim([1 31])
%         ylim([0 2.6])
%         set(gca,'xticklabel',1989:2019,'fontsize',10);
% yyaxis right
% plot(plt_chl,'g-','linew',2);
% ylim([0 15]);
% xtickangle(45)
% % print(fig,strcat(['yearly_chl_Kd_koem_st_',nn(L(i),:)]),'-dpng')
% end


%% SS vs. CHL
for i = 1:16
clearvars ss r_coef plt_chl coe_*
    ss = yr_ss_raw(i,:);
    plt_chl = yr_chl_raw(i,:);
    coe_ss=ss(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_ss)),coe_ss(~isnan(coe_ss)))    
fig = figure; hold on;
        plot(ss,'ko','linew',2);
        plot(ss,'k','linew',2);
        title(['1997-2018 yearly OBS Chl. vs. SS st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        text(9,35, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 45])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_SS_koem_st_',char(name_tag{i})]),'-dpng')
end

% 3p
clearvars ss r_coef plt_chl coe_*
    ss = nanmean(yr_ss_raw(6:8,21:end),1);
    plt_chl = nanmean(yr_chl_raw(6:8,21:end),1);
    coe_ss=ss(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_ss)),coe_ss(~isnan(coe_ss)))    
fig = figure; hold on;
        plot(21:size(yr_ss_raw,2),ss,'b*','linew',2);
        plot(21:size(yr_ss_raw,2),ss,'r','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. SS st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        text(5,35, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 45])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(21:size(yr_ss_raw,2),plt_chl,'g-','linew',2);
ylim([0 15]);
xtickangle(45)

clearvars ss r_coef plt_chl coe_*
    ss = nanmean(yr_ss_raw(7:8,21:end-1),1);
    plt_chl = nanmean(yr_chl_raw(7:8,21:end-1),1);
    coe_ss=ss(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_ss)),coe_ss(~isnan(coe_ss)))    
fig = figure; hold on;
        plot(21:size(yr_ss_raw,2)-1,ss,'b*','linew',2);
        plot(21:size(yr_ss_raw,2)-1,ss,'r','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. SS st. kang et al 2020']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        text(5,35, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 45])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(21:size(yr_ss_raw,2)-1,plt_chl,'g-','linew',2);
ylim([0 15]);
xtickangle(45)
% print(fig,strcat(['yearly_chl_SS_koem_st_',char(name_tag{i})]),'-dpng')

% secchi
clearvars secc r_coef plt_chl coe_*
    secc = nanmean(yr_secchi_raw(6:8,21:end-1),1);
    plt_chl = nanmean(yr_chl_raw(6:8,21:end-1),1);
    coe_secc=secc(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_secc)),coe_secc(~isnan(coe_secc)))    
fig = figure; hold on;
        plot(21:size(yr_secchi_raw,2)-1,secc,'b*','linew',2);
        plot(21:size(yr_secchi_raw,2)-1,secc,'r','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. secchi st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('secchi (m^-1)','fontsize',13)
        text(5,2, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 3])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(21:size(yr_secchi_raw,2)-1,plt_chl,'g-','linew',2);
ylim([0 15]);
xtickangle(45)



clearvars secc r_coef plt_chl coe_*
    secc = nanmean(yr_secchi_raw(7:8,21:end-1),1);
    plt_chl = nanmean(yr_chl_raw(7:8,21:end-1),1);
    coe_secc=secc(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_secc)),coe_secc(~isnan(coe_secc)))    
fig = figure; hold on;
        plot(21:size(yr_secchi_raw,2)-1,secc,'b*','linew',2);
        plot(21:size(yr_secchi_raw,2)-1,secc,'r','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. Secchi st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        text(5,35, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 3])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(21:size(yr_secchi_raw,2)-1,plt_chl,'g-','linew',2);
ylim([0 15]);
xtickangle(45)

plot(secc, plt_chl,'bo')

sp_pic = [5,7:9,13];
clearvars secc r_coef plt_chl coe_*
    secc = nanmean(yr_secchi_raw(sp_pic,21:end-1),1);
    plt_chl = nanmean(yr_chl_raw(6:8,21:end-1),1);
    coe_secc=secc(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_secc)),coe_secc(~isnan(coe_secc)))    
fig = figure; hold on;
        plot(21:size(yr_secchi_raw,2)-1,secc,'b*','linew',2);
        plot(21:size(yr_secchi_raw,2)-1,secc,'r','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. Secchi st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        text(5,35, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([1 3])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(21:size(yr_secchi_raw,2)-1,plt_chl,'g-','linew',2);
ylim([1 8]);
xtickangle(45)


sp_pic = [5,7:9,13];
clearvars secc r_coef plt_chl coe_*
    secc = nanmean(yr_secchi_raw(sp_pic,21:end-1),1);
    plt_chl = nanmean(yr_chl_raw(6:8,21:end-1),1);
    coe_secc=secc(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_secc)),coe_secc(~isnan(coe_secc)))    
fig = figure; hold on;
        plot(21:size(yr_secchi_raw,2)-1,1.7 ./secc,'b*','linew',2);
        plot(21:size(yr_secchi_raw,2)-1,1.7 ./secc,'r','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. Secchi st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        text(5,35, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
%         ylim([1 3])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(21:size(yr_secchi_raw,2)-1,plt_chl,'g-','linew',2);
ylim([1 8]);
xtickangle(45)

%% PO4
for i = 1:16
clearvars po4 r_coef plt_chl coe_*
    po4 = yr_po4_raw(i,:);
    plt_chl = yr_chl_raw(i,:);
    coe_po4=po4(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_po4)),coe_po4(~isnan(coe_po4)))    
fig = figure; hold on;
        plot(po4,'ko','linew',2);
        plot(po4,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. PO4 st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('PO4 (uM)','fontsize',13)
        text(9,.5, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        if i==1
%             ylim([0 8])
            ylim([0 5])
        else 
           ylim([0 2.5])
        end
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_PO4_koem_st_',char(name_tag{i})]),'-dpng')
end

%% DIN(no3+Nh4) : PO4 vs. chl
for i = 1:16
clearvars NP po4 din r_coef plt_chl coe_*
    po4 = yr_po4_raw(i,:);
    din = yr_din_raw(i,:);
    NP = din./po4;
    plt_chl = yr_chl_raw(i,:);
    coe_NP=NP(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_NP)),coe_NP(~isnan(coe_NP)))    
fig = figure; hold on;
        plot(NP,'ko','linew',2);
        plot(NP,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. NP ratio st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('NP ratio ','fontsize',13)
        text(9,20, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 45])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_NP_koem_st_',char(name_tag{i})]),'-dpng')
end


%% SI : PO4 vs. chl
for i = 1:16
clearvars siP po4 si r_coef plt_chl coe_*
    po4 = yr_po4_raw(i,:);
    si = yr_si_raw(i,:);
    siP = (si)./po4;
    plt_chl = yr_chl_raw(i,:);
    coe_siP=siP(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_siP)),coe_siP(~isnan(coe_siP)))    
fig = figure; hold on;
        plot(siP,'ko','linew',2);
        plot(siP,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. siP ratio st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('siP ratio ','fontsize',13)
        text(9,20, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 80])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_siP_koem_st_',char(name_tag{i})]),'-dpng')
end

%% Nsi (DIN : Silicate) vs. chl
for i = 1:16
clearvars Nsi din si r_coef plt_chl coe_*
    din = yr_din_raw(i,:);
    si = yr_si_raw(i,:);
    Nsi = din./si;
    plt_chl = yr_chl_raw(i,:);
    coe_Nsi=Nsi(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_Nsi)),coe_Nsi(~isnan(coe_Nsi)))    
fig = figure; hold on;
        plot(Nsi,'ko','linew',2);
        plot(Nsi,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. Nsi ratio st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('Nsi ratio ','fontsize',13)
        text(9,1, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 2])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_Nsi_koem_st_',char(name_tag{i})]),'-dpng')
end

for i = 2003:2019
t_text{i-2002}={num2str(i)}
end
%% siN (Silicate:DIN) vs. chl
for i = 1:16
clearvars Nsi din si r_coef plt_chl coe_*
    din = yr_din_raw(i,:);
    si = yr_si_raw(i,:);
    Nsi = si./din;
    plt_chl = yr_chl_raw(i,:);
    clearvars NP po4 r_coef  coe_*
    po4 = yr_po4_raw(i,:);
    din = yr_din_raw(i,:);
    NP = din./po4;  
fig = figure; hold on;
        plot(Nsi(1:15),NP(1:15),'c.','linew',2);
        plot(Nsi(16:end),NP(16:end),'k.','linew',2);
        for j= 1:16
        text(Nsi(j+15),NP(j+15),t_text{j})
        end
        title(['1989-2019 yearly OBS NP. vs. siN ratio st.', char(name_tag{i})]);
        xlabel('siN ratio','fontsize',13)
        ylabel('NP ratio ','fontsize',13)
        xline(1,'b--','linew',2);
        yline(16,'g--','linew',2);
        plot(0:0.01:10,(0/16:0.01/16:10/16).^-1,'r--','linew',2);
        grid on
        set(gca,'fontsize',13)
%         ylim([0 100])
%         xlim([0 10])
%         set(gca, 'YScale', 'log')
%         set(gca, 'XScale', 'log')        
        ylim([0 50])
        xlim([0 6])
print(fig,strcat(['yearly_NP_siN_koem_st_',char(name_tag{i})]),'-dpng')
end


%% siN (Silicate:DIN) vs. NP ratio
for i = 1:16
clearvars Nsi din si r_coef plt_chl coe_*
    din = yr_din_raw(i,:);
    si = yr_si_raw(i,:);
    Nsi = si./din;
    plt_chl = yr_chl_raw(i,:);
    coe_Nsi=Nsi(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_Nsi)),coe_Nsi(~isnan(coe_Nsi)))    
fig = figure; hold on;
        plot(Nsi,'ko','linew',2);
        plot(Nsi,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. siN ratio st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('siN ratio ','fontsize',13)
        text(9,3.5, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 5])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_siN_koem_st_',char(name_tag{i})]),'-dpng')
end

%% NH4
for i = 1:16
clearvars nh4 r_coef plt_chl coe_*
    nh4 = yr_nh4_raw(i,:);
    plt_chl = yr_chl_raw(i,:);
    coe_nh4=nh4(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_nh4)),coe_nh4(~isnan(coe_nh4)))    
fig = figure; hold on;
        plot(nh4,'ko','linew',2);
        plot(nh4,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. nh4 st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('nh4 (uM)','fontsize',13)
        text(9,5, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 8])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_nh4_koem_st_',char(name_tag{i})]),'-dpng')
end

%% DIN
for i = 1:16
clearvars din r_coef plt_chl coe_*
    din = yr_din_raw(i,:);
    plt_chl = yr_chl_raw(i,:);
    coe_din=din(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_din)),coe_din(~isnan(coe_din)))    
fig = figure; hold on;
        plot(din,'ko','linew',2);
        plot(din,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. din st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('DIN (uM)','fontsize',13)
        text(25,18, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 23.5])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_din_koem_st_',char(name_tag{i})]),'-dpng')
end

%% no3
for i = 1:16
clearvars no3 r_coef plt_chl coe_*
    no3 = yr_no3_raw(i,:);
    plt_chl = yr_chl_raw(i,:);
    coe_no3=no3(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_no3)),coe_no3(~isnan(coe_no3)))    
fig = figure; hold on;
        plot(no3,'ko','linew',2);
        plot(no3,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. no3 st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('no3 (uM)','fontsize',13)
        text(9,18, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 20])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_no3_koem_st_',char(name_tag{i})]),'-dpng')
end

%% SI
for i = 1:16
clearvars si r_coef plt_chl coe_*
    si = yr_si_raw(i,:);
    plt_chl = yr_chl_raw(i,:);
    coe_si=si(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_si)),coe_si(~isnan(coe_si)))    
fig = figure; hold on;
        plot(si,'ko','linew',2);
        plot(si,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. si st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('si (uM)','fontsize',13)
        text(9,20, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_si_koem_st_',char(name_tag{i})]),'-dpng')
end

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(1.7 ./yr_secchi,'b*','linew',2);
        plot(1.7 ./yr_secchi,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. Kd']);
        xlabel('time','fontsize',13)
        ylabel('Kd (m^-1)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)

% chl vs. transp
for i = 1:16
clearvars trans r_coef plt_chl coe_*
    trans = yr_trans;
    plt_chl = yr_chl_raw(i,:);
    coe_trans=trans(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_trans)),coe_trans(~isnan(coe_trans)))    
fig = figure; hold on;
        plot(trans,'ko','linew',2);
        plot(trans,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. trans st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('transport (m^3/s)','fontsize',13)
        text(5,4, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 inf])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_trans_koem_st_',char(name_tag{i})]),'-dpng')
end

% TEMP vs. chl
for i = 1:16
clearvars temp r_coef plt_chl coe_*
    temp= yr_temp_raw(i,:);
    plt_chl = yr_chl_raw(i,:);
    coe_temp=temp(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_temp)),coe_temp(~isnan(coe_temp)))    
fig = figure; hold on;
        plot(temp,'ko','linew',2);
        plot(temp,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. temp st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('temp (^oC)','fontsize',13)
        text(9,19, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([14 20])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_temp_koem_st_',char(name_tag{i})]),'-dpng')
end

% SS vs. chl
for i = 1:16
clearvars ss r_coef plt_chl coe_*
    ss = yr_ss_raw(i,:);
    plt_chl = yr_chl_raw(i,:);
    coe_ss=ss(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_ss)),coe_ss(~isnan(coe_ss)))    
fig = figure; hold on;
        plot(ss,'ko','linew',2);
        plot(ss,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. ss st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('ss (mg/L)','fontsize',13)
        text(9,4, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_ss_koem_st_',char(name_tag{i})]),'-dpng')
end

% SS vs. Kd
for i = 1:16
clearvars ss r_coef secc coe_*
    ss = yr_ss_raw(i,:);
    secc = 1.7 ./ yr_secchi_raw(i,:);
    coe_ss=ss(~isnan(secc)); coe_secc=secc(~isnan(secc));
    r_coef= corrcoef(coe_secc(~isnan(coe_ss)),coe_ss(~isnan(coe_ss)))    
fig = figure; hold on;
        plot(ss,'ko','linew',2);
        plot(ss,'k','linew',2);
        title(['1989-2019 yearly KOEM Kd. vs. SS st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        text(9,4, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(secc,'b-','linew',2);
ylim([0 2.5]);
xtickangle(45)
print(fig,strcat(['yearly_Kd_ss_koem_st_',char(name_tag{i})]),'-dpng')
end

% SS vs secchi
for i = 1:16
clearvars ss r_coef secc coe_*
    ss = yr_ss_raw(i,:);
    secc = yr_secchi_raw(i,:);
    coe_ss=ss(~isnan(secc)); coe_secc=secc(~isnan(secc));
    r_coef= corrcoef(coe_secc(~isnan(coe_ss)),coe_ss(~isnan(coe_ss)))    
fig = figure; hold on;
        plot(ss,'ko','linew',2);
        plot(ss,'k','linew',2);
        title(['1989-2019 yearly KOEM secchi. vs. SS st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        text(9,4, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(secc,'b-','linew',2);
ylim([0 5]);
xtickangle(45)
print(fig,strcat(['yearly_secchi_ss_koem_st_',char(name_tag{i})]),'-dpng')
end

% trans vs. SS
for i = 1:16
clearvars trans r_coef secc coe_*
    ss = yr_ss_raw(i,:);
    trans = yr_trans;
    coe_ss=ss(~isnan(trans)); coe_trans=trans(~isnan(trans));
    r_coef= corrcoef(coe_trans(~isnan(coe_ss)),coe_ss(~isnan(coe_ss)))    
fig = figure; hold on;
        plot(ss,'ko','linew',2);
        plot(ss,'k','linew',2);
        title(['1989-2019 yearly KOEM transp. vs. SS st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        text(5,4, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 30])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(trans,'b-','linew',2);
ylim([0 inf]);
xtickangle(45)
print(fig,strcat(['yearly_transp_ss_koem_st_',char(name_tag{i})]),'-dpng')
end



% Kd
for i = 1:16
clearvars secc r_coef plt_chl coe_* Kd_c idx_*
    secc = yr_secchi_raw(i,:);
    idx_1=find(secc < 2.20);
    idx_2=find(secc >= 2.20);
    Kd_c(idx_1) = 1.16./((secc(idx_1)).^0.62);
    Kd_c(idx_2) = exp((0.15-log(secc(idx_2)).*0.62) .* (1.68-log(secc(idx_2)))./0.89  ...
    + (-0.48-log(secc(idx_2)).*0.72).* (log(secc(idx_2))-0.79)./0.89);
    Kd_c(Kd_c==0)=NaN;
    plt_chl = yr_chl_raw(i,:);
    coe_secc=Kd_c(~isnan(plt_chl)); coe_chl=plt_chl(~isnan(plt_chl));
    r_coef= corrcoef(coe_chl(~isnan(coe_secc)),coe_secc(~isnan(coe_secc)))    
fig = figure; hold on;
        plot(Kd_c,'ko','linew',2);
        plot(Kd_c,'k','linew',2);
        title(['1989-2019 yearly OBS Chl. vs. Kd st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        text(9,2.8, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 1.4])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(plt_chl,'b-','linew',2);
ylim([0 15]);
xtickangle(45)
print(fig,strcat(['yearly_chl_Kdc_koem_st_',char(name_tag{i})]),'-dpng')
end

% transp vs. kd
for i = 1:16
clearvars secc r_coef trans coe_*
    secc = 1.7 ./ yr_secchi_raw(i,:);
     trans = yr_trans;
    coe_secc=secc(~isnan(trans)); coe_trans=trans(~isnan(trans));
    r_coef= corrcoef(coe_trans(~isnan(coe_secc)),coe_secc(~isnan(coe_secc)))    
fig = figure; hold on;
        plot(secc,'ko','linew',2);
        plot(secc,'k','linew',2);
        title(['1989-2019 yearly OBS transp. vs. Kd st.', char(name_tag{i})]);
        xlabel('time','fontsize',13)
        ylabel('Kd (m^-1)','fontsize',13)
        text(1,.2, ['r = ',num2str(r_coef(1,2),'%0.2f')],'fontsize',15)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        ylim([0 2.5])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
yyaxis right
plot(trans,'b-','linew',2);
ylim([0 inf]);
xtickangle(45)
print(fig,strcat(['yearly_transp_Kd_koem_st_',char(name_tag{i})]),'-dpng')
end


%% SS
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_ss,'b*','linew',2);
        plot(yr_ss,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)






%% Si
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_si,'b*','linew',2);
        plot(yr_si,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. Si']);
        xlabel('time','fontsize',13)
        ylabel('Si (mmol si / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end

%% secchi
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_secchi,'b*','linew',2);
        plot(yr_secchi,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. secchi']);
        xlabel('time','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(1.7 ./yr_secchi,'b*','linew',2);
        plot(1.7 ./yr_secchi,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. Kd']);
        xlabel('time','fontsize',13)
        ylabel('Kd (m^-1)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)

%% chl
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_chl,'b*','linew',2);
        plot(yr_chl,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2);
xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end

%% NH4
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_nh4,'b*','linew',2);
        plot(yr_nh4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. nh4']);
        xlabel('time','fontsize',13)
        ylabel('NH4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_no3,'b*','linew',2);
        plot(yr_no3,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. no3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)
% for i = 6:12:360; xline(i,'color','m'); end


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_po4,'b*','linew',2);
        plot(yr_po4,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1989-2019 yearly OBS transp vs. PO4']);
        xlabel('time','fontsize',13)
        ylabel('PO4 (mmol P / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([1 31])
        set(gca,'xticklabel',1989:2019,'fontsize',10);
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
yyaxis right
plot(yr_trans,'g-','linew',2); xtickangle(45)


%% KOEM - 9 st timeseries

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_temp_9p,'b*','linew',2);
        plot(yr_temp_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS temperature']);
        xlabel('time','fontsize',13)
        ylabel('temperature (^oC)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_nh4_9p,'b*','linew',2);
        plot(yr_nh4_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS NH4-N']);
        xlabel('time','fontsize',13)
        ylabel('NH4-N (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 5])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_no3_9p,'b*','linew',2);
        plot(yr_no3_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS NO3-N']);
        xlabel('time','fontsize',13)
        ylabel('NO3-N (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 17])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_din_9p,'b*','linew',2);
        plot(yr_din_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS DIN']);
        xlabel('time','fontsize',13)
        ylabel('DIN (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 17])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_po4_9p,'b*','linew',2);
        plot(yr_po4_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS PO4-P']);
        xlabel('time','fontsize',13)
        ylabel('PO4-P (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 2.0])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_si_9p,'b*','linew',2);
        plot(yr_si_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS SI']);
        xlabel('time','fontsize',13)
        ylabel('si (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_temp_b_9p,'b*','linew',2);
        plot(yr_temp_b_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS bottom temperature']);
        xlabel('time','fontsize',13)
        ylabel('temperature (^oC)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_nh4_b_9p,'b*','linew',2);
        plot(yr_nh4_b_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS bottom NH4-N']);
        xlabel('time','fontsize',13)
        ylabel('NH4-N (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 5])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_no3_b_9p,'b*','linew',2);
        plot(yr_no3_b_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS bottom NO3-N']);
        xlabel('time','fontsize',13)
        ylabel('NO3-N (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 17])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_din_b_9p,'b*','linew',2);
        plot(yr_din_b_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS bottom DIN']);
        xlabel('time','fontsize',13)
        ylabel('DIN (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 17])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_po4_b_9p,'b*','linew',2);
        plot(yr_po4_b_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS bottom PO4-P']);
        xlabel('time','fontsize',13)
        ylabel('PO4-P (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        ylim([0 2.0])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_si_9p,'b*','linew',2);
        plot(yr_si_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS bottom SI']);
        xlabel('time','fontsize',13)
        ylabel('si (mmol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_salt_9p,'b*','linew',2);
        plot(yr_salt_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS salt']);
        xlabel('time','fontsize',13)
        ylabel('salt (PSU)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_salt_b_9p,'b*','linew',2);
        plot(yr_salt_b_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS bottom salt']);
        xlabel('time','fontsize',13)
        ylabel('salt (PSU)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_chl_9p,'b*','linew',2);
        plot(yr_chl_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS Chl']);
        xlabel('time','fontsize',13)
        ylabel('Chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_chl_b_9p,'b*','linew',2);
        plot(yr_chl_b_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS bottom Chl']);
        xlabel('time','fontsize',13)
        ylabel('Chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_secchi_9p,'b*','linew',2);
        plot(yr_secchi_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS secchi depth']);
        xlabel('time','fontsize',13)
        ylabel('secchi depth (m)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 3.2])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_ss_9p,'b*','linew',2);
        plot(yr_ss_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 20])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(yr_ss_b_9p,'b*','linew',2);
        plot(yr_ss_b_9p,'r','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997-2019 yearly OBS bottom SS']);
        xlabel('time','fontsize',13)
        ylabel('SS (mg/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',1:31);
        xlim([9 31])
%         ylim([0 20])
        set(gca,'xticklabel',1989:2020,'fontsize',10);
xtickangle(45)