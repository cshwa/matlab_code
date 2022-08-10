close all; clear all; clc; 
cd D:\장기생태\Dynamic\06_river\data\Namgang_dam


t1 = datetime(2000,01,01,00,00,00); %>> 21,00,00 (GMT)
t2 = datetime(2020,12,31,00,00,00); %>> 22,00,00  (GMT)
t = t1:days(1):t2;
time = datenum(t);
clearvars total_date total_mmdd
for k=1:length(time)
total_date{k,1}=datestr(time(k),'yyyy-mm-dd'); 
total_mmdd{k,1}=datestr(time(k),'mm-dd');
end

% t_day_cut1=365*length(1989:2006);
% t_day_cut2=365*length(1989:2015);
% 2007-01-01~2015-12-31: 2558:5844
% 2016-01-01~2019-12-31: 5845:7305

m=0;
for k=2558:5844
m=m+1;
first_mmdd{m,1}=datestr(time(k),'mm-dd');
end
m=0;
for k=5845:7305
m=m+1;
second_mmdd{m,1}=datestr(time(k),'mm-dd');
end


clearvars ref_365_mmdd
t1 = datetime(2001,01,01,00,00,00); %>> 21,00,00 (GMT)
t2 = datetime(2001,12,31,00,00,00); %>> 22,00,00  (GMT)
t = t1:days(1):t2;
time = datenum(t);
ref_365_mmdd=datestr(time,'mm-dd');

clearvars indx_365
for i = 1:length(ref_365_mmdd)
   if  sum(strcmp(ref_365_mmdd(i,:),total_mmdd)) ~= 0
       indx_365{i} = find([strcmp(ref_365_mmdd(i,:),total_mmdd)] == 1);     
   end
end

clearvars indx_365_1st
for i = 1:length(ref_365_mmdd)
   if  sum(strcmp(ref_365_mmdd(i,:),total_mmdd)) ~= 0
       indx_365_1st{i} = find([strcmp(ref_365_mmdd(i,:),first_mmdd)] == 1);     
   end
end

clearvars indx_365_2nd
for i = 1:length(ref_365_mmdd)
   if  sum(strcmp(ref_365_mmdd(i,:),total_mmdd)) ~= 0
       indx_365_2nd{i} = find([strcmp(ref_365_mmdd(i,:),second_mmdd)] == 1);     
   end
end


for i = 2000:2020
    clearvars raw txt temp_*
[raw txt]=xlsread(strcat('namgang_dam.xlsx'),num2str(i));

% reference time
t1 = datetime(i,01,01,00,00,00); %>> 21,00,00 (GMT)
t2 = datetime(i,12,31,00,00,00); %>> 22,00,00  (GMT)
t = t1:days(1):t2;
time = datenum(t);
clearvars ref_date
for k=1:length(time)
ref_date{k,1}=datestr(time(k),'yyyy-mm-dd'); 
end

if i > 2018
    % data's date
    date=txt(2:end,1);

    %% pick matched name with tag
% for j = 1:length(ref_date)
%    if  sum(strcmp(ref_date(j,:), date)) ~= 0
%        indx{j} = find([strcmp(ref_date(j,:), date)] == 1)     
%    end
% end
clearvars indx
for j = 1:length(date)
   if  sum(strcmp(date{j},ref_date)) ~= 0
       indx(j) = find([strcmp(date{j},ref_date)] == 1);     
   end
end

    temp_prec = NaN(length(ref_date),1);
    temp_inflow = NaN(length(ref_date),1);
    temp_inflow_self = NaN(length(ref_date),1);
    temp_outflow_total = NaN(length(ref_date),1);
    temp_outflow_power = NaN(length(ref_date),1);
    temp_outflow_illryu = NaN(length(ref_date),1);
    temp_outflow_power2 = NaN(length(ref_date),1);
    temp_outflow_gawha = NaN(length(ref_date),1);
    temp_outflow_sacheon = NaN(length(ref_date),1);
    temp_outflow_namgang = NaN(length(ref_date),1);
    temp_outflow_jinju = NaN(length(ref_date),1);
    temp_outflow_namdong = NaN(length(ref_date),1);
    temp_outflow_power3 = NaN(length(ref_date),1);
   

    temp_prec(indx) = [ raw(:,5) ];
    temp_inflow(indx) = [ raw(:,6) ];
    temp_inflow_self(indx) = [ raw(:,7) ];
    temp_outflow_total(indx) = [ raw(:,8) ];
    temp_outflow_power(indx) = [ raw(:,9) ];
    temp_outflow_illryu(indx) = [ raw(:,10) ];
    temp_outflow_power2(indx) = [ raw(:,11) ];
    temp_outflow_gawha(indx) = [ raw(:,12) ];
    temp_outflow_sacheon(indx) = [ raw(:,13) ];
    temp_outflow_namgang(indx) = [ raw(:,14) ];
    temp_outflow_jinju(indx) = [ raw(:,15) ];
    temp_outflow_namdong(indx) = [ raw(:,16) ];
    temp_outflow_power3(indx) = [ raw(:,17) ];
    
    if i==2019        
        merg_prec = [merg_prec; temp_prec;];
        merg_inflow = [merg_inflow; temp_inflow;];
        merg_inflow_self = [merg_inflow_self; temp_inflow_self;];
        merg_outflow_total = [merg_outflow_total; temp_outflow_total;];
        merg_outflow_power = [merg_outflow_power; temp_outflow_power;];
        merg_outflow_illryu = [merg_outflow_illryu; temp_outflow_illryu;];
        merg_outflow_power2 = [ temp_outflow_power2 ];
        merg_outflow_gawha = [merg_outflow_gawha; temp_outflow_gawha;];
        merg_outflow_sacheon = [ temp_outflow_sacheon ];
        merg_outflow_namgang = [merg_outflow_namgang; temp_outflow_namgang;];
        merg_outflow_jinju = [merg_outflow_jinju; temp_outflow_jinju;];
        merg_outflow_namdong = [merg_outflow_namdong; temp_outflow_namdong;];
        merg_outflow_power3 = [ temp_outflow_power3 ];
    else
    
        merg_prec = [merg_prec; temp_prec;];
        merg_inflow = [merg_inflow; temp_inflow;];
        merg_inflow_self = [merg_inflow_self; temp_inflow_self;];
        merg_outflow_total = [merg_outflow_total; temp_outflow_total;];
        merg_outflow_power = [merg_outflow_power; temp_outflow_power;];
        merg_outflow_illryu = [merg_outflow_illryu; temp_outflow_illryu;];
        merg_outflow_power2 = [merg_outflow_power2; temp_outflow_power2;];
        merg_outflow_gawha = [merg_outflow_gawha; temp_outflow_gawha;];
        merg_outflow_sacheon = [merg_outflow_sacheon; temp_outflow_sacheon;];
        merg_outflow_namgang = [merg_outflow_namgang; temp_outflow_namgang;];
        merg_outflow_jinju = [merg_outflow_jinju; temp_outflow_jinju;];
        merg_outflow_namdong = [merg_outflow_namdong; temp_outflow_namdong;];
        merg_outflow_power3 = [merg_outflow_power3; temp_outflow_power3;];
    end
else
    % data's date
    date=txt(3:end,1);
    
    %% pick matched name with tag
% for j = 1:length(ref_date)
%    if  sum(strcmp(ref_date(j,:), date)) ~= 0
%        indx{j} = find([strcmp(ref_date(j,:), date)] == 1)     
%    end
% end
clearvars indx
for j = 1:length(date)
   if  sum(strcmp(date{j},ref_date)) ~= 0
       indx(j) = find([strcmp(date{j},ref_date)] == 1);     
   end
end

    temp_prec = NaN(length(ref_date),1);
    temp_inflow = NaN(length(ref_date),1);
    temp_inflow_self = NaN(length(ref_date),1);
    temp_outflow_total = NaN(length(ref_date),1);
    temp_outflow_power = NaN(length(ref_date),1);
    temp_outflow_illryu = NaN(length(ref_date),1);
    temp_outflow_gawha = NaN(length(ref_date),1);
    temp_outflow_namgang = NaN(length(ref_date),1);
    temp_outflow_jinju = NaN(length(ref_date),1);
    temp_outflow_namdong = NaN(length(ref_date),1);

    temp_prec(indx) = [ raw(:,5) ];
    temp_inflow(indx) = [ raw(:,6) ];
    temp_inflow_self(indx) = [ raw(:,7) ];
    temp_outflow_total(indx) = [ raw(:,8) ];
    temp_outflow_power(indx) = [ raw(:,9) ];
    temp_outflow_illryu(indx) = [ raw(:,10) ];
    temp_outflow_gawha(indx) = [ raw(:,11) ];
    temp_outflow_namgang(indx) = [ raw(:,12) ];
    temp_outflow_jinju(indx) = [ raw(:,13) ];
    temp_outflow_namdong(indx) = [ raw(:,14) ];


if i == 2000
    merg_prec = temp_prec;
    merg_inflow = temp_inflow;
    merg_inflow_self = temp_inflow_self;
    merg_outflow_total = temp_outflow_total;
    merg_outflow_power = temp_outflow_power;
    merg_outflow_illryu = temp_outflow_illryu;
    merg_outflow_gawha = temp_outflow_gawha;
    merg_outflow_namgang = temp_outflow_namgang;
    merg_outflow_jinju = temp_outflow_jinju;
    merg_outflow_namdong = temp_outflow_namdong;
else
    merg_prec = [merg_prec; temp_prec;];
    merg_inflow = [merg_inflow; temp_inflow;];
    merg_inflow_self = [merg_inflow_self; temp_inflow_self;];
    merg_outflow_total = [merg_outflow_total; temp_outflow_total;];
    merg_outflow_power = [merg_outflow_power; temp_outflow_power;];
    merg_outflow_illryu = [merg_outflow_illryu; temp_outflow_illryu;];
    merg_outflow_gawha = [merg_outflow_gawha; temp_outflow_gawha;];
    merg_outflow_namgang = [merg_outflow_namgang; temp_outflow_namgang;];
    merg_outflow_jinju = [merg_outflow_jinju; temp_outflow_jinju;];
    merg_outflow_namdong = [merg_outflow_namdong; temp_outflow_namdong;];
end
end
end

%% missing value interped
x_interped = 1:length(total_date);
        merg_prec(isnan(merg_prec)) = interp1(x_interped(~isnan(merg_prec)),merg_prec(~isnan(merg_prec)),x_interped(isnan(merg_prec)));
        merg_inflow(isnan(merg_inflow)) = interp1(x_interped(~isnan(merg_inflow)),merg_inflow(~isnan(merg_inflow)),x_interped(isnan(merg_inflow)));
        merg_inflow_self(isnan(merg_inflow_self)) = interp1(x_interped(~isnan(merg_inflow_self)),merg_inflow_self(~isnan(merg_inflow_self)),x_interped(isnan(merg_inflow_self)));
        merg_outflow_total(isnan(merg_outflow_total)) = interp1(x_interped(~isnan(merg_outflow_total)),merg_outflow_total(~isnan(merg_outflow_total)),x_interped(isnan(merg_outflow_total)));
        merg_outflow_power(isnan(merg_outflow_power)) = interp1(x_interped(~isnan(merg_outflow_power)),merg_outflow_power(~isnan(merg_outflow_power)),x_interped(isnan(merg_outflow_power)));
        merg_outflow_illryu(isnan(merg_outflow_illryu)) = interp1(x_interped(~isnan(merg_outflow_illryu)),merg_outflow_illryu(~isnan(merg_outflow_illryu)),x_interped(isnan(merg_outflow_illryu)));
        merg_outflow_power2(isnan(merg_outflow_power2)) = interp1(x_interped(~isnan(merg_outflow_power2)),merg_outflow_power2(~isnan(merg_outflow_power2)),x_interped(isnan(merg_outflow_power2)));
        merg_outflow_gawha(isnan(merg_outflow_gawha)) = interp1(x_interped(~isnan(merg_outflow_gawha)),merg_outflow_gawha(~isnan(merg_outflow_gawha)),x_interped(isnan(merg_outflow_gawha)));
        merg_outflow_sacheon(isnan(merg_outflow_sacheon)) = interp1(x_interped(~isnan(merg_outflow_sacheon)),merg_outflow_sacheon(~isnan(merg_outflow_sacheon)),x_interped(isnan(merg_outflow_sacheon)));
        merg_outflow_namgang(isnan(merg_outflow_namgang)) = interp1(x_interped(~isnan(merg_outflow_namgang)),merg_outflow_namgang(~isnan(merg_outflow_namgang)),x_interped(isnan(merg_outflow_namgang)));
        merg_outflow_jinju(isnan(merg_outflow_jinju)) = interp1(x_interped(~isnan(merg_outflow_jinju)),merg_outflow_jinju(~isnan(merg_outflow_jinju)),x_interped(isnan(merg_outflow_jinju)));
        merg_outflow_namdong(isnan(merg_outflow_namdong)) = interp1(x_interped(~isnan(merg_outflow_namdong)),merg_outflow_namdong(~isnan(merg_outflow_namdong)),x_interped(isnan(merg_outflow_namdong)));
        merg_outflow_power3(isnan(merg_outflow_power3)) = interp1(x_interped(~isnan(merg_outflow_power3)),merg_outflow_power3(~isnan(merg_outflow_power3)),x_interped(isnan(merg_outflow_power3)));

       merg_inflow_self(1) = merg_inflow_self(2);
       merg_outflow_total(1) = merg_outflow_total(2);
       merg_outflow_power(1) = merg_outflow_power(2);
       merg_outflow_illryu(1) = merg_outflow_illryu(2);
%        merg_outflow_power2(1) = merg_outflow_power2(2);
       merg_outflow_gawha(1) = merg_outflow_gawha(2);
%        merg_outflow_sacheon(1) = merg_outflow_sacheon(2);
       merg_outflow_namgang(1) = merg_outflow_namgang(2);
       merg_outflow_jinju(1) = merg_outflow_jinju(2);
       merg_outflow_namdong(1) = merg_outflow_namdong(2);
%        merg_outflow_power3(1) = merg_outflow_power3(2);

% correlation coefficients 
r_dis_illryu_gawha=corrcoef(merg_outflow_illryu,merg_outflow_gawha); % correlation coefficients between illryu and gawha
r_prec_gawha=corrcoef(merg_prec,merg_outflow_gawha);
r_prec_illryu=corrcoef(merg_prec,merg_outflow_illryu);
r_prec_totalout=corrcoef(merg_prec,merg_outflow_total);
r_prec_inflow=corrcoef(merg_prec,merg_inflow);


%climate
for i = 1:length(indx_365)
climate_out_illryu(i)=nanmean(merg_outflow_illryu(indx_365{i}));
climate_out_gawha(i)=nanmean(merg_outflow_gawha(indx_365{i}));
end

% 2007-01-01~2015-12-31: 2558:5844
% 2016-01-01~2019-01-01: 5845:7305

merg_outflow_illryu_1st = merg_outflow_illryu(2558:5844); % 2007-01-01~2015-12-31: 2558:5844
merg_outflow_illryu_2nd = merg_outflow_illryu(5845:end); % 2016-01-01~2020-12-31: 5845:7305
merg_outflow_gawha_1st = merg_outflow_gawha(2558:5844); % 2007-01-01~2015-12-31: 2558:5844
merg_outflow_gawha_2nd = merg_outflow_gawha(5845:end); % 2016-01-01~2020-12-31: 5845:7305

%climate
for i = 1:length(indx_365_1st)
climate_out_illryu_1st(i)=nanmean(merg_outflow_illryu_1st(indx_365_1st{i}));
climate_out_gawha_1st(i)=nanmean(merg_outflow_gawha_1st(indx_365_1st{i}));
end

for i = 1:length(indx_365_2nd)
climate_out_illryu_2nd(i)=nanmean(merg_outflow_illryu_2nd(indx_365_2nd{i}));
climate_out_gawha_2nd(i)=nanmean(merg_outflow_gawha_2nd(indx_365_2nd{i}));
end


figure; hold on;
plot(climate_out_gawha,'r');
plot(climate_out_illryu,'b');
r_dis_illryu_gawha=corrcoef(climate_out_illryu,climate_out_gawha); % correlation coefficients between illryu and gawha

% 1st and 2nd climate
figure; hold on;
plot(climate_out_gawha_2nd,'r');
plot(climate_out_illryu_2nd,'b');
r_dis_illryu_gawha_2nd=corrcoef(climate_out_illryu_2nd,climate_out_gawha_2nd) % correlation coefficients between illryu and gawha

figure; hold on;
plot(climate_out_gawha_1st,'r');
plot(climate_out_illryu_1st,'b');
r_dis_illryu_gawha_1st=corrcoef(climate_out_illryu_1st,climate_out_gawha_1st) % correlation coefficients between illryu and gawha

save('gawha_gate_discharge_and_compare_with_namgang_2020.mat');

figure; 
plot(temp_inflow)
plot(merg_inflow)
hold on; plot(merg_outflow_total,'r')

%% outflow compare
figure; hold on;
plot(merg_outflow_total,'b');
plot(merg_outflow_power + merg_outflow_illryu + merg_outflow_gawha + merg_outflow_namgang + merg_outflow_jinju + merg_outflow_namdong,'r');

figure; hold on;
plot(merg_outflow_gawha,'r');
plot(merg_outflow_illryu,'b');
plot(merg_outflow_namgang,'g');
plot(merg_outflow_jinju,'k'); 
plot(merg_outflow_namdong,'c'); 
plot(merg_outflow_power,'color',[.5 .5 .5]);


figure; 
plot(merg_outflow_total - (merg_outflow_power + merg_outflow_illryu + merg_outflow_gawha + merg_outflow_namgang + merg_outflow_jinju + merg_outflow_namdong),'c');

figure; plot(merg_inflow-merg_outflow_total)


%% Compare with songjung dis
close all; clear; clc; 
%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
load('D:\장기생태\Dynamic\06_river\data\sj_1980to1996\songjung_discharge_1980to2018.mat');  % pre_merg_dis is discharge
%climate
sj = load('D:\장기생태\Dynamic\06_river\data\sj_1980to1996\songjung_climate_days_bio_to07to16_high_airT.mat');
nam = load('D:\장기생태\Dynamic\06_river\data\Namgang_dam\gawha_gate_discharge_and_compare_with_namgang.mat');

% load excel 2019 trans
[raw txt]=xlsread('D:\장기생태\Dynamic\06_river\환경과학원\수리·수문·기상_유량_구례_20190101-20191231.xls','유량');

t_year = 1989:2018

clearvars tempo_trans sj_trans_out_365
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
% tempo_temp=sumjin_re_w_c{t_year(order_i)-1989};
tempo_trans = dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
%     temperature=tempo_temp;
    sj_trans_out_365=tempo_trans;
else
%     temperature=cat(1,temperature,tempo_temp);
    if length(tempo_trans) == 366
        tempo_trans(60)=[];
    end
    sj_trans_out_365=cat(1,sj_trans_out_365,tempo_trans);
end
end
sj_trans_out_365(end+1:end+365) = raw(:,4);

clearvars tempo_trans sj_trans_out
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
sj_trans_out(end+1:end+365) = raw(:,4);


t_day_cut0=365*length(1989:1996);
t_day_cut1=365*length(1989:2006);
t_day_cut2=365*length(1989:2015);
sj_day_trans_1=sj_trans_out_365(t_day_cut0+1:t_day_cut1);
sj_day_trans_2=sj_trans_out_365(t_day_cut1+1:t_day_cut2);
sj_day_trans_3=sj_trans_out_365(t_day_cut2+1:end);

for i = 1:365 %days  
clim_day_trans_1(i) = nanmean(sj_day_trans_1(i:365:end));
clim_day_trans_2(i) = nanmean(sj_day_trans_2(i:365:end));
clim_day_trans_3(i) = nanmean(sj_day_trans_3(i:365:end));
end

figure; hold on;
plot(nam.climate_out_gawha_1st,'r');
plot(clim_day_trans_2,'b');

figure; hold on;
plot(nam.climate_out_gawha_2nd,'r');
plot(clim_day_trans_3,'b');

r_dis_sj_gawha_1st=corrcoef(clim_day_trans_2,nam.climate_out_gawha_1st) % correlation coefficients between illryu and gawha
r_dis_sj_gawha_2nd=corrcoef(clim_day_trans_3,nam.climate_out_gawha_2nd) % correlation coefficients between illryu and gawha


