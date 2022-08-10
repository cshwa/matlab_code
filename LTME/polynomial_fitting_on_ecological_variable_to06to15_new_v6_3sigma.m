close all; clear; clc; 
cd D:\������\Dynamic\06_river\ȯ����п�
[raw txt]=xlsread('������������_����_fix.xls','�˻����','');
% [raw txt]=xlsread('������������_��õõ.xls','�˻����','');

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

r_tp_txt=[r_txt_ud(2:end,12)];
for i = 1:length(r_tp_txt)
    if strcmp(r_tp_txt{i,1},'') == 1
       r_tp(i) = NaN;
    elseif strcmp(r_tp_txt{i,1},'') == 0
       r_tp(i) = str2num(char(r_tp_txt{i,1}));
    end
end


r_ss_txt=[r_txt_ud(2:end,15)];
for i = 1:length(r_ss_txt)
    if strcmp(r_ss_txt{i,1},'') == 1
       r_ss(i) = NaN;
    elseif strcmp(r_ss_txt{i,1},'') == 0
       r_ss(i) = str2num(char(r_ss_txt{i,1}));
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

r_po4_txt=[r_txt_ud(1:end-1,23)];
for i = 2:length(r_po4_txt)
    if strcmp(r_po4_txt{i,1},'') == 1
       r_po4(i) = NaN;
    elseif strcmp(r_po4_txt{i,1},'') == 0
       r_po4(i) = str2num(char(r_po4_txt{i,1}));
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
        raw_po4(i) = NaN;
        raw_ss(i) = NaN;
        raw_tn(i) = NaN;
        raw_tp(i) = NaN;
        r_raw_date_c{i} = '';
    else
        raw_do(i) = nanmean(r_do(indx_raw{i}));
        raw_chl(i) = nanmean(r_chl(indx_raw{i}));
        raw_no3(i) = nanmean(r_no3(indx_raw{i}));
        raw_nh4(i) = nanmean(r_nh4(indx_raw{i}));
        raw_po4(i) = nanmean(r_po4(indx_raw{i}));
        raw_ss(i) = nanmean(r_ss(indx_raw{i}));
        raw_tn(i) = nanmean(r_tn(indx_raw{i}));
        raw_tp(i) = nanmean(r_tp(indx_raw{i}));
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

%% divide the regime (in case, 1996-12, 2006-12, 2015-12)
cut_0 = find(strcmp(yymmdd_txt_c, {'1996.12.31'}) ~= 0) % 2006-12
cut_1 = find(strcmp(yymmdd_txt_c, {'2006.12.31'}) ~= 0) % 2006-12
cut_2 = find(strcmp(yymmdd_txt_c, {'2015.12.31'}) ~= 0) % 2015-12

%% cut 1989~1996, % 1997 ~ 2006, 2007 ~ 2015, 2016 ~ 2019
regime_do_1 = raw_do(cut_0+1:cut_1);
regime_do_2 = raw_do(cut_1+1:cut_2);
regime_do_3 = raw_do(cut_2+1:end);

regime_tn_1 = raw_tn(cut_0+1:cut_1);
regime_tn_2 = raw_tn(cut_1+1:cut_2);
regime_tn_3 = raw_tn(cut_2+1:end);

regime_tp_1 = raw_tp(cut_0+1:cut_1);
regime_tp_2 = raw_tp(cut_1+1:cut_2);
regime_tp_3 = raw_tp(cut_2+1:end);

regime_chl_1 = raw_chl(cut_0+1:cut_1);
regime_chl_2 = raw_chl(cut_1+1:cut_2);
regime_chl_3 = raw_chl(cut_2+1:end);

regime_ss_1 = raw_ss(cut_0+1:cut_1);
regime_ss_2 = raw_ss(cut_1+1:cut_2);
regime_ss_3 = raw_ss(cut_2+1:end);

regime_po4_1 = raw_po4(cut_0+1:cut_1);
regime_po4_2 = raw_po4(cut_1+1:cut_2);
regime_po4_3 = raw_po4(cut_2+1:end);

regime_no3_1 = raw_no3(cut_0+1:cut_1);
regime_no3_2 = raw_no3(cut_1+1:cut_2);
regime_no3_3 = raw_no3(cut_2+1:end);

regime_nh4_1 = raw_nh4(cut_0+1:cut_1);
regime_nh4_2 = raw_nh4(cut_1+1:cut_2);
regime_nh4_3 = raw_nh4(cut_2+1:end);

%% extract over 3sig
sig = 3; %% sigma
for varlist = { 'do','chl','ss','po4','no3','nh4','tn','tp'}
        clearvars varname 
        varname = char(varlist);
    for i = 1:3 % num of regime
        clearvars data data_ref
        eval(['data_ref = ','regime_',varname,'_',num2str(i),';']);
        eval(['data = ','regime_',varname,'_',num2str(i),';']);
        data(data_ref > nanmean(data_ref) + sig*nanstd(data_ref)) =NaN;
        data(data_ref < nanmean(data_ref) - sig*nanstd(data_ref)) =NaN;
        eval(['regime_',varname,'_',num2str(i),'_s = data;']);  % extracted data OUT
        disp(['regime_',varname,'_',num2str(i),'_s'])
    end
end

%% extract regime mean
% regime_do = regime_do - regm_do;
% regime_chl = regime_chl - regm_chl;
% regime_no3 = regime_no3 - regm_no3;
% regime_nh4 = regime_nh4 - regm_nh4;
% 
% regime_do_af = regime_do_af - (regm_do_af - regm_do) - regm_do;
% regime_chl_af = regime_chl_af - (regm_chl_af - regm_chl) - regm_chl;
% regime_no3_af = regime_no3_af - (regm_no3_af - regm_no3) - regm_no3;
% regime_nh4_af = regime_nh4_af - (regm_nh4_af - regm_nh4) - regm_nh4;

%% combine
com_do = [regime_do_1_s'; regime_do_2_s'; regime_do_3_s';] *0.7*44.661;
com_chl = [regime_chl_1_s'; regime_chl_2_s'; regime_chl_3_s';];
com_no3 = [regime_no3_1_s'; regime_no3_2_s'; regime_no3_3_s';] .*1000 ./14.006720;
com_nh4 = [regime_nh4_1_s'; regime_nh4_2_s'; regime_nh4_3_s';] .*1000 ./14.006720;
com_po4 = [regime_po4_1_s'; regime_po4_2_s'; regime_po4_3_s';] .*1000 ./30.973762;
com_ss = [regime_ss_1_s'; regime_ss_2_s'; regime_ss_3_s';];
com_tn = [regime_tn_1_s'; regime_tn_2_s'; regime_tn_3_s';];
com_tp = [regime_tp_1_s'; regime_tp_2_s'; regime_tp_3_s';];

return

clearvars regime_*
regime_do = com_do;
regime_chl = com_chl;
regime_no3 = com_no3;
regime_nh4 = com_nh4;
regime_po4 = com_po4;
regime_ss = com_ss;
regime_tn = com_tn;
regime_tp = com_tp;

%% monthly mean (yy.mm) form after over 3sigma value extraction.

cut_indx= find(t_indx ==  cut_0); % from 1997-01-01
t_indx_97 = t_indx - t_indx(cut_indx); t_indx_97(t_indx_97 <= 0 ) = [];
for i=1:length(t_indx_97) % time for num. of months 
if i ==1
   monthly_do(i) = nanmean(regime_do(1:t_indx_97(i)));
   monthly_chl(i) = nanmean(regime_chl(1:t_indx_97(i)));
   monthly_no3(i) = nanmean(regime_no3(1:t_indx_97(i)));
   monthly_nh4(i) = nanmean(regime_nh4(1:t_indx_97(i)));
   monthly_po4(i) = nanmean(regime_po4(1:t_indx_97(i)));
   monthly_ss(i) = nanmean(regime_ss(1:t_indx_97(i)));
   monthly_tn(i) = nanmean(regime_tn(1:t_indx_97(i)));
   monthly_tp(i) = nanmean(regime_tp(1:t_indx_97(i)));
else
   monthly_do(i) = nanmean(regime_do(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_chl(i) = nanmean(regime_chl(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_no3(i) = nanmean(regime_no3(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_nh4(i) = nanmean(regime_nh4(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_po4(i) = nanmean(regime_po4(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_ss(i) = nanmean(regime_ss(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_tn(i) = nanmean(regime_tn(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_tp(i) = nanmean(regime_tp(t_indx_97(i-1)+1:t_indx_97(i)));
end
end

% fill missing value
mon_in_chl = monthly_chl;
mon_in_no3 = monthly_no3;
mon_in_nh4 = monthly_nh4;
mon_in_do = monthly_do;
mon_in_po4 = monthly_po4;
mon_in_ss = monthly_ss;

t=1:length(t_indx_97);
mon_in_chl(isnan(mon_in_chl))=interp1(t(~isnan(mon_in_chl)),mon_in_chl(~isnan(mon_in_chl)),t(isnan(mon_in_chl)));
mon_in_no3(isnan(mon_in_no3))=interp1(t(~isnan(mon_in_no3)),mon_in_no3(~isnan(mon_in_no3)),t(isnan(mon_in_no3)));
mon_in_nh4(isnan(mon_in_nh4))=interp1(t(~isnan(mon_in_nh4)),mon_in_nh4(~isnan(mon_in_nh4)),t(isnan(mon_in_nh4)));
mon_in_do(isnan(mon_in_do))=interp1(t(~isnan(mon_in_do)),mon_in_do(~isnan(mon_in_do)),t(isnan(mon_in_do)));
mon_in_po4(isnan(mon_in_po4))=interp1(t(~isnan(mon_in_po4)),mon_in_po4(~isnan(mon_in_po4)),t(isnan(mon_in_po4)));
mon_in_ss(isnan(mon_in_ss))=interp1(t(~isnan(mon_in_ss)),mon_in_ss(~isnan(mon_in_ss)),t(isnan(mon_in_ss)));

save('songjung_yymm_monthly_data_to06to15_3sig.mat','mon*','t_indx*','yymmdd_txt_c');
save('songjung_yymm_raw_data_to06to15_3sig.mat','com*','t_indx*','yymmdd_txt_c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NAMGANG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\������\Dynamic\06_river\ȯ����п�
% [raw txt]=xlsread('������������_����_fix.xls','�˻����','');
[raw txt]=xlsread('������������_������1.xls','�˻����','');

dash_c = '.';
% r_txt_ud = flipud(txt);
r_txt_ud = txt;
r_date_txt=[char(r_txt_ud(2:end,5))];

r_ss_txt=[r_txt_ud(2:end,15)];
for i = 1:length(r_ss_txt)
    if strcmp(r_ss_txt{i,1},'') == 1
       r_ss(i) = NaN;
    elseif strcmp(r_ss_txt{i,1},'') == 0
       r_ss(i) = str2num(char(r_ss_txt{i,1}));
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

r_po4_txt=[r_txt_ud(1:end-1,23)];
for i = 2:length(r_po4_txt)
    if strcmp(r_po4_txt{i,1},'') == 1
       r_po4(i) = NaN;
    elseif strcmp(r_po4_txt{i,1},'') == 0
       r_po4(i) = str2num(char(r_po4_txt{i,1}));
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
        raw_po4(i) = NaN; 
        raw_ss(i) = NaN; 
        r_raw_date_c{i} = '';
    else
        raw_do(i) = nanmean(r_do(indx_raw{i}));
        raw_chl(i) = nanmean(r_chl(indx_raw{i}));
        raw_no3(i) = nanmean(r_no3(indx_raw{i}));
        raw_nh4(i) = nanmean(r_nh4(indx_raw{i}));
        raw_po4(i) = nanmean(r_po4(indx_raw{i}));
        raw_ss(i) = nanmean(r_ss(indx_raw{i}));
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
%% divide the regime (in case, 1996-12, 2006-12, 2015-12)
cut_0 = find(strcmp(yymmdd_txt_c, {'1996.12.31'}) ~= 0) % 2006-12
cut_1 = find(strcmp(yymmdd_txt_c, {'2006.12.31'}) ~= 0) % 2006-12
cut_2 = find(strcmp(yymmdd_txt_c, {'2015.12.31'}) ~= 0) % 2015-12

%% cut 1989~1996, % 1997 ~ 2006, 2007 ~ 2015, 2016 ~ 2019
regime_do_1 = raw_do(cut_0+1:cut_1);
regime_do_2 = raw_do(cut_1+1:cut_2);
regime_do_3 = raw_do(cut_2+1:end);

regime_chl_1 = raw_chl(cut_0+1:cut_1);
regime_chl_2 = raw_chl(cut_1+1:cut_2);
regime_chl_3 = raw_chl(cut_2+1:end);

regime_ss_1 = raw_ss(cut_0+1:cut_1);
regime_ss_2 = raw_ss(cut_1+1:cut_2);
regime_ss_3 = raw_ss(cut_2+1:end);

regime_po4_1 = raw_po4(cut_0+1:cut_1);
regime_po4_2 = raw_po4(cut_1+1:cut_2);
regime_po4_3 = raw_po4(cut_2+1:end);

regime_no3_1 = raw_no3(cut_0+1:cut_1);
regime_no3_2 = raw_no3(cut_1+1:cut_2);
regime_no3_3 = raw_no3(cut_2+1:end);

regime_nh4_1 = raw_nh4(cut_0+1:cut_1);
regime_nh4_2 = raw_nh4(cut_1+1:cut_2);
regime_nh4_3 = raw_nh4(cut_2+1:end);

%% extract over 3sig
sig = 3; %% sigma
for varlist = { 'do','chl','ss','po4','no3','nh4'}
        clearvars varname 
        varname = char(varlist);
    for i = 1:3 % num of regime
        clearvars data data_ref
        eval(['data_ref = ','regime_',varname,'_',num2str(i),';']);
        eval(['data = ','regime_',varname,'_',num2str(i),';']);
        data(data_ref > nanmean(data_ref) + sig*nanstd(data_ref)) =NaN;
        data(data_ref < nanmean(data_ref) - sig*nanstd(data_ref)) =NaN;
        eval(['regime_',varname,'_',num2str(i),'_s = data;']);  % extracted data OUT
        disp(['regime_',varname,'_',num2str(i),'_s'])
    end
end

%% extract regime mean
% regime_do = regime_do - regm_do;
% regime_chl = regime_chl - regm_chl;
% regime_no3 = regime_no3 - regm_no3;
% regime_nh4 = regime_nh4 - regm_nh4;
% 
% regime_do_af = regime_do_af - (regm_do_af - regm_do) - regm_do;
% regime_chl_af = regime_chl_af - (regm_chl_af - regm_chl) - regm_chl;
% regime_no3_af = regime_no3_af - (regm_no3_af - regm_no3) - regm_no3;
% regime_nh4_af = regime_nh4_af - (regm_nh4_af - regm_nh4) - regm_nh4;

%% combine
com_do = [regime_do_1_s'; regime_do_2_s'; regime_do_3_s';] *0.7*44.661;
com_chl = [regime_chl_1_s'; regime_chl_2_s'; regime_chl_3_s';];
com_no3 = [regime_no3_1_s'; regime_no3_2_s'; regime_no3_3_s';] .*1000 ./14.006720;
com_nh4 = [regime_nh4_1_s'; regime_nh4_2_s'; regime_nh4_3_s';] .*1000 ./14.006720;
com_po4 = [regime_po4_1_s'; regime_po4_2_s'; regime_po4_3_s';] .*1000 ./30.973762;
com_ss = [regime_ss_1_s'; regime_ss_2_s'; regime_ss_3_s';];

return

clearvars regime_*
regime_do = com_do;
regime_chl = com_chl;
regime_no3 = com_no3;
regime_nh4 = com_nh4;
regime_po4 = com_po4;
regime_ss = com_ss;


%% monthly mean (yy.mm) form after over 3sigma value extraction.

cut_indx= find(t_indx ==  cut_0); % from 1997-01-01
t_indx_97 = t_indx - t_indx(cut_indx); t_indx_97(t_indx_97 <= 0 ) = [];
for i=1:length(t_indx_97) % time for num. of months 
if i ==1
   monthly_do(i) = nanmean(regime_do(1:t_indx_97(i)));
   monthly_chl(i) = nanmean(regime_chl(1:t_indx_97(i)));
   monthly_no3(i) = nanmean(regime_no3(1:t_indx_97(i)));
   monthly_nh4(i) = nanmean(regime_nh4(1:t_indx_97(i)));
   monthly_po4(i) = nanmean(regime_po4(1:t_indx_97(i)));
   monthly_ss(i) = nanmean(regime_ss(1:t_indx_97(i)));
else
   monthly_do(i) = nanmean(regime_do(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_chl(i) = nanmean(regime_chl(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_no3(i) = nanmean(regime_no3(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_nh4(i) = nanmean(regime_nh4(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_po4(i) = nanmean(regime_po4(t_indx_97(i-1)+1:t_indx_97(i)));
   monthly_ss(i) = nanmean(regime_ss(t_indx_97(i-1)+1:t_indx_97(i)));
end
end

% fill missing value
mon_in_chl = monthly_chl;
mon_in_no3 = monthly_no3;
mon_in_nh4 = monthly_nh4;
mon_in_do = monthly_do;
mon_in_po4 = monthly_po4;
mon_in_ss = monthly_ss;

t=1:length(t_indx_97);
mon_in_chl(isnan(mon_in_chl))=interp1(t(~isnan(mon_in_chl)),mon_in_chl(~isnan(mon_in_chl)),t(isnan(mon_in_chl)));
mon_in_no3(isnan(mon_in_no3))=interp1(t(~isnan(mon_in_no3)),mon_in_no3(~isnan(mon_in_no3)),t(isnan(mon_in_no3)));
mon_in_nh4(isnan(mon_in_nh4))=interp1(t(~isnan(mon_in_nh4)),mon_in_nh4(~isnan(mon_in_nh4)),t(isnan(mon_in_nh4)));
mon_in_do(isnan(mon_in_do))=interp1(t(~isnan(mon_in_do)),mon_in_do(~isnan(mon_in_do)),t(isnan(mon_in_do)));
mon_in_po4(isnan(mon_in_po4))=interp1(t(~isnan(mon_in_po4)),mon_in_po4(~isnan(mon_in_po4)),t(isnan(mon_in_po4)));
mon_in_ss(isnan(mon_in_ss))=interp1(t(~isnan(mon_in_ss)),mon_in_ss(~isnan(mon_in_ss)),t(isnan(mon_in_ss)));

save('namgang_yymm_monthly_data_to06to15_3sig.mat','mon*','t_indx*','yymmdd_txt_c');
save('namgang_yymm_raw_data_to06to15_3sig.mat','com*','t_indx*','yymmdd_txt_c');