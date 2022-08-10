close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\환경과학원
% [raw txt]=xlsread('수질측정지점_구례_fix.xls','검색결과','');
[raw txt]=xlsread('수질측정지점_구례_1989-2018.xls','검색결과','');
[raw1 txt1]=xlsread('수질_일반측정망_구례_201901-201912.xls','수질(일반측정망)','');
[raw2 txt2]=xlsread('수질_일반측정망_구례_202001-202012.xls','수질(일반측정망)','');
[raw3 txt3]=xlsread('수질자료_구례(일)_2021.xlsx','수질측정망 일자료 조회','');

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
for i = 1989:2021
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
for i = 1:length(1989:2021)
    l=0
    for n = 1:12
        m = m+1;
    for j = 1:eom_d_raw(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        yymmdd_txt(k,:)=[num2str(i+1988) '.' num2str(n,'%02d') '.'  num2str(j,'%02d')];
        yymmdd_txt_slash(k,:)=[num2str(i+1988) '/' num2str(n,'%02d') '/'  num2str(j,'%02d')];
    end
    end
end


% make it to cell-array
for i = 1:length(yymmdd_txt)
    yymmdd_txt_c{i,1} = yymmdd_txt(i,:); % raw
    yymmdd_txt_slash_c{i,1} = yymmdd_txt_slash(i,:); % raw
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

%% 2019 (10958) ~ 2021
% pick matched date from water temp date
% r_date_txt=[char(r_txt_ud(2:end,5))];
r_date_txt_to21 = [char(txt1(3:end,3))];
r_date_txt_to21(end+1:end+length(txt2(3:end,3)),:)=[char(txt2(3:end,3))];
r_date_txt_to21(end+1:end+length(txt3(3:end,2)),:)=[char(txt3(3:end,2))];

for i = 1:length(r_date_txt_to21)
    r_raw_txt_c_to21{i,1} = r_date_txt_to21(i,:); % raw
end

% obs value
% %  DO : col 7
% %  Chl : col 35
% %  no3 : col 32
% %  nh4 : col 31
% %  po4 : col 34
% %  ss : col 10
% %  tn : col 11
% %  tp : col 12

do_19to21 = [raw1(:,7); raw2(:,7); raw3(:,2);];
chl_19to21 = [raw1(:,35); raw2(:,35); raw3(:,5);];
no3_19to21 = [raw1(:,32); raw2(:,32); raw3(:,12);];
nh4_19to21 = [raw1(:,31); raw2(:,31); raw3(:,11);];
po4_19to21 = [raw1(:,34); raw2(:,34); raw3(:,14);];
ss_19to21 = [raw1(:,10); raw2(:,10); raw3(:,15);];
tn_19to21 = [raw1(:,11); raw2(:,11); raw3(:,6);];
tp_19to21 = [raw1(:,12); raw2(:,12); raw3(:,7);];

for i = 10958 :length(yymmdd_txt_slash_c)
       indx_raw_19to21{i} = find([strcmp(yymmdd_txt_slash_c{i}, r_raw_txt_c_to21)] == 1);
end
% yymmdd_txt_slash_c{11679}

for i = 10958:length(yymmdd_txt_slash_c)
    if size(indx_raw_19to21{i},1) == 0
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
        raw_do(i) = nanmean(do_19to21(indx_raw_19to21{i}));
        raw_chl(i) = nanmean(chl_19to21(indx_raw_19to21{i}));
        raw_no3(i) = nanmean(no3_19to21(indx_raw_19to21{i}));
        raw_nh4(i) = nanmean(nh4_19to21(indx_raw_19to21{i}));
        raw_po4(i) = nanmean(po4_19to21(indx_raw_19to21{i}));
        raw_ss(i) = nanmean(ss_19to21(indx_raw_19to21{i}));
        raw_tn(i) = nanmean(tn_19to21(indx_raw_19to21{i}));
        raw_tp(i) = nanmean(tp_19to21(indx_raw_19to21{i}));
        r_raw_date_c{i} = r_date_txt_to21(indx_raw_19to21{i},6:end);  %remove year
    end
end

return

% make 1989~present
k=0
for i = 1989:2021
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

% %% divide the regime (in case, 1996-12, 2006-12, 2015-12)
% cut_0 = find(strcmp(yymmdd_txt_c, {'1992.12.31'}) ~= 0) % 2006-12
% cut_1 = find(strcmp(yymmdd_txt_c, {'2004.12.31'}) ~= 0) % 2006-12
% cut_2 = find(strcmp(yymmdd_txt_c, {'2006.12.31'}) ~= 0) % 2015-12
% cut_3 = find(strcmp(yymmdd_txt_c, {'2015.12.31'}) ~= 0) % 2015-12
% 
% %% cut 1989~1996, % 1997 ~ 2006, 2007 ~ 2015, 2016 ~ 2020
% regime_do_1 = raw_do(cut_0+1:cut_1);
% regime_do_2 = raw_do(cut_2+1:cut_3);                                                              
% regime_do_3 = raw_do(cut_3+1:end);


%% extract over 3sig
clearvars regime_*
sig = 3; %% sigma
for varlist = { 'do','chl','ss','po4','no3','nh4','tn','tp'}
        clearvars varname 
        varname = char(varlist);
%     for i = 1:1 % num of regime
        clearvars data data_ref
        eval(['data_ref = ','raw_',varname,';']);
        eval(['data = ','raw_',varname,';']);
        data(data_ref > nanmean(data_ref) + sig*nanstd(data_ref)) =NaN;
        data(data_ref < nanmean(data_ref) - sig*nanstd(data_ref)) =NaN;
        eval([varname,'_s = data;']);  % extracted data OUT
        disp([varname,'_s'])
%     end
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
com_do = do_s *0.7*44.661;
com_chl = chl_s;
com_no3 = no3_s .*1000 ./14.006720;
com_nh4 = nh4_s .*1000 ./14.006720;
com_po4 = po4_s .*1000 ./30.973762;
com_ss = ss_s;
com_tn = tn_s .*1000 ./14.006720;
com_tp = tp_s .*1000 ./30.973762;

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

% cut_indx= find(t_indx ==  cut_0 ); % from 1993-01-01
% t_indx_93 = t_indx - t_indx(cut_indx); t_indx_93(t_indx_93 <= 0 ) = [];
% t_indx_93(t_indx_93 > (cut_1 - cut_0) ) = []; % 1993 ~ 2004

for i=1:length(t_indx) % time for num. of months 
if i ==1
   monthly_do(i) = nanmean(do_s(1:t_indx(i)));
   monthly_chl(i) = nanmean(chl_s(1:t_indx(i)));
   monthly_no3(i) = nanmean(no3_s(1:t_indx(i)));
   monthly_nh4(i) = nanmean(nh4_s(1:t_indx(i)));
   monthly_po4(i) = nanmean(po4_s(1:t_indx(i)));
   monthly_ss(i) = nanmean(ss_s(1:t_indx(i)));
   monthly_tn(i) = nanmean(tn_s(1:t_indx(i)));
   monthly_tp(i) = nanmean(tp_s(1:t_indx(i)));
else
   monthly_do(i) = nanmean(do_s(t_indx(i-1)+1:t_indx(i)));
   monthly_chl(i) = nanmean(chl_s(t_indx(i-1)+1:t_indx(i)));
   monthly_no3(i) = nanmean(no3_s(t_indx(i-1)+1:t_indx(i)));
   monthly_nh4(i) = nanmean(nh4_s(t_indx(i-1)+1:t_indx(i)));
   monthly_po4(i) = nanmean(po4_s(t_indx(i-1)+1:t_indx(i)));
   monthly_ss(i) = nanmean(ss_s(t_indx(i-1)+1:t_indx(i)));
   monthly_tn(i) = nanmean(tn_s(t_indx(i-1)+1:t_indx(i)));
   monthly_tp(i) = nanmean(tp_s(t_indx(i-1)+1:t_indx(i)));
end
end

% fill missing value
mon_in_chl = monthly_chl;
mon_in_no3 = monthly_no3;
mon_in_nh4 = monthly_nh4;
mon_in_do = monthly_do;
mon_in_po4 = monthly_po4;
mon_in_ss = monthly_ss;
mon_in_tn = monthly_tn;
mon_in_tp = monthly_tp;

t=1:length(t_indx);
mon_in_chl(isnan(mon_in_chl))=interp1(t(~isnan(mon_in_chl)),mon_in_chl(~isnan(mon_in_chl)),t(isnan(mon_in_chl)));
mon_in_no3(isnan(mon_in_no3))=interp1(t(~isnan(mon_in_no3)),mon_in_no3(~isnan(mon_in_no3)),t(isnan(mon_in_no3)));
mon_in_nh4(isnan(mon_in_nh4))=interp1(t(~isnan(mon_in_nh4)),mon_in_nh4(~isnan(mon_in_nh4)),t(isnan(mon_in_nh4)));
mon_in_do(isnan(mon_in_do))=interp1(t(~isnan(mon_in_do)),mon_in_do(~isnan(mon_in_do)),t(isnan(mon_in_do)));
mon_in_po4(isnan(mon_in_po4))=interp1(t(~isnan(mon_in_po4)),mon_in_po4(~isnan(mon_in_po4)),t(isnan(mon_in_po4)));
mon_in_ss(isnan(mon_in_ss))=interp1(t(~isnan(mon_in_ss)),mon_in_ss(~isnan(mon_in_ss)),t(isnan(mon_in_ss)));
mon_in_tn(isnan(mon_in_tn))=interp1(t(~isnan(mon_in_tn)),mon_in_tn(~isnan(mon_in_tn)),t(isnan(mon_in_tn)));
mon_in_tp(isnan(mon_in_tp))=interp1(t(~isnan(mon_in_tp)),mon_in_tp(~isnan(mon_in_tp)),t(isnan(mon_in_tp)));

save('songjung_yymm_monthly_data_3sig_1989to2021.mat','mon*','t_indx*','yymmdd_txt_c');
save('songjung_yymm_raw_data_3sig_1989to2021.mat','com*','t_indx*','yymmdd_txt_c');
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