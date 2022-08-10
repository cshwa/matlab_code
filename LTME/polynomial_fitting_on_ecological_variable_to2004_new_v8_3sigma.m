close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\환경과학원
[raw txt]=xlsread('수질측정지점_구례_1989-2018.xls','검색결과','');
[raw2 txt2]=xlsread('수질_일반측정망_구례_201901-201912.xls','수질(일반측정망)','');

dash_c = '.';
% r_txt_ud = flipud(txt);
r_txt_ud = txt;
r_date_txt=[char(r_txt_ud(2:end,5))];

dash_c = '.';
% r_txt_ud = flipud(txt);
r_txt_ud2 = txt2;
r_date_txt2=[char(r_txt_ud2(3:end,3))];


r_ss_txt=[r_txt_ud(2:end,15)];
for i = 1:length(r_ss_txt)
    if strcmp(r_ss_txt{i,1},'') == 1
       r_ss(i) = NaN;
    elseif strcmp(r_ss_txt{i,1},'') == 0
       r_ss(i) = str2num(char(r_ss_txt{i,1}));
    end
end

r_ss_raw2=[raw2(1:end,10)];
r_ss(end+1:end+length(r_ss_raw2)) = r_ss_raw2;


r_do_txt=[r_txt_ud(2:end,9)];
for i = 1:length(r_do_txt)
    if strcmp(r_do_txt{i,1},'') == 1
       r_do(i) = NaN;
    elseif strcmp(r_do_txt{i,1},'') == 0
       r_do(i) = str2num(char(r_do_txt{i,1}));
    end
end

r_do_raw2=[raw2(1:end,7)];
r_do(end+1:end+length(r_do_raw2)) = r_do_raw2;


r_chl_txt=[r_txt_ud(2:end,16)];
for i = 1:length(r_chl_txt)
    if strcmp(r_chl_txt{i,1},'') == 1
       r_chl(i) = NaN;
    elseif strcmp(r_chl_txt{i,1},'') == 0
       r_chl(i) = str2num(char(r_chl_txt{i,1}));
    end
end

r_chl_raw2=[raw2(1:end,35)];
r_chl(end+1:end+length(r_chl_raw2)) = r_chl_raw2;


r_no3_txt=[r_txt_ud(2:end,20)];
for i = 1:length(r_no3_txt)
    if strcmp(r_no3_txt{i,1},'') == 1
       r_no3(i) = NaN;
    elseif strcmp(r_no3_txt{i,1},'') == 0
       r_no3(i) = str2num(char(r_no3_txt{i,1}));
    end
end

r_no3_raw2=[raw2(1:end,32)];
r_no3(end+1:end+length(r_no3_raw2)) = r_no3_raw2;


r_nh4_txt=[r_txt_ud(2:end,21)];
for i = 1:length(r_nh4_txt)
    if strcmp(r_nh4_txt{i,1},'') == 1
       r_nh4(i) = NaN;
    elseif strcmp(r_nh4_txt{i,1},'') == 0
       r_nh4(i) = str2num(char(r_nh4_txt{i,1}));
    end
end

r_nh4_raw2=[raw2(1:end,31)];
r_nh4(end+1:end+length(r_nh4_raw2)) = r_nh4_raw2;

r_po4_txt=[r_txt_ud(1:end-1,23)];
for i = 2:length(r_po4_txt)
    if strcmp(r_po4_txt{i,1},'') == 1
       r_po4(i) = NaN;
    elseif strcmp(r_po4_txt{i,1},'') == 0
       r_po4(i) = str2num(char(r_po4_txt{i,1}));
    end
end

r_po4_raw2=[raw2(1:end,34)];
r_po4(end+1:end+length(r_po4_raw2)) = r_po4_raw2;


r_toc_txt=[r_txt_ud(1:end-1,17)];
for i = 2:length(r_toc_txt)
    if strcmp(r_toc_txt{i,1},'') == 1
       r_toc(i) = NaN;
    elseif strcmp(r_toc_txt{i,1},'') == 0
       r_toc(i) = str2num(char(r_toc_txt{i,1}));
    end
end

r_toc_raw2=[raw2(1:end,13)];
r_toc(end+1:end+length(r_toc_raw2)) = r_toc_raw2;

for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
    r_raw_txt_c{i,1} = r_date_txt(i,:); % raw
end

k=0;
for i= length(r_date_txt)+1:length(r_date_txt)+length(r_date_txt2)
    k=k+1;
    r_date_txt_c{i,1} = [r_date_txt2(k,6:7),'.',r_date_txt2(k,9:10)]; % delete year
    r_raw_txt_c{i,1} = [r_date_txt2(k,1:4),'.',r_date_txt2(k,6:7),'.',r_date_txt2(k,9:10)]; % raw
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
        raw_toc(i) = NaN;
        raw_ss(i) = NaN;
        r_raw_date_c{i} = '';
    else
        raw_do(i) = nanmean(r_do(indx_raw{i}));
        raw_chl(i) = nanmean(r_chl(indx_raw{i}));
        raw_no3(i) = nanmean(r_no3(indx_raw{i}));
        raw_nh4(i) = nanmean(r_nh4(indx_raw{i}));
        raw_po4(i) = nanmean(r_po4(indx_raw{i}));
        raw_ss(i) = nanmean(r_ss(indx_raw{i}));
        raw_toc(i) = nanmean(r_toc(indx_raw{i}));
        r_raw_date_c{i} = r_raw_txt_c{indx_raw{i}}(6:10);  %remove year
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
regime_po4_pre = raw_po4(1:t_tick(end));
regime_ss_pre = raw_ss(1:t_tick(end));
regime_toc_pre = raw_toc(1:t_tick(end));
% 2004 ~ 2019
regime_do_af_pre = raw_do(t_tick(end)+1:end);
regime_chl_af_pre = raw_chl(t_tick(end)+1:end);
regime_no3_af_pre = raw_no3(t_tick(end)+1:end);
regime_nh4_af_pre = raw_nh4(t_tick(end)+1:end);
regime_po4_af_pre = raw_po4(t_tick(end)+1:end);
regime_ss_af_pre = raw_ss(t_tick(end)+1:end);
regime_toc_af_pre = raw_toc(t_tick(end)+1:end);

%% extract over 3sig
sig = 3; %% sigma
% ~2004
clearvars idx1
idx1 = find(isnan(regime_ss_pre) == 0);
regime_ss=regime_ss_pre;
regime_ss(find(regime_ss_pre > mean(regime_ss_pre(idx1)) + sig*std(regime_ss_pre(idx1))))=NaN;
regime_ss(find(regime_ss_pre < mean(regime_ss_pre(idx1)) - sig*std(regime_ss_pre(idx1))))=NaN;
regm_ss = nanmean(regime_ss)

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

clearvars idx1
idx1 = find(isnan(regime_po4_pre) == 0);
regime_po4=regime_po4_pre;
regime_po4(find(regime_po4_pre > mean(regime_po4_pre(idx1)) + sig*std(regime_po4_pre(idx1))))=NaN;
regime_po4(find(regime_po4_pre < mean(regime_po4_pre(idx1)) - sig*std(regime_po4_pre(idx1))))=NaN;
regm_po4 = nanmean(regime_po4)

clearvars idx1
idx1 = find(isnan(regime_toc_pre) == 0);
regime_toc=regime_toc_pre;
regime_toc(find(regime_toc_pre > mean(regime_toc_pre(idx1)) + sig*std(regime_toc_pre(idx1))))=NaN;
regime_toc(find(regime_toc_pre < mean(regime_toc_pre(idx1)) - sig*std(regime_toc_pre(idx1))))=NaN;
regm_toc = nanmean(regime_toc)

% 2004~
clearvars idx1
idx1 = find(isnan(regime_ss_af_pre) == 0);
regime_ss_af=regime_ss_af_pre;
regime_ss_af(find(regime_ss_af_pre > mean(regime_ss_af_pre(idx1)) + sig*std(regime_ss_af_pre(idx1))))=NaN;
regime_ss_af(find(regime_ss_af_pre < mean(regime_ss_af_pre(idx1)) - sig*std(regime_ss_af_pre(idx1))))=NaN;
regm_ss_af = nanmean(regime_ss_af)

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

clearvars idx1
idx1 = find(isnan(regime_po4_af_pre) == 0);
regime_po4_af=regime_po4_af_pre;
regime_po4_af(find(regime_po4_af_pre > mean(regime_po4_af_pre(idx1)) + sig*std(regime_po4_af_pre(idx1))))=NaN;
regime_po4_af(find(regime_po4_af_pre < mean(regime_po4_af_pre(idx1)) - sig*std(regime_po4_af_pre(idx1))))=NaN;
regm_po4_af = nanmean(regime_po4_af)

clearvars idx1
idx1 = find(isnan(regime_toc_af_pre) == 0);
regime_toc_af=regime_toc_af_pre;
regime_toc_af(find(regime_toc_af_pre > mean(regime_toc_af_pre(idx1)) + sig*std(regime_toc_af_pre(idx1))))=NaN;
regime_toc_af(find(regime_toc_af_pre < mean(regime_toc_af_pre(idx1)) - sig*std(regime_toc_af_pre(idx1))))=NaN;
regm_toc_af = nanmean(regime_toc_af)

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
com_do = [regime_do'; regime_do_af';] *0.7*44.661;
com_chl = [regime_chl'; regime_chl_af';];
com_no3 = [regime_no3'; regime_no3_af';] .*1000 ./14.006720;
com_nh4 = [regime_nh4'; regime_nh4_af';] .*1000 ./14.006720;
com_po4 = [regime_po4'; regime_po4_af';] .*1000 ./30.973762;
com_toc = [regime_toc'; regime_toc_af';]; % mg/L
com_ss = [regime_ss'; regime_ss_af';]; 



clearvars regime_*
regime_do = com_do;
regime_chl = com_chl;
regime_no3 = com_no3;
regime_nh4 = com_nh4;
regime_po4 = com_po4;
regime_toc = com_toc;
regime_ss = com_ss;


%% monthly mean (yy.mm) form after over 3sigma value extraction.

for i=1:length(t_indx) % time for num. of months 
if i ==1
   monthly_do(i) = nanmean(regime_do(1:t_indx(i)));
   monthly_chl(i) = nanmean(regime_chl(1:t_indx(i)));
   monthly_no3(i) = nanmean(regime_no3(1:t_indx(i)));
   monthly_nh4(i) = nanmean(regime_nh4(1:t_indx(i)));
   monthly_po4(i) = nanmean(regime_po4(1:t_indx(i)));
   monthly_toc(i) = nanmean(regime_toc(1:t_indx(i)));
   monthly_ss(i) = nanmean(regime_ss(1:t_indx(i)));
else
   monthly_do(i) = nanmean(regime_do(t_indx(i-1)+1:t_indx(i)));
   monthly_chl(i) = nanmean(regime_chl(t_indx(i-1)+1:t_indx(i)));
   monthly_no3(i) = nanmean(regime_no3(t_indx(i-1)+1:t_indx(i)));
   monthly_nh4(i) = nanmean(regime_nh4(t_indx(i-1)+1:t_indx(i)));
   monthly_po4(i) = nanmean(regime_po4(t_indx(i-1)+1:t_indx(i)));
   monthly_toc(i) = nanmean(regime_toc(t_indx(i-1)+1:t_indx(i)));
   monthly_ss(i) = nanmean(regime_ss(t_indx(i-1)+1:t_indx(i)));
end
end

% fill missing value
mon_in_chl = monthly_chl;
mon_in_no3 = monthly_no3;
mon_in_nh4 = monthly_nh4;
mon_in_do = monthly_do;
mon_in_po4 = monthly_po4;
mon_in_toc = monthly_toc;
mon_in_ss = monthly_ss;

t=1:length(t_indx);
mon_in_chl(isnan(mon_in_chl))=interp1(t(~isnan(mon_in_chl)),mon_in_chl(~isnan(mon_in_chl)),t(isnan(mon_in_chl)));
mon_in_no3(isnan(mon_in_no3))=interp1(t(~isnan(mon_in_no3)),mon_in_no3(~isnan(mon_in_no3)),t(isnan(mon_in_no3)));
mon_in_nh4(isnan(mon_in_nh4))=interp1(t(~isnan(mon_in_nh4)),mon_in_nh4(~isnan(mon_in_nh4)),t(isnan(mon_in_nh4)));
mon_in_do(isnan(mon_in_do))=interp1(t(~isnan(mon_in_do)),mon_in_do(~isnan(mon_in_do)),t(isnan(mon_in_do)));
mon_in_po4(isnan(mon_in_po4))=interp1(t(~isnan(mon_in_po4)),mon_in_po4(~isnan(mon_in_po4)),t(isnan(mon_in_po4)));
mon_in_toc(isnan(mon_in_toc))=interp1(t(~isnan(mon_in_toc)),mon_in_toc(~isnan(mon_in_toc)),t(isnan(mon_in_toc)));
mon_in_ss(isnan(mon_in_ss))=interp1(t(~isnan(mon_in_ss)),mon_in_ss(~isnan(mon_in_ss)),t(isnan(mon_in_ss)));
return

save('songjung_yymm_monthly_data_89to19_3sig.mat','mon*','t_indx','yymmdd_txt_c');
save('songjung_yymm_raw_data_89to19_3sig.mat','com*','t_indx','yymmdd_txt_c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NAMGANG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\환경과학원
% [raw txt]=xlsread('수질측정지점_구례_fix.xls','검색결과','');
[raw txt]=xlsread('수질측정지점_남강댐1.xls','검색결과','');
[raw2 txt2]=xlsread('남강댐1_수질_일반측정망_201901-201912.xls','수질(일반측정망)','');

dash_c = '.';
% r_txt_ud = flipud(txt);
r_txt_ud = txt;
r_date_txt=[char(r_txt_ud(2:end,5))];

dash_c = '.';
% r_txt_ud = flipud(txt);
r_txt_ud2 = txt2;
r_date_txt2=[char(r_txt_ud2(3:end,3))];

r_ss_txt=[r_txt_ud(2:end,15)];
for i = 1:length(r_ss_txt)
    if strcmp(r_ss_txt{i,1},'') == 1
       r_ss(i) = NaN;
    elseif strcmp(r_ss_txt{i,1},'') == 0
       r_ss(i) = str2num(char(r_ss_txt{i,1}));
    end
end

r_ss_raw2=[raw2(1:end,10)];
r_ss(end+1:end+length(r_ss_raw2)) = r_ss_raw2;

r_do_txt=[r_txt_ud(2:end,9)];
for i = 1:length(r_do_txt)
    if strcmp(r_do_txt{i,1},'') == 1
       r_do(i) = NaN;
    elseif strcmp(r_do_txt{i,1},'') == 0
       r_do(i) = str2num(char(r_do_txt{i,1}));
    end
end

r_do_raw2=[raw2(1:end,7)];
r_do(end+1:end+length(r_do_raw2)) = r_do_raw2;

r_chl_txt=[r_txt_ud(2:end,16)];
for i = 1:length(r_chl_txt)
    if strcmp(r_chl_txt{i,1},'') == 1
       r_chl(i) = NaN;
    elseif strcmp(r_chl_txt{i,1},'') == 0
       r_chl(i) = str2num(char(r_chl_txt{i,1}));
    end
end

r_chl_raw2=[raw2(1:end,35)];
r_chl(end+1:end+length(r_chl_raw2)) = r_chl_raw2;

r_no3_txt=[r_txt_ud(2:end,20)];
for i = 1:length(r_no3_txt)
    if strcmp(r_no3_txt{i,1},'') == 1
       r_no3(i) = NaN;
    elseif strcmp(r_no3_txt{i,1},'') == 0
       r_no3(i) = str2num(char(r_no3_txt{i,1}));
    end
end

r_no3_raw2=[raw2(1:end,32)];
r_no3(end+1:end+length(r_no3_raw2)) = r_no3_raw2;

r_nh4_txt=[r_txt_ud(2:end,21)];
for i = 1:length(r_nh4_txt)
    if strcmp(r_nh4_txt{i,1},'') == 1
       r_nh4(i) = NaN;
    elseif strcmp(r_nh4_txt{i,1},'') == 0
       r_nh4(i) = str2num(char(r_nh4_txt{i,1}));
    end
end

r_nh4_raw2=[raw2(1:end,31)];
r_nh4(end+1:end+length(r_nh4_raw2)) = r_nh4_raw2;


r_po4_txt=[r_txt_ud(1:end-1,23)];
for i = 2:length(r_po4_txt)
    if strcmp(r_po4_txt{i,1},'') == 1
       r_po4(i) = NaN;
    elseif strcmp(r_po4_txt{i,1},'') == 0
       r_po4(i) = str2num(char(r_po4_txt{i,1}));
    end
end

r_po4_raw2=[raw2(1:end,34)];
r_po4(end+1:end+length(r_po4_raw2)) = r_po4_raw2;


for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
    r_raw_txt_c{i,1} = r_date_txt(i,:); % raw
end

k=0;
for i= length(r_date_txt)+1:length(r_date_txt)+length(r_date_txt2)
    k=k+1;
    r_date_txt_c{i,1} = [r_date_txt2(k,6:7),'.',r_date_txt2(k,9:10)]; % delete year
    r_raw_txt_c{i,1} = [r_date_txt2(k,1:4),'.',r_date_txt2(k,6:7),'.',r_date_txt2(k,9:10)]; % raw
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
        r_raw_date_c{i} = r_raw_txt_c{indx_raw{i}(1)}(6:10);  %remove year
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
regime_po4_pre = raw_po4(1:t_tick(end));
regime_ss_pre = raw_ss(1:t_tick(end));

% 2004 ~ 2019
regime_do_af_pre = raw_do(t_tick(end)+1:end);
regime_chl_af_pre = raw_chl(t_tick(end)+1:end);
regime_no3_af_pre = raw_no3(t_tick(end)+1:end);
regime_nh4_af_pre = raw_nh4(t_tick(end)+1:end);
regime_po4_af_pre = raw_po4(t_tick(end)+1:end);
regime_ss_af_pre = raw_ss(t_tick(end)+1:end);

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

clearvars idx1
idx1 = find(isnan(regime_po4_pre) == 0);
regime_po4=regime_po4_pre;
regime_po4(find(regime_po4_pre > mean(regime_po4_pre(idx1)) + sig*std(regime_po4_pre(idx1))))=NaN;
regime_po4(find(regime_po4_pre < mean(regime_po4_pre(idx1)) - sig*std(regime_po4_pre(idx1))))=NaN;
regm_po4 = nanmean(regime_po4)

clearvars idx1
idx1 = find(isnan(regime_ss_pre) == 0);
regime_ss=regime_ss_pre;
regime_ss(find(regime_ss_pre > mean(regime_ss_pre(idx1)) + sig*std(regime_ss_pre(idx1))))=NaN;
regime_ss(find(regime_ss_pre < mean(regime_ss_pre(idx1)) - sig*std(regime_ss_pre(idx1))))=NaN;
regm_ss = nanmean(regime_ss)

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

clearvars idx1
idx1 = find(isnan(regime_po4_af_pre) == 0);
regime_po4_af=regime_po4_af_pre;
regime_po4_af(find(regime_po4_af_pre > mean(regime_po4_af_pre(idx1)) + sig*std(regime_po4_af_pre(idx1))))=NaN;
regime_po4_af(find(regime_po4_af_pre < mean(regime_po4_af_pre(idx1)) - sig*std(regime_po4_af_pre(idx1))))=NaN;
regm_po4_af = nanmean(regime_po4_af)

clearvars idx1
idx1 = find(isnan(regime_ss_af_pre) == 0);
regime_ss_af=regime_ss_af_pre;
regime_ss_af(find(regime_ss_af_pre > mean(regime_ss_af_pre(idx1)) + sig*std(regime_ss_af_pre(idx1))))=NaN;
regime_ss_af(find(regime_ss_af_pre < mean(regime_ss_af_pre(idx1)) - sig*std(regime_ss_af_pre(idx1))))=NaN;
regm_ss_af = nanmean(regime_ss_af)

%% extract regime mean
% regime_do = regime_do - regm_do;
% regime_chl = regime_chl - regm_chl;
% regime_no3 = regime_no3 - regm_no3;
% regime_nh4 = regime_nh4 - regm_nh4;
% regime_po4 = regime_po4 - regm_po4;
% regime_ss = regime_ss - regm_ss;
% 
% regime_do_af = regime_do_af - (regm_do_af - regm_do) - regm_do;
% regime_chl_af = regime_chl_af - (regm_chl_af - regm_chl) - regm_chl;
% regime_no3_af = regime_no3_af - (regm_no3_af - regm_no3) - regm_no3;
% regime_nh4_af = regime_nh4_af - (regm_nh4_af - regm_nh4) - regm_nh4;
% regime_po4_af = regime_po4_af - (regm_po4_af - regm_po4) - regm_po4;
% regime_ss_af = regime_ss_af - (regm_ss_af - regm_ss) - regm_ss;

%% combine
com_do = [regime_do'; regime_do_af';] *0.7*44.661;
com_chl = [regime_chl'; regime_chl_af';];
com_no3 = [regime_no3'; regime_no3_af';] .*1000 ./14.006720;
com_nh4 = [regime_nh4'; regime_nh4_af';] .*1000 ./14.006720;
com_po4 = [regime_po4'; regime_po4_af';] .*1000 ./30.973762;
com_ss = [regime_ss'; regime_ss_af';];

return

clearvars regime_*
regime_do = com_do;
regime_chl = com_chl;
regime_no3 = com_no3;
regime_nh4 = com_nh4;
regime_po4 = com_po4;
regime_ss = com_ss;


%% monthly mean (yy.mm) form after over 3sigma value extraction.

for i=1:length(t_indx) % time for num. of months 
if i ==1
   monthly_do(i) = nanmean(regime_do(1:t_indx(i)));
   monthly_chl(i) = nanmean(regime_chl(1:t_indx(i)));
   monthly_no3(i) = nanmean(regime_no3(1:t_indx(i)));
   monthly_nh4(i) = nanmean(regime_nh4(1:t_indx(i)));
   monthly_po4(i) = nanmean(regime_po4(1:t_indx(i)));
   monthly_ss(i) = nanmean(regime_ss(1:t_indx(i)));
else
   monthly_do(i) = nanmean(regime_do(t_indx(i-1)+1:t_indx(i)));
   monthly_chl(i) = nanmean(regime_chl(t_indx(i-1)+1:t_indx(i)));
   monthly_no3(i) = nanmean(regime_no3(t_indx(i-1)+1:t_indx(i)));
   monthly_nh4(i) = nanmean(regime_nh4(t_indx(i-1)+1:t_indx(i)));
   monthly_po4(i) = nanmean(regime_po4(t_indx(i-1)+1:t_indx(i)));
   monthly_ss(i) = nanmean(regime_ss(t_indx(i-1)+1:t_indx(i)));
end
end

% fill missing value
mon_in_chl = monthly_chl;
mon_in_no3 = monthly_no3;
mon_in_nh4 = monthly_nh4;
mon_in_do = monthly_do;
mon_in_po4 = monthly_po4;
mon_in_ss = monthly_ss;

t=1:length(t_indx);
mon_in_chl(isnan(mon_in_chl))=interp1(t(~isnan(mon_in_chl)),mon_in_chl(~isnan(mon_in_chl)),t(isnan(mon_in_chl)));
mon_in_no3(isnan(mon_in_no3))=interp1(t(~isnan(mon_in_no3)),mon_in_no3(~isnan(mon_in_no3)),t(isnan(mon_in_no3)));
mon_in_nh4(isnan(mon_in_nh4))=interp1(t(~isnan(mon_in_nh4)),mon_in_nh4(~isnan(mon_in_nh4)),t(isnan(mon_in_nh4)));
mon_in_do(isnan(mon_in_do))=interp1(t(~isnan(mon_in_do)),mon_in_do(~isnan(mon_in_do)),t(isnan(mon_in_do)));
mon_in_po4(isnan(mon_in_po4))=interp1(t(~isnan(mon_in_po4)),mon_in_po4(~isnan(mon_in_po4)),t(isnan(mon_in_po4)));
mon_in_ss(isnan(mon_in_ss))=interp1(t(~isnan(mon_in_ss)),mon_in_ss(~isnan(mon_in_ss)),t(isnan(mon_in_ss)));

save('namgang_yymm_monthly_data_89to19_3sig.mat','mon*','t_indx','yymmdd_txt_c');
save('namgang_yymm_raw_data_89to19_3sig.mat','com*','t_indx','yymmdd_txt_c');