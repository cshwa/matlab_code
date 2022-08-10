close all; clear; clc;   % -v3
cd C:\Users\user\Desktop\장기생태_하수_여수월내
% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port

name_tag = {'광양산업단지폐수종말처리장', '광양중앙하수종말처리장','광양하수종말처리장',...
    '광영하수종말처리장','여수월내산업단지폐수종말처리장','여수율촌산업단지폐수종말처리장',...
    '여수중흥산업단지폐수종말처리장','여수하수종말처리장','진월하수종말처리장' } 

% combining the tag and outter point excluding
size_tag = length(name_tag);

%% pick the row on the excel which has same name with tag
% col 1: {'처리시설명'}
% col 3: {'운영일자'}
% col 5: {'방류(연계)유량 및 방류(연계)농도_유량(㎥/일)_물리적'}
% col 10: {'방류(연계)유량 및 방류(연계)농도_SS(mg/ℓ)'}
% col 11: {'방류(연계)유량 및 방류(연계)농도_TN(mg/ℓ)'}
% col 12: {'방류(연계)유량 및 방류(연계)농도_TP(mg/ℓ)'}

%광양산업단지폐수종말처리장
[raw_g_12 txt_g_12]=xlsread('2012_전국오염원조사_일_9.xlsx',name_tag{1},'');
[raw_g_13 txt_g_13]=xlsread('2013_전국오염원조사_일_9.xlsx',name_tag{1},'');
[raw_g_14 txt_g_14]=xlsread('2014_전국오염원조사_일_9.xlsx',name_tag{1},'');
[raw_g_15 txt_g_15]=xlsread('2015_전국오염원조사_일_9.xlsx',name_tag{1},'');
[raw_g_16 txt_g_16]=xlsread('2016_전국오염원조사_일_9.xlsx',name_tag{1},'');
[raw_g_17 txt_g_17]=xlsread('2017_전국오염원조사_일_9.xlsx',name_tag{1},'');
[raw_g_18 txt_g_18]=xlsread('2018_전국오염원조사_일_9.xlsx',name_tag{1},'');
[raw_g_19 txt_g_19]=xlsread('2019_전국오염원조사_일_9.xlsx',name_tag{1},'');

%광양중앙하수종말처리장
[raw_gc_12 txt_gc_12]=xlsread('2012_전국오염원조사_일_9.xlsx',name_tag{2},'');
[raw_gc_13 txt_gc_13]=xlsread('2013_전국오염원조사_일_9.xlsx',name_tag{2},'');
[raw_gc_14 txt_gc_14]=xlsread('2014_전국오염원조사_일_9.xlsx',name_tag{2},'');
[raw_gc_15 txt_gc_15]=xlsread('2015_전국오염원조사_일_9.xlsx',name_tag{2},'');
[raw_gc_16 txt_gc_16]=xlsread('2016_전국오염원조사_일_9.xlsx',name_tag{2},'');
[raw_gc_17 txt_gc_17]=xlsread('2017_전국오염원조사_일_9.xlsx',name_tag{2},'');
[raw_gc_18 txt_gc_18]=xlsread('2018_전국오염원조사_일_9.xlsx',name_tag{2},'');
[raw_gc_19 txt_gc_19]=xlsread('2019_전국오염원조사_일_9.xlsx',name_tag{2},'');

% 광양하수종말처리장
[raw_gh_12 txt_gh_12]=xlsread('2012_전국오염원조사_일_9.xlsx',name_tag{3},'');
[raw_gh_13 txt_gh_13]=xlsread('2013_전국오염원조사_일_9.xlsx',name_tag{3},'');
[raw_gh_14 txt_gh_14]=xlsread('2014_전국오염원조사_일_9.xlsx',name_tag{3},'');
[raw_gh_15 txt_gh_15]=xlsread('2015_전국오염원조사_일_9.xlsx',name_tag{3},'');
[raw_gh_16 txt_gh_16]=xlsread('2016_전국오염원조사_일_9.xlsx',name_tag{3},'');
[raw_gh_17 txt_gh_17]=xlsread('2017_전국오염원조사_일_9.xlsx',name_tag{3},'');
[raw_gh_18 txt_gh_18]=xlsread('2018_전국오염원조사_일_9.xlsx',name_tag{3},'');
[raw_gh_19 txt_gh_19]=xlsread('2019_전국오염원조사_일_9.xlsx',name_tag{3},'');

% 광영하수종말처리장
[raw_gh2_12 txt_gh2_12]=xlsread('2012_전국오염원조사_일_9.xlsx',name_tag{4},'');
[raw_gh2_13 txt_gh2_13]=xlsread('2013_전국오염원조사_일_9.xlsx',name_tag{4},'');
[raw_gh2_14 txt_gh2_14]=xlsread('2014_전국오염원조사_일_9.xlsx',name_tag{4},'');
[raw_gh2_15 txt_gh2_15]=xlsread('2015_전국오염원조사_일_9.xlsx',name_tag{4},'');
[raw_gh2_16 txt_gh2_16]=xlsread('2016_전국오염원조사_일_9.xlsx',name_tag{4},'');
[raw_gh2_17 txt_gh2_17]=xlsread('2017_전국오염원조사_일_9.xlsx',name_tag{4},'');
[raw_gh2_18 txt_gh2_18]=xlsread('2018_전국오염원조사_일_9.xlsx',name_tag{4},'');
[raw_gh2_19 txt_gh2_19]=xlsread('2019_전국오염원조사_일_9.xlsx',name_tag{4},'');

%여수월내산업단지폐수종말처리장
[raw_w_12 txt_w_12]=xlsread('2012_전국오염원조사_일_9.xlsx',name_tag{5},'');
[raw_w_13 txt_w_13]=xlsread('2013_전국오염원조사_일_9.xlsx',name_tag{5},'');
[raw_w_14 txt_w_14]=xlsread('2014_전국오염원조사_일_9.xlsx',name_tag{5},'');
[raw_w_15 txt_w_15]=xlsread('2015_전국오염원조사_일_9.xlsx',name_tag{5},'');
[raw_w_16 txt_w_16]=xlsread('2016_전국오염원조사_일_9.xlsx',name_tag{5},'');
[raw_w_17 txt_w_17]=xlsread('2017_전국오염원조사_일_9.xlsx',name_tag{5},'');
[raw_w_18 txt_w_18]=xlsread('2018_전국오염원조사_일_9.xlsx',name_tag{5},'');
[raw_w_19 txt_w_19]=xlsread('2019_전국오염원조사_일_9.xlsx',name_tag{5},'');

%여수율촌산업단지폐수종말처리장
[raw_y_12 txt_y_12]=xlsread('2012_전국오염원조사_일_9.xlsx',name_tag{6},'');
[raw_y_13 txt_y_13]=xlsread('2013_전국오염원조사_일_9.xlsx',name_tag{6},'');
[raw_y_14 txt_y_14]=xlsread('2014_전국오염원조사_일_9.xlsx',name_tag{6},'');
[raw_y_15 txt_y_15]=xlsread('2015_전국오염원조사_일_9.xlsx',name_tag{6},'');
[raw_y_16 txt_y_16]=xlsread('2016_전국오염원조사_일_9.xlsx',name_tag{6},'');
[raw_y_17 txt_y_17]=xlsread('2017_전국오염원조사_일_9.xlsx',name_tag{6},'');
[raw_y_18 txt_y_18]=xlsread('2018_전국오염원조사_일_9.xlsx',name_tag{6},'');
[raw_y_19 txt_y_19]=xlsread('2019_전국오염원조사_일_9.xlsx',name_tag{6},'');

%여수중흥산업단지폐수종말처리장
[raw_j_12 txt_j_12]=xlsread('2012_전국오염원조사_일_9.xlsx',name_tag{7},'');
[raw_j_13 txt_j_13]=xlsread('2013_전국오염원조사_일_9.xlsx',name_tag{7},'');
[raw_j_14 txt_j_14]=xlsread('2014_전국오염원조사_일_9.xlsx',name_tag{7},'');
[raw_j_15 txt_j_15]=xlsread('2015_전국오염원조사_일_9.xlsx',name_tag{7},'');
[raw_j_16 txt_j_16]=xlsread('2016_전국오염원조사_일_9.xlsx',name_tag{7},'');
[raw_j_17 txt_j_17]=xlsread('2017_전국오염원조사_일_9.xlsx',name_tag{7},'');
[raw_j_18 txt_j_18]=xlsread('2018_전국오염원조사_일_9.xlsx',name_tag{7},'');
[raw_j_19 txt_j_19]=xlsread('2019_전국오염원조사_일_9.xlsx',name_tag{7},'');

%여수하수종말처리장
[raw_ye_12 txt_ye_12]=xlsread('2012_전국오염원조사_일_9.xlsx',name_tag{8},'');
[raw_ye_13 txt_ye_13]=xlsread('2013_전국오염원조사_일_9.xlsx',name_tag{8},'');
[raw_ye_14 txt_ye_14]=xlsread('2014_전국오염원조사_일_9.xlsx',name_tag{8},'');
[raw_ye_15 txt_ye_15]=xlsread('2015_전국오염원조사_일_9.xlsx',name_tag{8},'');
[raw_ye_16 txt_ye_16]=xlsread('2016_전국오염원조사_일_9.xlsx',name_tag{8},'');
[raw_ye_17 txt_ye_17]=xlsread('2017_전국오염원조사_일_9.xlsx',name_tag{8},'');
[raw_ye_18 txt_ye_18]=xlsread('2018_전국오염원조사_일_9.xlsx',name_tag{8},'');
[raw_ye_19 txt_ye_19]=xlsread('2019_전국오염원조사_일_9.xlsx',name_tag{8},'');

%진월하수종말처리장
[raw_jw_12 txt_jw_12]=xlsread('2012_전국오염원조사_일_9.xlsx',name_tag{9},'');
[raw_jw_13 txt_jw_13]=xlsread('2013_전국오염원조사_일_9.xlsx',name_tag{9},'');
[raw_jw_14 txt_jw_14]=xlsread('2014_전국오염원조사_일_9.xlsx',name_tag{9},'');
[raw_jw_15 txt_jw_15]=xlsread('2015_전국오염원조사_일_9.xlsx',name_tag{9},'');
[raw_jw_16 txt_jw_16]=xlsread('2016_전국오염원조사_일_9.xlsx',name_tag{9},'');
[raw_jw_17 txt_jw_17]=xlsread('2017_전국오염원조사_일_9.xlsx',name_tag{9},'');
[raw_jw_18 txt_jw_18]=xlsread('2018_전국오염원조사_일_9.xlsx',name_tag{9},'');
[raw_jw_19 txt_jw_19]=xlsread('2019_전국오염원조사_일_9.xlsx',name_tag{9},'');

% [raw_12 txt_12]=xlsread('2012_전국오염원조사_일_9.xlsx',name_tag{1},'');

%광양산업단지폐수종말처리장
merge_g_discharge = [raw_g_12(:,4); raw_g_13(:,4); raw_g_14(:,4); raw_g_15(:,4); raw_g_16(:,4); raw_g_17(:,4); raw_g_18(:,4); raw_g_19(:,6);];
merge_g_tn = [raw_g_12(:,8); raw_g_13(:,8); raw_g_14(:,8); raw_g_15(:,8); raw_g_16(:,8); raw_g_17(:,8); raw_g_18(:,8); raw_g_19(:,10);];
merge_g_tp = [raw_g_12(:,9); raw_g_13(:,9); raw_g_14(:,9); raw_g_15(:,9); raw_g_16(:,9); raw_g_17(:,9); raw_g_18(:,9); raw_g_19(:,11);];
merge_g_ss = [raw_g_12(:,7); raw_g_13(:,7); raw_g_14(:,7); raw_g_15(:,7); raw_g_16(:,7); raw_g_17(:,7); raw_g_18(:,7); raw_g_19(:,9);];
raw_g_19_pre=num2str(raw_g_19(:,2));
raw_g_19_re_yy=num2str(raw_g_19_pre(:,1:4));
raw_g_19_re_mm=num2str(raw_g_19_pre(:,5:6));
raw_g_19_re_dd=num2str(raw_g_19_pre(:,7:8));

%광양중앙하수종말처리장
merge_gc_discharge = [raw_gc_12(:,4); raw_gc_13(:,4); raw_gc_14(:,4); raw_gc_15(:,4); raw_gc_16(:,4); raw_gc_17(:,4); raw_gc_18(:,4); raw_gc_19(:,6);];
merge_gc_tn = [raw_gc_12(:,8); raw_gc_13(:,8); raw_gc_14(:,8); raw_gc_15(:,8); raw_gc_16(:,8); raw_gc_17(:,8); raw_gc_18(:,8); raw_gc_19(:,10);];
merge_gc_tp = [raw_gc_12(:,9); raw_gc_13(:,9); raw_gc_14(:,9); raw_gc_15(:,9); raw_gc_16(:,9); raw_gc_17(:,9); raw_gc_18(:,9); raw_gc_19(:,11);];
merge_gc_ss = [raw_gc_12(:,7); raw_gc_13(:,7); raw_gc_14(:,7); raw_gc_15(:,7); raw_gc_16(:,7); raw_gc_17(:,7); raw_gc_18(:,7); raw_gc_19(:,9);];
raw_gc_19_pre=num2str(raw_gc_19(:,2));
raw_gc_19_re_yy=num2str(raw_gc_19_pre(:,1:4));
raw_gc_19_re_mm=num2str(raw_gc_19_pre(:,5:6));
raw_gc_19_re_dd=num2str(raw_gc_19_pre(:,7:8));

%광양하수종말처리장
merge_gh_discharge = [raw_gh_12(:,4); raw_gh_13(:,4); raw_gh_14(:,4); raw_gh_15(:,4); raw_gh_16(:,4); raw_gh_17(:,4); raw_gh_18(:,4); raw_gh_19(:,6-1);];
merge_gh_tn = [raw_gh_12(:,8); raw_gh_13(:,8); raw_gh_14(:,8); raw_gh_15(:,8); raw_gh_16(:,8); raw_gh_17(:,8); raw_gh_18(:,8); raw_gh_19(:,10-1);];
merge_gh_tp = [raw_gh_12(:,9); raw_gh_13(:,9); raw_gh_14(:,9); raw_gh_15(:,9); raw_gh_16(:,9); raw_gh_17(:,9); raw_gh_18(:,9); raw_gh_19(:,11-1);];
merge_gh_ss = [raw_gh_12(:,7); raw_gh_13(:,7); raw_gh_14(:,7); raw_gh_15(:,7); raw_gh_16(:,7); raw_gh_17(:,7); raw_gh_18(:,7); raw_gh_19(:,9-1);];
raw_gh_19_pre=num2str(raw_gh_19(:,2-1));
raw_gh_19_re_yy=num2str(raw_gh_19_pre(:,1:4));
raw_gh_19_re_mm=num2str(raw_gh_19_pre(:,5:6));
raw_gh_19_re_dd=num2str(raw_gh_19_pre(:,7:8));

%광영하수종말처리장
merge_gh2_discharge = [raw_gh2_12(:,4); raw_gh2_13(:,4); raw_gh2_14(:,4); raw_gh2_15(:,4); raw_gh2_16(:,4); raw_gh2_17(:,4); raw_gh2_18(:,4); raw_gh2_19(:,6);];
merge_gh2_tn = [raw_gh2_12(:,8); raw_gh2_13(:,8); raw_gh2_14(:,8); raw_gh2_15(:,8); raw_gh2_16(:,8); raw_gh2_17(:,8); raw_gh2_18(:,8); raw_gh2_19(:,10);];
merge_gh2_tp = [raw_gh2_12(:,9); raw_gh2_13(:,9); raw_gh2_14(:,9); raw_gh2_15(:,9); raw_gh2_16(:,9); raw_gh2_17(:,9); raw_gh2_18(:,9); raw_gh2_19(:,11);];
merge_gh2_ss = [raw_gh2_12(:,7); raw_gh2_13(:,7); raw_gh2_14(:,7); raw_gh2_15(:,7); raw_gh2_16(:,7); raw_gh2_17(:,7); raw_gh2_18(:,7); raw_gh2_19(:,9);];
raw_gh2_19_pre=num2str(raw_gh2_19(:,2));
raw_gh2_19_re_yy=num2str(raw_gh2_19_pre(:,1:4));
raw_gh2_19_re_mm=num2str(raw_gh2_19_pre(:,5:6));
raw_gh2_19_re_dd=num2str(raw_gh2_19_pre(:,7:8));

%여수월내산업단지폐수종말처리장
merge_w_discharge = [raw_w_12(:,4); raw_w_13(:,4); raw_w_14(:,4); raw_w_15(:,4); raw_w_16(:,4); raw_w_17(:,4); raw_w_18(:,4); raw_w_19(:,6);];
merge_w_tn = [raw_w_12(:,8); raw_w_13(:,8); raw_w_14(:,8); raw_w_15(:,8); raw_w_16(:,8); raw_w_17(:,8); raw_w_18(:,8); raw_w_19(:,10);];
merge_w_tp = [raw_w_12(:,9); raw_w_13(:,9); raw_w_14(:,9); raw_w_15(:,9); raw_w_16(:,9); raw_w_17(:,9); raw_w_18(:,9); raw_w_19(:,11);];
merge_w_ss = [raw_w_12(:,7); raw_w_13(:,7); raw_w_14(:,7); raw_w_15(:,7); raw_w_16(:,7); raw_w_17(:,7); raw_w_18(:,7); raw_w_19(:,9);];
raw_w_19_pre=num2str(raw_w_19(:,2));
raw_w_19_re_yy=num2str(raw_w_19_pre(:,1:4));
raw_w_19_re_mm=num2str(raw_w_19_pre(:,5:6));
raw_w_19_re_dd=num2str(raw_w_19_pre(:,7:8));

% 여수율촌산업단지폐수종말처리장
merge_y_discharge = [raw_y_12(:,4); raw_y_13(:,4); raw_y_14(:,4); raw_y_15(:,4); raw_y_16(:,4); raw_y_17(:,4); raw_y_18(:,4); raw_y_19(:,6-1);];
merge_y_tn = [raw_y_12(:,8); raw_y_13(:,8); raw_y_14(:,8); raw_y_15(:,8); raw_y_16(:,8); raw_y_17(:,8); raw_y_18(:,8); raw_y_19(:,10-1);];
merge_y_tp = [raw_y_12(:,9); raw_y_13(:,9); raw_y_14(:,9); raw_y_15(:,9); raw_y_16(:,9); raw_y_17(:,9); raw_y_18(:,9); raw_y_19(:,11-1);];
merge_y_ss = [raw_y_12(:,7); raw_y_13(:,7); raw_y_14(:,7); raw_y_15(:,7); raw_y_16(:,7); raw_y_17(:,7); raw_y_18(:,7); raw_y_19(:,9-1);];
raw_y_19_pre=num2str(raw_y_19(:,1));
raw_y_19_re_yy=num2str(raw_y_19_pre(:,1:4));
raw_y_19_re_mm=num2str(raw_y_19_pre(:,5:6));
raw_y_19_re_dd=num2str(raw_y_19_pre(:,7:8));

%여수중흥산업단지폐수종말처리장
merge_j_discharge = [raw_j_12(:,4); raw_j_13(:,4); raw_j_14(:,4); raw_j_15(:,4); raw_j_16(:,4); raw_j_17(:,4); raw_j_18(:,4); raw_j_19(:,6);];
merge_j_tn = [raw_j_12(:,8); raw_j_13(:,8); raw_j_14(:,8); raw_j_15(:,8); raw_j_16(:,8); raw_j_17(:,8); raw_j_18(:,8); raw_j_19(:,10);];
merge_j_tp = [raw_j_12(:,9); raw_j_13(:,9); raw_j_14(:,9); raw_j_15(:,9); raw_j_16(:,9); raw_j_17(:,9); raw_j_18(:,9); raw_j_19(:,11);];
merge_j_ss = [raw_j_12(:,7); raw_j_13(:,7); raw_j_14(:,7); raw_j_15(:,7); raw_j_16(:,7); raw_j_17(:,7); raw_j_18(:,7); raw_j_19(:,9);];
raw_j_19_pre=num2str(raw_j_19(:,2));
raw_j_19_re_yy=num2str(raw_j_19_pre(:,1:4));
raw_j_19_re_mm=num2str(raw_j_19_pre(:,5:6));
raw_j_19_re_dd=num2str(raw_j_19_pre(:,7:8));

%여수하수종말처리장
merge_ye_discharge = [raw_ye_12(:,4); raw_ye_13(:,4); raw_ye_14(:,4); raw_ye_15(:,4); raw_ye_16(:,4); raw_ye_17(:,4); raw_ye_18(:,4); raw_ye_19(:,6);];
merge_ye_tn = [raw_ye_12(:,8); raw_ye_13(:,8); raw_ye_14(:,8); raw_ye_15(:,8); raw_ye_16(:,8); raw_ye_17(:,8); raw_ye_18(:,8); raw_ye_19(:,10);];
merge_ye_tp = [raw_ye_12(:,9); raw_ye_13(:,9); raw_ye_14(:,9); raw_ye_15(:,9); raw_ye_16(:,9); raw_ye_17(:,9); raw_ye_18(:,9); raw_ye_19(:,11);];
merge_ye_ss = [raw_ye_12(:,7); raw_ye_13(:,7); raw_ye_14(:,7); raw_ye_15(:,7); raw_ye_16(:,7); raw_ye_17(:,7); raw_ye_18(:,7); raw_ye_19(:,9);];
raw_ye_19_pre=num2str(raw_ye_19(:,2));
raw_ye_19_re_yy=num2str(raw_ye_19_pre(:,1:4));
raw_ye_19_re_mm=num2str(raw_ye_19_pre(:,5:6));
raw_ye_19_re_dd=num2str(raw_ye_19_pre(:,7:8));

%진월하수종말처리장
merge_jw_discharge = [raw_jw_12(:,4); raw_jw_13(:,4); raw_jw_14(:,4); raw_jw_15(:,4); raw_jw_16(:,4); raw_jw_17(:,4); raw_jw_18(:,4); raw_jw_19(:,6-1);];
merge_jw_tn = [raw_jw_12(:,8); raw_jw_13(:,8); raw_jw_14(:,8); raw_jw_15(:,8); raw_jw_16(:,8); raw_jw_17(:,8); raw_jw_18(:,8); raw_jw_19(:,10-1);];
merge_jw_tp = [raw_jw_12(:,9); raw_jw_13(:,9); raw_jw_14(:,9); raw_jw_15(:,9); raw_jw_16(:,9); raw_jw_17(:,9); raw_jw_18(:,9); raw_jw_19(:,11-1);];
merge_jw_ss = [raw_jw_12(:,7); raw_jw_13(:,7); raw_jw_14(:,7); raw_jw_15(:,7); raw_jw_16(:,7); raw_jw_17(:,7); raw_jw_18(:,7); raw_jw_19(:,9-1);];
raw_jw_19_pre=num2str(raw_jw_19(:,2-1));
raw_jw_19_re_yy=num2str(raw_jw_19_pre(:,1:4));
raw_jw_19_re_mm=num2str(raw_jw_19_pre(:,5:6));
raw_jw_19_re_dd=num2str(raw_jw_19_pre(:,7:8));

%광양산업단지폐수종말처리장
clearvars raw_g_19_re
for i = 1:size(raw_g_19_pre,1)
raw_g_19_re{i,1} = {[raw_g_19_re_yy(i,:),'-',raw_g_19_re_mm(i,:),'-' ,raw_g_19_re_dd(i,:)]};
end

%광양중앙하수종말처리장
clearvars raw_gc_19_re
for i = 1:size(raw_gc_19_pre,1)
raw_gc_19_re{i,1} = {[raw_gc_19_re_yy(i,:),'-',raw_gc_19_re_mm(i,:),'-' ,raw_gc_19_re_dd(i,:)]};
end

%광양하수종말처리장
clearvars raw_gh_19_re
for i = 1:size(raw_gh_19_pre,1)
raw_gh_19_re{i,1} = {[raw_gh_19_re_yy(i,:),'-',raw_gh_19_re_mm(i,:),'-' ,raw_gh_19_re_dd(i,:)]};
end

%광영하수종말처리장
clearvars raw_gh2_19_re
for i = 1:size(raw_gh2_19_pre,1)
raw_gh2_19_re{i,1} = {[raw_gh2_19_re_yy(i,:),'-',raw_gh2_19_re_mm(i,:),'-' ,raw_gh2_19_re_dd(i,:)]};
end

%여수월내산업단지폐수종말처리장
clearvars raw_w_19_re
for i = 1:size(raw_w_19_pre,1)
raw_w_19_re{i,1} = {[raw_w_19_re_yy(i,:),'-',raw_w_19_re_mm(i,:),'-' ,raw_w_19_re_dd(i,:)]};
end

% 여수율촌산업단지폐수종말처리장
clearvars raw_y_19_re
for i = 1:size(raw_y_19_pre,1)
raw_y_19_re{i,1} = {[raw_y_19_re_yy(i,:),'-',raw_y_19_re_mm(i,:),'-' ,raw_y_19_re_dd(i,:)]};
end

%여수중흥산업단지폐수종말처리장
clearvars raw_j_19_re
for i = 1:size(raw_j_19_pre,1)
raw_j_19_re{i,1} = {[raw_j_19_re_yy(i,:),'-',raw_j_19_re_mm(i,:),'-' ,raw_j_19_re_dd(i,:)]};
end

%여수하수종말처리장
clearvars raw_ye_19_re
for i = 1:size(raw_ye_19_pre,1)
raw_ye_19_re{i,1} = {[raw_ye_19_re_yy(i,:),'-',raw_ye_19_re_mm(i,:),'-' ,raw_ye_19_re_dd(i,:)]};
end

%진월하수종말처리장
clearvars raw_jw_19_re
for i = 1:size(raw_jw_19_pre,1)
raw_jw_19_re{i,1} = {[raw_jw_19_re_yy(i,:),'-',raw_jw_19_re_mm(i,:),'-' ,raw_jw_19_re_dd(i,:)]};
end

%% merge date

%광양산업단지폐수종말처리장
clearvars merge_g_date
merge_g_date = [txt_g_12(2:end,3); txt_g_13(2:end,7); txt_g_14(2:end,7); txt_g_15(2:end,7); txt_g_16(2:end,7); txt_g_17(2:end,7); txt_g_18(2:end,7); raw_g_19_re(:,1);];

%광양중앙하수종말처리장
clearvars merge_gc_date
merge_gc_date = [txt_gc_12(2:end,3); txt_gc_13(2:end,7); txt_gc_14(2:end,7); txt_gc_15(2:end,7); txt_gc_16(2:end,7); txt_gc_17(2:end,7); txt_gc_18(2:end,7); raw_gc_19_re(:,1);];

% 광양하수종말처리장
clearvars merge_gh_date
merge_gh_date = [txt_gh_12(2:end,3); txt_gh_13(2:end,7); txt_gh_14(2:end,7); txt_gh_15(2:end,7); txt_gh_16(2:end,7); txt_gh_17(2:end,7); txt_gh_18(2:end,7); raw_gh_19_re(:,1);];

% 광영하수종말처리장
clearvars merge_gh2_date
merge_gh2_date = [txt_gh2_12(2:end,3); txt_gh2_13(2:end,7); txt_gh2_14(2:end,7); txt_gh2_15(2:end,7); txt_gh2_16(2:end,7); txt_gh2_17(2:end,7); txt_gh2_18(2:end,7); raw_gh2_19_re(:,1);];

%여수월내산업단지폐수종말처리장
clearvars merge_w_date
merge_w_date = [txt_w_12(2:end,3); txt_w_13(2:end,7); txt_w_14(2:end,7); txt_w_15(2:end,7); txt_w_16(2:end,7); txt_w_17(2:end,7); txt_w_18(2:end,7); raw_w_19_re(:,1);];
% 2922 days => 2012:2019

% 여수율촌산업단지폐수종말처리장
merge_y_date = [txt_y_12(2:end,3); txt_y_13(2:end,7); txt_y_14(2:end,7); txt_y_15(2:end,7); txt_y_16(2:end,7); txt_y_17(2:end,7); txt_y_18(2:end,7); raw_y_19_re(:,1);];

%여수중흥산업단지폐수종말처리장
clearvars merge_j_date
merge_j_date = [txt_j_12(2:end,3); txt_j_13(2:end,7); txt_j_14(2:end,7); txt_j_15(2:end,7); txt_j_16(2:end,7); txt_j_17(2:end,7); txt_j_18(2:end,7); raw_j_19_re(:,1);];

%여수하수종말처리장
clearvars merge_ye_date
merge_ye_date = [txt_ye_12(2:end,3); txt_ye_13(2:end,7); txt_ye_14(2:end,7); txt_ye_15(2:end,7); txt_ye_16(2:end,7); txt_ye_17(2:end,7); txt_ye_18(2:end,7); raw_ye_19_re(:,1);];

%진월하수종말처리장
clearvars merge_jw_date
merge_jw_date = [txt_jw_12(2:end,3); txt_jw_13(2:end,7); txt_jw_14(2:end,7); txt_jw_15(2:end,7); txt_jw_16(2:end,7); txt_jw_17(2:end,7); txt_jw_18(2:end,7); raw_jw_19_re(:,1);];


%% pick matched name with tag
% %port
% for i = 1:length(name_tag)
%    if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
%        indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)     
%    end
% end

%% make date to be 'yymm' form
% 광양산업단지폐수종말처리장
clearvars temp
for i = 1:length(merge_g_date)
temp = char(merge_g_date{i});
if size(temp) ~= 7
    temp=temp(1,1:7);
end
merge_g_yymm{i,1} = temp;
end

% 광양중앙하수종말처리장
clearvars temp
for i = 1:length(merge_gc_date)
temp = char(merge_gc_date{i});
if size(temp) ~= 7
    temp=temp(1,1:7);
end
merge_gc_yymm{i,1} = temp;
end

% 광양하수종말처리장
clearvars temp
for i = 1:length(merge_gh_date)
temp = char(merge_gh_date{i});
if size(temp) ~= 7
    temp=temp(1,1:7);
end
merge_gh_yymm{i,1} = temp;
end

% 광영하수종말처리장
for i = 1:length(merge_gh2_date)
temp = char(merge_gh2_date{i});
if size(temp) ~= 7
    temp=temp(1,1:7);
end
merge_gh2_yymm{i,1} = temp;
end

% 여수월내산업단지폐수종말처리장
clearvars temp
for i = 1:length(merge_w_date)
temp = char(merge_w_date{i});
if size(temp) ~= 7
    temp=temp(1,1:7);
end
merge_yymm{i,1} = temp;
end

% 여수율촌산업단지폐수종말처리장
clearvars temp_y
for i = 1:length(merge_y_date)
temp_y = char(merge_y_date{i});
if size(temp_y) ~= 7
    temp_y=temp_y(1,1:7);
end
merge_y_yymm{i,1} = temp_y;
end

% 여수중흥산업단지폐수종말처리장
clearvars temp
for i = 1:length(merge_j_date)
temp = char(merge_j_date{i});
if size(temp) ~= 7
    temp=temp(1,1:7);
end
merge_j_yymm{i,1} = temp;
end

% 여수하수종말처리장
clearvars temp
for i = 1:length(merge_ye_date)
temp = char(merge_ye_date{i});
if size(temp) ~= 7
    temp=temp(1,1:7);
end
merge_ye_yymm{i,1} = temp;
end

% 진월하수종말처리장
clearvars temp
for i = 1:length(merge_jw_date)
temp = char(merge_jw_date{i});
if size(temp) ~= 7
    temp=temp(1,1:7);
end
merge_jw_yymm{i,1} = temp;
end


% %% make date to be 'mm' form
% for i = 1:length(merge_date)
% temp = char(merge_date{i});
% temp = temp(1,6:7);
% merge_mm{i,1} = temp;
% end

%% make 1997 to 2018 'yymm' form
k=0
for i = 1997:2020
    for j = 1:12
        k=k+1;
        ref_date{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end

% %% make 1997 to 2018 'yymm' form
% for j = 1:12
%  ref_date_mm{j,1} = [num2str(j,'%02d')];
% end

%% matched date 'yymm' form

%광양산업단지폐수종말처리장
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_g_yymm)) ~= 0
           indx_date_g{i} = find([strcmp(ref_date{i}, merge_g_yymm)] == 1);     
       end
    end


%광양중앙하수종말처리장
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_gc_yymm)) ~= 0
           indx_date_gc{i} = find([strcmp(ref_date{i}, merge_gc_yymm)] == 1);     
       end
    end


%광양하수종말처리장
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_gh_yymm)) ~= 0
           indx_date_gh{i} = find([strcmp(ref_date{i}, merge_gh_yymm)] == 1);     
       end
    end
    
% 광영하수종말처리장
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_gh2_yymm)) ~= 0
           indx_date_gh2{i} = find([strcmp(ref_date{i}, merge_gh2_yymm)] == 1);     
       end
    end
    

%여수월내산업단지폐수종말처리장
% for j = 1:length(indx) % st. axis
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_yymm)) ~= 0
           indx_date{i} = find([strcmp(ref_date{i}, merge_yymm)] == 1);     
       end
    end
% end

%여수율촌산업단지폐수종말처리장
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_y_yymm)) ~= 0
           indx_date_y{i} = find([strcmp(ref_date{i}, merge_y_yymm)] == 1);     
       end
    end

%여수중흥산업단지폐수종말처리장
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_j_yymm)) ~= 0
           indx_date_j{i} = find([strcmp(ref_date{i}, merge_j_yymm)] == 1);     
       end
    end

%여수하수종말처리장
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_ye_yymm)) ~= 0
           indx_date_ye{i} = find([strcmp(ref_date{i}, merge_ye_yymm)] == 1);     
       end
    end

%진월하수종말처리장
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_jw_yymm)) ~= 0
           indx_date_jw{i} = find([strcmp(ref_date{i}, merge_jw_yymm)] == 1);     
       end
    end


% % matched date 'mm' form
% for j = 1:length(indx) % st. axis
%     for i = 1:length(ref_date_mm) % date axis
%        if  sum(strcmp(ref_date_mm{i}, merge_mm(indx{j}))) ~= 0
%            indx_date_mm{j,i} = find([strcmp(ref_date_mm{i}, merge_mm(indx{j}))] == 1);     
%        end
%     end
% end

% merge_w_discharge
% merge_w_tn
% merge_w_tp
% merge_w_ss 

%광양산업단지폐수종말처리장
clearvars discharge_g_yymm
for i = 1:length(indx_date_g)
        discharge_g_yymm(i) = nanmean(merge_g_discharge(indx_date_g{i}));
end

 clearvars tn_g_yymm
for i = 1:length(indx_date_g)
        tn_g_yymm(i) = nanmean(merge_g_tn(indx_date_g{i}));
end
 
clearvars tp_g_yymm
for i = 1:length(indx_date_g)
        tp_g_yymm(i) = nanmean(merge_g_tp(indx_date_g{i}));
end
 
clearvars ss_g_yymm
for i = 1:length(indx_date_g)
        ss_g_yymm(i) = nanmean(merge_g_ss(indx_date_g{i}));
end

%광양중앙하수종말처리장
clearvars discharge_gc_yymm
for i = 1:length(indx_date_gc)
        discharge_gc_yymm(i) = nanmean(merge_gc_discharge(indx_date_gc{i}));
end

 clearvars tn_gc_yymm
for i = 1:length(indx_date_gc)
        tn_gc_yymm(i) = nanmean(merge_gc_tn(indx_date_gc{i}));
end
 
clearvars tp_gc_yymm
for i = 1:length(indx_date_gc)
        tp_gc_yymm(i) = nanmean(merge_gc_tp(indx_date_gc{i}));
end
 
clearvars ss_gc_yymm
for i = 1:length(indx_date_gc)
        ss_gc_yymm(i) = nanmean(merge_gc_ss(indx_date_gc{i}));
end


%광양하수종말처리장
clearvars discharge_gh_yymm
for i = 1:length(indx_date_gh)
        discharge_gh_yymm(i) = nanmean(merge_gh_discharge(indx_date_gh{i}));
end

 clearvars tn_gh_yymm
for i = 1:length(indx_date_gh)
        tn_gh_yymm(i) = nanmean(merge_gh_tn(indx_date_gh{i}));
end
 
clearvars tp_gh_yymm
for i = 1:length(indx_date_gh)
        tp_gh_yymm(i) = nanmean(merge_gh_tp(indx_date_gh{i}));
end
 
clearvars ss_gh_yymm
for i = 1:length(indx_date_gh)
        ss_gh_yymm(i) = nanmean(merge_gh_ss(indx_date_gh{i}));
end

%광영하수종말처리장
clearvars discharge_gh2_yymm
for i = 1:length(indx_date_gh2)
        discharge_gh2_yymm(i) = nanmean(merge_gh2_discharge(indx_date_gh2{i}));
end

 clearvars tn_gh2_yymm
for i = 1:length(indx_date_gh2)
        tn_gh2_yymm(i) = nanmean(merge_gh2_tn(indx_date_gh2{i}));
end
 
clearvars tp_gh2_yymm
for i = 1:length(indx_date_gh2)
        tp_gh2_yymm(i) = nanmean(merge_gh2_tp(indx_date_gh2{i}));
end
 
clearvars ss_gh2_yymm
for i = 1:length(indx_date_gh2)
        ss_gh2_yymm(i) = nanmean(merge_gh2_ss(indx_date_gh2{i}));
end

%여수월내산업단지폐수종말처리장
%discharge_yymm
clearvars discharge_yymm
for i = 1:length(indx_date)
        discharge_yymm(i) = nanmean(merge_w_discharge(indx_date{i}));
end

 clearvars tn_yymm
for i = 1:length(indx_date)
        tn_yymm(i) = nanmean(merge_w_tn(indx_date{i}));
end
 
clearvars tp_yymm
for i = 1:length(indx_date)
        tp_yymm(i) = nanmean(merge_w_tp(indx_date{i}));
end
 
clearvars ss_yymm
for i = 1:length(indx_date)
        ss_yymm(i) = nanmean(merge_w_ss(indx_date{i}));
end

% 여수율촌산업단지폐수종말처리장
clearvars discharge_y_yymm
for i = 1:length(indx_date_y)
        discharge_y_yymm(i) = nanmean(merge_y_discharge(indx_date_y{i}));
end

 clearvars tn_y_yymm
for i = 1:length(indx_date_y)
        tn_y_yymm(i) = nanmean(merge_y_tn(indx_date_y{i}));
end
 
clearvars tp_y_yymm
for i = 1:length(indx_date_y)
        tp_y_yymm(i) = nanmean(merge_y_tp(indx_date_y{i}));
end
 
clearvars ss_y_yymm
for i = 1:length(indx_date_y)
        ss_y_yymm(i) = nanmean(merge_y_ss(indx_date_y{i}));
end

%여수중흥산업단지폐수종말처리장
clearvars discharge_j_yymm
for i = 1:length(indx_date_j)
        discharge_j_yymm(i) = nanmean(merge_j_discharge(indx_date_j{i}));
end

 clearvars tn_j_yymm
for i = 1:length(indx_date_j)
        tn_j_yymm(i) = nanmean(merge_j_tn(indx_date_j{i}));
end
 
clearvars tp_j_yymm
for i = 1:length(indx_date_j)
        tp_j_yymm(i) = nanmean(merge_j_tp(indx_date_j{i}));
end
 
clearvars ss_j_yymm
for i = 1:length(indx_date_j)
        ss_j_yymm(i) = nanmean(merge_j_ss(indx_date_j{i}));
end

%여수하수종말처리장
clearvars discharge_ye_yymm
for i = 1:length(indx_date_ye)
        discharge_ye_yymm(i) = nanmean(merge_ye_discharge(indx_date_ye{i}));
end

 clearvars tn_ye_yymm
for i = 1:length(indx_date_ye)
        tn_ye_yymm(i) = nanmean(merge_ye_tn(indx_date_ye{i}));
end
 
clearvars tp_ye_yymm
for i = 1:length(indx_date_ye)
        tp_ye_yymm(i) = nanmean(merge_ye_tp(indx_date_ye{i}));
end
 
clearvars ss_ye_yymm
for i = 1:length(indx_date_ye)
        ss_ye_yymm(i) = nanmean(merge_ye_ss(indx_date_ye{i}));
end

%진월하수종말처리장
clearvars discharge_jw_yymm
for i = 1:length(indx_date_jw)
        discharge_jw_yymm(i) = nanmean(merge_jw_discharge(indx_date_jw{i}));
end

 clearvars tn_jw_yymm
for i = 1:length(indx_date_jw)
        tn_jw_yymm(i) = nanmean(merge_jw_tn(indx_date_jw{i}));
end
 
clearvars tp_jw_yymm
for i = 1:length(indx_date_jw)
        tp_jw_yymm(i) = nanmean(merge_jw_tp(indx_date_jw{i}));
end
 
clearvars ss_jw_yymm
for i = 1:length(indx_date_jw)
        ss_jw_yymm(i) = nanmean(merge_jw_ss(indx_date_jw{i}));
end

cd D:\장기생태\Dynamic\KOEM
% name_tag{sp_gy}
% 愿묒뼇?빆, 愿묒뼇4, 愿묒뼇3, 愿묒뼇2, 愿묒뼇1, 愿묒뼇5, ?뿬?닔2, ?뿬?닔3, ?뿬?닔1
% [res I]=sort([4,3,2,1,5]);
P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;

yoonja=load('yoonjakangs_koem_data_monthly.mat'); 
no3_sur=[yoonja.no3_sur(1,:) NaN(1,12)]; %1997 to 2020
no3_bot=[yoonja.no3_bot(1,:) NaN(1,12)]; %1997 to 2020
nh4_sur=[yoonja.nh4_sur(1,:) NaN(1,12)]; %1997 to 2020
nh4_bot=[yoonja.nh4_bot(1,:) NaN(1,12)]; %1997 to 2020
po4_sur=[yoonja.po4_sur(1,:) NaN(1,12)]; %1997 to 2020
po4_bot=[yoonja.po4_bot(1,:) NaN(1,12)]; %1997 to 2020

cd D:\장기생태\Dynamic\06_river\하수종말처리장
load sewer_from_LTME.mat  %2007~2017, means 121:252 
% dis_sewer, tn_sewer, tp_sewer

%광양산업단지폐수종말처리장
dis_g_2set = discharge_g_yymm; dis_g_2set(121:252)=dis_sewer(1,:); dis_g_2set(277:288)=NaN;
tn_g_2set = tn_g_yymm; tn_g_2set(121:252)=tn_sewer(1,:); tn_g_2set(277:288)=NaN;
tp_g_2set = tp_g_yymm; tp_g_2set(121:252)=tp_sewer(1,:); tp_g_2set(277:288)=NaN;

%광양중앙하수종말처리장
dis_gc_2set = discharge_gc_yymm; dis_gc_2set(121:252)=dis_sewer(2,:); dis_gc_2set(277:288)=NaN;
tn_gc_2set = tn_gc_yymm; tn_gc_2set(121:252)=tn_sewer(2,:); tn_gc_2set(277:288)=NaN;
tp_gc_2set = tp_gc_yymm; tp_gc_2set(121:252)=tp_sewer(2,:); tp_gc_2set(277:288)=NaN;

%광양하수종말처리장
dis_gh_2set = discharge_gh_yymm; dis_gh_2set(121:252)=dis_sewer(3,:); dis_gh_2set(277:288)=NaN;
tn_gh_2set = tn_gh_yymm; tn_gh_2set(121:252)=tn_sewer(3,:); tn_gh_2set(277:288)=NaN;
tp_gh_2set = tp_gh_yymm; tp_gh_2set(121:252)=tp_sewer(3,:); tp_gh_2set(277:288)=NaN;

%광영하수종말처리장
dis_gh2_2set = discharge_gh2_yymm; dis_gh2_2set(121:252)=dis_sewer(4,:); dis_gh2_2set(277:288)=NaN;
tn_gh2_2set = tn_gh2_yymm; tn_gh2_2set(121:252)=tn_sewer(4,:); tn_gh2_2set(277:288)=NaN;
tp_gh2_2set = tp_gh2_yymm; tp_gh2_2set(121:252)=tp_sewer(4,:); tp_gh2_2set(277:288)=NaN;

%여수월내산업단지폐수종말처리장
dis_2set = discharge_yymm; dis_2set(121:252)=dis_sewer(5,:); dis_2set(277:288)=NaN;
tn_2set = tn_yymm; tn_2set(121:252)=tn_sewer(5,:); tn_2set(277:288)=NaN;
tp_2set = tp_yymm; tp_2set(121:252)=tp_sewer(5,:); tp_2set(277:288)=NaN;

%여수율촌산업단지폐수종말처리장
dis_y_2set = discharge_y_yymm; dis_y_2set(121:252)=dis_sewer(6,:); dis_y_2set(277:288)=NaN;
tn_y_2set = tn_y_yymm; tn_y_2set(121:252)=tn_sewer(6,:); tn_y_2set(277:288)=NaN;
tp_y_2set = tp_y_yymm; tp_y_2set(121:252)=tp_sewer(6,:); tp_y_2set(277:288)=NaN;

%여수중흥산업단지폐수종말처리장
dis_j_2set = discharge_j_yymm; dis_j_2set(121:252)=dis_sewer(7,:); dis_j_2set(277:288)=NaN;
tn_j_2set = tn_j_yymm; tn_j_2set(121:252)=tn_sewer(7,:); tn_j_2set(277:288)=NaN;
tp_j_2set = tp_j_yymm; tp_j_2set(121:252)=tp_sewer(7,:); tp_j_2set(277:288)=NaN;

%여수하수종말처리장
dis_ye_2set = discharge_ye_yymm; dis_ye_2set(121:252)=dis_sewer(8,:); dis_ye_2set(277:288)=NaN;
tn_ye_2set = tn_ye_yymm; tn_ye_2set(121:252)=tn_sewer(8,:); tn_ye_2set(277:288)=NaN;
tp_ye_2set = tp_ye_yymm; tp_ye_2set(121:252)=tp_sewer(8,:); tp_ye_2set(277:288)=NaN;

%진월하수종말처리장
dis_jw_2set = discharge_jw_yymm; dis_jw_2set(121:252)=dis_sewer(9,:); dis_jw_2set(277:288)=NaN;
tn_jw_2set = tn_jw_yymm; tn_jw_2set(121:252)=tn_sewer(9,:); tn_jw_2set(277:288)=NaN;
tp_jw_2set = tp_jw_yymm; tp_jw_2set(121:252)=tp_sewer(9,:); tp_jw_2set(277:288)=NaN;

%check same data or not
% figure; hold on; 
% plot(1:276,discharge_yymm,'b'); plot(121:252, dis_sewer(5,:),'r');

%check same data or not
% figure; hold on; 
% plot(1:276,tn_yymm,'b'); plot(121:252, tn_sewer(5,:),'r');

%check same data or not
% figure; hold on; 
% plot(1:276,tp_yymm,'b'); plot(121:252, tp_sewer(5,:),'r');

cd C:\Users\user\Desktop\장기생태_하수_여수월내
% river near yeosu wallne
wallne=load('yeosu_wallne_data.mat');
wallne_tn=[wallne.tn./1000.*14; NaN]; %%  mM/m^3 to mg/L ;
wallne_tp=[wallne.tp./1000.*30.973762; NaN];

cd D:\장기생태\Dynamic\06_river\하수종말처리장\여수_중흥
% river near yeosu wallne
clearvars jungheung*
jungheung=load('yeosu_jungheung_data.mat');
jungheung_tn=[NaN(length(wallne_tn)-length(jungheung.tn),1); jungheung.tn./1000.*14;]; %%  mM/m^3 to mg/L ;
jungheung_tp=[NaN(length(wallne_tn)-length(jungheung.tn),1); jungheung.tp./1000.*30.973762;];

% dis vs
figure; hold on; 
plot(dis_g_2set,'b'); plot(dis_gc_2set,'g'); plot(dis_gh_2set,'k');
plot(dis_gh2_2set,'color',[148/255 0 211/255]); plot(dis_2set,'r'); plot(dis_y_2set,'c');
 plot(dis_j_2set,'m');  plot(dis_ye_2set,'color',[255/255 192/255 203/255]);  plot(dis_jw_2set,'color',[255/255 165/255 0]);
xticklabels(1997:2020);
xticks(1:12:288);
xlim([121 288]); grid on;
xtickangle(45);
ylabel('kg/m^3')
legend(name_tag);
set(gca,'fontsize',12);

figure; hold on; 
plot(dis_g_2set,'b'); plot(dis_gh2_2set,'color',[148/255 0 211/255]); plot(dis_y_2set,'c');  
% plot(dis_ye_2set,'color',[255/255 192/255 203/255]); 
plot(dis_jw_2set,'color',[255/255 165/255 0]);
xticklabels(1997:2020);
xticks(1:12:288);
xlim([121 288]); grid on;
xtickangle(45);
ylabel('kg/m^3')
legend('광양산업단지폐수종말처리장','광영하수종말처리장','여수율촌산업단지폐수종말처리장','진월하수종말처리장');
set(gca,'fontsize',12);


% dis_2set
plot(wallne_tp,dis_2set,'.')
plot(wallne_tn,dis_2set,'.')

plot(nh4_sur,tn_2set)

% po4 vs
figure; hold on; 
plot(po4_sur./1000,tp_2set,'.'); grid on;
xlabel('KOEM PO4-P (mg/L)'); ylabel('월내산업단지(폐수처리장) TP (mg/L)'); 
set(gca,'fontsize',13);

%
clearvars nonanidx
tp_2set=tp_2set';
nonanidx=find(isnan(wallne_tp + tp_2set)==0);
p_w=polyfit(wallne_tp(nonanidx),tp_2set(nonanidx),1);
figure; hold on; 
plot(wallne_tp,tp_2set,'.'); grid on;
xlabel('여수-월내하천 TP (mg/L)'); ylabel('월내산업단지(폐수처리장) TP (mg/L)'); 
set(gca,'fontsize',13);
plot(wallne_tp,polyval(p_w,wallne_tp),'r');
xlim([0 4])
ylim([0 4])
text(2.5,2,['y = ',num2str(p_w(1),'%0.2f'),'x +',num2str(p_w(2),'%0.2f')],'Color','r')

%
clearvars nonanidx
tp_j_2set=tp_j_2set';
nonanidx=find(isnan(jungheung_tp + tp_j_2set)==0);
p_j=polyfit(jungheung_tp(nonanidx),tp_j_2set(nonanidx),1);
figure; hold on; 
plot(jungheung_tp,tp_j_2set,'.'); grid on;
xlabel('여수-중흥하천 TP (mg/L)'); ylabel('중흥산업단지(폐수처리장) TP (mg/L)'); 
set(gca,'fontsize',13);
xlim([0 4])
ylim([0 4])
plot(jungheung_tp,polyval(p_j,jungheung_tp),'r');
text(2.5,2,['y = ',num2str(p_j(1),'%0.2f'),'x +',num2str(p_j(2),'%0.2f')],'Color','r')


% TN vs
figure; hold on; 
plot(nh4_sur./1000,tn_2set,'.'); grid on;
xlabel('KOEM NH4-N (mg/L)'); ylabel('월내산업단지(폐수처리장) TP (mg/L)'); 
set(gca,'fontsize',13);

%
clearvars nonanidx
tn_2set=tn_2set';
nonanidx=find(isnan(wallne_tn + tn_2set)==0);
p_w_n=polyfit(wallne_tn(nonanidx),tn_2set(nonanidx),1);
figure; hold on; 
plot(wallne_tn,tn_2set,'.'); grid on;
xlabel('여수-월내하천 TN (mg/L)'); ylabel('월내산업단지(폐수처리장) TN (mg/L)'); 
set(gca,'fontsize',13);
xlim([0 50])
ylim([0 50])
plot(wallne_tn,polyval(p_w_n,wallne_tn),'r');
text(30,10,['y = ',num2str(p_w_n(1),'%0.2f'),'x +',num2str(p_w_n(2),'%0.2f')],'Color','r')

%
clearvars nonanidx
tn_j_2set=tn_j_2set';
nonanidx=find(isnan(jungheung_tn + tn_j_2set)==0);
p_j_n=polyfit(jungheung_tn(nonanidx),tn_j_2set(nonanidx),1);
figure; hold on; 
plot(jungheung_tn,tn_j_2set,'.'); grid on;
xlabel('여수-중흥하천 TN (mg/L)'); ylabel('중흥산업단지(폐수처리장) TP (mg/L)'); 
set(gca,'fontsize',13);
xlim([0 50])
ylim([0 50])
plot(jungheung_tn,polyval(p_j_n,jungheung_tn),'r');
text(30,10,['y = ',num2str(p_j_n(1),'%0.2f'),'x +',num2str(p_j_n(2),'%0.2f')],'Color','r')

% 2nd
clearvars nonanidx p_j2
% tn_j_2set=tn_j_2set';
jungheung_tn2=jungheung_tn;
jungheung_tn2(199:end)=NaN;
nonanidx=find(isnan(jungheung_tn2 + tn_j_2set)==0);
p_j_n2=polyfit(jungheung_tn2(nonanidx),tn_j_2set(nonanidx),1);
figure; hold on; 
plot(jungheung_tn2,tn_j_2set,'.'); grid on;
xlabel('여수-중흥하천 TN (mg/L)'); ylabel('중흥산업단지(폐수처리장) TN (mg/L)'); 
set(gca,'fontsize',13);
xlim([0 50])
ylim([0 50])
plot(jungheung_tn,polyval(p_j_n2,jungheung_tn),'r');
text(30,10,['y = ',num2str(p_j_n2(1),'%0.2f'),'x +',num2str(p_j_n2(2),'%0.2f')],'Color','r')


corr(wallne_tn,tp_2set)


figure; hold on; 
plot(tn_g_2set,'b'); plot(tn_gc_2set,'g'); plot(tn_gh_2set,'k');
plot(tn_gh2_2set,'color',[148/255 0 211/255]); plot(tn_2set,'r'); plot(tn_y_2set,'c');
 plot(tn_j_2set,'m');  plot(tn_ye_2set,'color',[255/255 192/255 203/255]);  plot(tn_jw_2set,'color',[255/255 165/255 0]);
xticklabels(1997:2020);
xticks(1:12:288);
xlim([121 288]); grid on;
xtickangle(45);
ylabel('TN (mg/L)')
legend(name_tag);
set(gca,'fontsize',12);

figure; hold on; 
plot(tp_g_2set,'b'); plot(tp_gc_2set,'g'); plot(tp_gh_2set,'k');
plot(tp_gh2_2set,'color',[148/255 0 211/255]); plot(tp_2set,'r'); % plot(tp_y_2set,'c');
 plot(tp_j_2set,'m');  plot(tp_ye_2set,'color',[255/255 192/255 203/255]);  plot(tp_jw_2set,'color',[255/255 165/255 0]);
xticklabels(1997:2020);
xticks(1:12:288);
xlim([121 288]); grid on;
xtickangle(45);
ylabel('TP (mg/L)')
% legend(name_tag);
set(gca,'fontsize',12);

figure; hold on; 
plot(wallne_tp); plot(tp_2set,'r'); %plot(po4_sur./1000, 'g+');
plot(polyval(p_w,wallne_tp),'g');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TP (mg/L)')
legend('여수-월내하천','월내산업단지','recon.');
set(gca,'fontsize',12);

figure; hold on; 
plot(wallne_tn); plot(tn_2set,'r'); %plot(po4_sur./1000, 'g+');
plot(polyval(p_w_n,wallne_tn),'g');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TN (mg/L)')
legend('여수-월내하천','월내산업단지','recon.');
set(gca,'fontsize',12);


figure; hold on; 
plot(jungheung_tp); plot(tp_j_2set,'r'); %plot(po4_sur./1000, 'g+');
plot(polyval(p_j,jungheung_tp),'g'); 
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TP (mg/L)')
legend('여수-중흥하천','중흥산업단지','recon.');
set(gca,'fontsize',12);

figure; hold on; 
plot(jungheung_tn); plot(tn_j_2set,'r'); %plot(po4_sur./1000, 'g+');
plot(polyval(p_j_n,jungheung_tn),'g'); plot(polyval(p_j_n2,jungheung_tn2),'color',[255/255 192/255 203/255]);
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TN (mg/L)')
legend('여수-중흥하천','중흥산업단지','recon.','recon.(~13)');
set(gca,'fontsize',12);

%% RECON. PLT
clearvars tp_2set_re
tp_2set_re = tp_2set;
wallne_tp_recon = polyval(p_w,wallne_tp);
tp_2set_re(isnan(tp_2set))=wallne_tp_recon(isnan(tp_2set));

figure; hold on; 
plot(wallne_tp);  %plot(po4_sur./1000, 'g+');
plot(tp_2set_re,'g');plot(tp_2set,'r');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TP (mg/L)')
legend('여수-월내하천','recon.+ 월내산업단지','월내산업단지');
set(gca,'fontsize',12);

clearvars tn_2set_re
tn_2set_re = tn_2set;
wallne_tn_recon = polyval(p_w_n,wallne_tn);
tn_2set_re(isnan(tn_2set))=wallne_tn_recon(isnan(tn_2set));

figure; hold on; 
plot(wallne_tn);  %plot(po4_sur./1000, 'g+');
plot(tn_2set_re,'g'); plot(tn_2set,'r');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TN (mg/L)')
legend('여수-월내하천','recon.+ 월내산업단지','월내산업단지');
set(gca,'fontsize',12);


clearvars tp_j_2set_re
tp_j_2set_re = tp_j_2set;
jungheung_tp_recon = polyval(p_j,jungheung_tp);
tp_j_2set_re(isnan(tp_j_2set))=jungheung_tp_recon(isnan(tp_j_2set));

figure; hold on; 
plot(jungheung_tp);  %plot(po4_sur./1000, 'g+');
plot(tp_j_2set_re,'g'); plot(tp_j_2set,'r');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TP (mg/L)')
legend('여수-중흥하천','recon.+ 중흥산업단지','중흥산업단지');
set(gca,'fontsize',12);

clearvars tn_j_2set_re
tn_j_2set_re = tn_j_2set;
jungheung_tn_recon = polyval(p_j,jungheung_tn);
tn_j_2set_re(isnan(tn_j_2set))=jungheung_tn_recon(isnan(tn_j_2set));

figure; hold on; 
plot(jungheung_tn); %plot(po4_sur./1000, 'g+');
plot(tn_j_2set_re,'g'); plot(tn_j_2set,'r');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TN (mg/L)')
legend('여수-중흥하천','recon.+ 중흥산업단지','중흥산업단지');
set(gca,'fontsize',12);

%% RECON. PLT 2
clearvars tp_2set_re2
tp_2set_re2 = tp_2set;
tp_2set_re2(isnan(tp_2set))=wallne_tp(isnan(tp_2set));

figure; hold on; 
plot(wallne_tp); 
plot(tp_2set_re2,'g');plot(tp_2set,'r');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TP (mg/L)')
legend('여수-월내하천','recon.+ 월내산업단지','월내산업단지');
set(gca,'fontsize',12);

clearvars tn_2set_re2
tn_2set_re2 = tn_2set;
tn_2set_re2(isnan(tn_2set))=wallne_tn(isnan(tn_2set));

figure; hold on; 
plot(wallne_tn);  %plot(po4_sur./1000, 'g+');
plot(tn_2set_re2,'g'); plot(tn_2set,'r');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TN (mg/L)')
legend('여수-월내하천','recon.+ 월내산업단지','월내산업단지');
set(gca,'fontsize',12);


clearvars tp_j_2set_re2
tp_j_2set_re2 = tp_j_2set;
tp_j_2set_re2(isnan(tp_j_2set))=jungheung_tp(isnan(tp_j_2set));

figure; hold on; 
plot(jungheung_tp);  %plot(po4_sur./1000, 'g+');
plot(tp_j_2set_re2,'g'); plot(tp_j_2set,'r');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TP (mg/L)')
legend('여수-중흥하천','recon.+ 중흥산업단지','중흥산업단지');
set(gca,'fontsize',12);

clearvars tn_j_2set_re2
tn_j_2set_re2 = tn_j_2set;
tn_j_2set_re2(isnan(tn_j_2set))=jungheung_tn(isnan(tn_j_2set));
tn_j_2set_re2(277:end)=NaN;

figure; hold on; 
plot(jungheung_tn); %plot(po4_sur./1000, 'g+');
plot(tn_j_2set_re2,'g'); plot(tn_j_2set,'r');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
xtickangle(45);
ylabel('TN (mg/L)')
legend('여수-중흥하천','recon.+ 중흥산업단지','중흥산업단지');
set(gca,'fontsize',12);

save('remake_sewer_data_monthly_to06to15.mat');
% tp_2set_re2, tn_2set_re2
% tp_j_2set_re2, tn_j_2set_re2

%% regime shift test

tp_w_nonan =find(isnan(tp_2set_re2)==0);
tp_2set_re2(tp_w_nonan)
ref_date(tp_w_nonan)

tn_w_nonan =find(isnan(tn_2set_re2)==0);
tn_2set_re2(tn_w_nonan)

%
tp_j_nonan =find(isnan(tp_j_2set_re2)==0);
tp_j_2set_re2(tp_j_nonan)
ref_date(tp_j_nonan)

tn_j_nonan =find(isnan(tn_j_2set_re2)==0);
tn_j_2set_re2(tn_j_nonan)
ref_date(tn_j_nonan)

%% regime climate
% tp_2set_re2, tn_2set_re2
% tp_j_2set_re2, tn_j_2set_re2

% tp_wallne = 1999-12, 2008-10
ttp_filter=NaN(length(tp_2set_re2),1);
ttp2_filter=NaN(length(tp_2set_re2),1);
ttp3_filter=NaN(length(tp_2set_re2),1);

find_ttp_wallne = {'2006-12'};
find_ttp2_wallne = {'2015-12'};
find_ttp_date=find(strcmp(ref_date,find_ttp_wallne) == 1);
find_ttp2_date=find(strcmp(ref_date,find_ttp2_wallne) == 1);

ttp_filter(1:find_ttp_date) = 1; ttp_filter(find_ttp_date+1:end) = NaN; % pass to 1999-11 (1st)
ttp2_filter(find_ttp_date+1:find_ttp2_date) = 1; ttp2_filter(find_ttp2_date+1:end) = NaN; % pass from 1999-12 to 2008-09 (2nd)
ttp3_filter(1:find_ttp2_date) = NaN; ttp3_filter(find_ttp2_date+1:end) = 1; % pass from 2008-10 to end (3rd)

% tn_wallne = 1999-10, 2003-11
ttn_filter=NaN(length(tn_2set_re2),1);
ttn2_filter=NaN(length(tn_2set_re2),1);
ttn3_filter=NaN(length(tn_2set_re2),1);

find_ttn_wallne = {'2006-12'};
find_ttn2_wallne = {'2015-12'};
find_ttn_date=find(strcmp(ref_date,find_ttn_wallne) == 1);
find_ttn2_date=find(strcmp(ref_date,find_ttn2_wallne) == 1);

ttn_filter(1:find_ttn_date) = 1; ttn_filter(find_ttn_date+1:end) = NaN; % pass to 1999-09 (1st)
ttn2_filter(find_ttn_date+1:find_ttn2_date) = 1; ttn2_filter(find_ttn2_date+1:end) = NaN; % pass from 1999-10 to 2003-09 (2nd)
ttn3_filter(1:find_ttn2_date) = NaN; ttn3_filter(find_ttn2_date+1:end) = 1; % pass from 2003-10 to end (3rd)


% tp_jungheung = 2003-05, 2009-01
ttp_j_filter=NaN(length(tn_j_2set_re2),1);
ttp2_j_filter=NaN(length(tn_j_2set_re2),1);
ttp3_j_filter=NaN(length(tn_j_2set_re2),1);

find_ttp_jungheung = {'2006-12'};
find_ttp2_jungheung = {'2015-12'};
find_ttp_date=find(strcmp(ref_date,find_ttp_jungheung) == 1);
find_ttp2_date=find(strcmp(ref_date,find_ttp2_jungheung) == 1);

ttp_j_filter(1:find_ttp_date) = 1; ttp_j_filter(find_ttp_date+1:end) = NaN; % pass to 2003-04 (1st)
ttp2_j_filter(find_ttp_date+1:find_ttp2_date) = 1; ttp2_j_filter(find_ttp2_date+1:end) = NaN; % pass from 2003-05 to 2008-12 (2nd)
ttp3_j_filter(1:find_ttp2_date) = NaN; ttp3_j_filter(find_ttp2_date+1:end) = 1; % pass from 2009-01 to end (3rd)


% tn_jungheung = 2003-05, 2009-01
ttn_j_filter=NaN(length(tn_j_2set_re2),1);
ttn2_j_filter=NaN(length(tn_j_2set_re2),1);
ttn3_j_filter=NaN(length(tn_j_2set_re2),1);

find_ttn_jungheung = {'2006-12'};
find_ttn2_jungheung = {'2015-12'};
find_ttn_date=find(strcmp(ref_date,find_ttn_jungheung) == 1);
find_ttn2_date=find(strcmp(ref_date,find_ttn2_jungheung) == 1);

ttn_j_filter(1:find_ttn_date) = 1; ttn_j_filter(find_ttn_date+1:end) = NaN; % pass to 2003-04 (1st)
ttn2_j_filter(find_ttn_date+1:find_ttn2_date) = 1; ttn2_j_filter(find_ttn2_date+1:end) = NaN; % pass from 2003-05 to 2008-12 (2nd)
ttn3_j_filter(1:find_ttn2_date) = NaN; ttn3_j_filter(find_ttn2_date+1:end) = 1; % pass from 2009-01 to end (3rd)


% regime
% tp_2set_re2, tn_2set_re2
% tp_j_2set_re2, tn_j_2set_re2
clearvars regm_*
regm_tp_w_1st = tp_2set_re2 .* ttp_filter; 
regm_tp_w_2nd = tp_2set_re2 .* ttp2_filter; 
regm_tp_w_3rd  = tp_2set_re2 .* ttp3_filter;

regm_tn_w_1st = tn_2set_re2 .* ttn_filter; 
regm_tn_w_2nd = tn_2set_re2 .* ttn2_filter; 
regm_tn_w_3rd  = tn_2set_re2 .* ttn3_filter;

regm_tp_j_1st = tp_j_2set_re2 .* ttp_j_filter; 
regm_tp_j_2nd = tp_j_2set_re2 .* ttp2_j_filter; 
regm_tp_j_3rd = tp_j_2set_re2 .* ttp3_j_filter; 

regm_tn_j_1st = tn_j_2set_re2 .* ttn_j_filter; 
regm_tn_j_2nd = tn_j_2set_re2 .* ttn2_j_filter; 
regm_tn_j_3rd  = tn_j_2set_re2 .* ttn3_j_filter;

%spurious data extract on 2010
dis_g_2set(157:168)=NaN;
dis_gc_2set(157:168)=NaN;
dis_gh_2set(157:168)=NaN;
dis_gh2_2set(157:168)=NaN;
dis_2set(157:168)=NaN;
dis_y_2set(157:168)=NaN;
dis_j_2set(157:168)=NaN;
dis_ye_2set(157:168)=NaN;
dis_jw_2set(157:168)=NaN;

for i=1:3 %sigma
regm_tn_w_1st_3s = regm_tn_w_1st;
regm_tn_w_1st_3s(regm_tn_w_1st_3s > nanmean(regm_tn_w_1st) + i*nanstd(regm_tn_w_1st)) =NaN;
regm_tn_w_1st_3s(regm_tn_w_1st_3s < nanmean(regm_tn_w_1st) - i*nanstd(regm_tn_w_1st)) =NaN;
regm_tn_w_1st_ns(i,:)= regm_tn_w_1st_3s;

regm_tn_w_2nd_3s = regm_tn_w_2nd;
regm_tn_w_2nd_3s(regm_tn_w_2nd_3s > nanmean(regm_tn_w_2nd) + i*nanstd(regm_tn_w_2nd)) =NaN;
regm_tn_w_2nd_3s(regm_tn_w_2nd_3s < nanmean(regm_tn_w_2nd) - i*nanstd(regm_tn_w_2nd)) =NaN;
regm_tn_w_2nd_ns(i,:)= regm_tn_w_2nd_3s;

regm_tn_w_3rd_3s = regm_tn_w_3rd;
regm_tn_w_3rd_3s(regm_tn_w_3rd_3s > nanmean(regm_tn_w_3rd) + i*nanstd(regm_tn_w_3rd)) =NaN;
regm_tn_w_3rd_3s(regm_tn_w_3rd_3s < nanmean(regm_tn_w_3rd) - i*nanstd(regm_tn_w_3rd)) =NaN;
regm_tn_w_3rd_ns(i,:)= regm_tn_w_3rd_3s;

regm_tp_w_1st_3s = regm_tp_w_1st;
regm_tp_w_1st_3s(regm_tp_w_1st_3s > nanmean(regm_tp_w_1st) + i*nanstd(regm_tp_w_1st)) =NaN;
regm_tp_w_1st_3s(regm_tp_w_1st_3s < nanmean(regm_tp_w_1st) - i*nanstd(regm_tp_w_1st)) =NaN;
regm_tp_w_1st_ns(i,:)= regm_tp_w_1st_3s;

regm_tp_w_2nd_3s = regm_tp_w_2nd;
regm_tp_w_2nd_3s(regm_tp_w_2nd_3s > nanmean(regm_tp_w_2nd) + i*nanstd(regm_tp_w_2nd)) =NaN;
regm_tp_w_2nd_3s(regm_tp_w_2nd_3s < nanmean(regm_tp_w_2nd) - i*nanstd(regm_tp_w_2nd)) =NaN;
regm_tp_w_2nd_ns(i,:)= regm_tp_w_2nd_3s;

regm_tp_w_3rd_3s = regm_tp_w_3rd;
regm_tp_w_3rd_3s(regm_tp_w_3rd_3s > nanmean(regm_tp_w_3rd) + i*nanstd(regm_tp_w_3rd)) =NaN;
regm_tp_w_3rd_3s(regm_tp_w_3rd_3s < nanmean(regm_tp_w_3rd) - i*nanstd(regm_tp_w_3rd)) =NaN;
regm_tp_w_3rd_ns(i,:)= regm_tp_w_3rd_3s;

regm_tn_j_1st_3s = regm_tn_j_1st;
regm_tn_j_1st_3s(regm_tn_j_1st_3s > nanmean(regm_tn_j_1st) + i*nanstd(regm_tn_j_1st)) =NaN;
regm_tn_j_1st_3s(regm_tn_j_1st_3s < nanmean(regm_tn_j_1st) - i*nanstd(regm_tn_j_1st)) =NaN;
regm_tn_j_1st_ns(i,:)= regm_tn_j_1st_3s;

regm_tn_j_2nd_3s = regm_tn_j_2nd;
regm_tn_j_2nd_3s(regm_tn_j_2nd_3s > nanmean(regm_tn_j_2nd) + i*nanstd(regm_tn_j_2nd)) =NaN;
regm_tn_j_2nd_3s(regm_tn_j_2nd_3s < nanmean(regm_tn_j_2nd) - i*nanstd(regm_tn_j_2nd)) =NaN;
regm_tn_j_2nd_ns(i,:)= regm_tn_j_2nd_3s;

regm_tn_j_3rd_3s = regm_tn_j_3rd;
regm_tn_j_3rd_3s(regm_tn_j_3rd_3s > nanmean(regm_tn_j_3rd) + i*nanstd(regm_tn_j_3rd)) =NaN;
regm_tn_j_3rd_3s(regm_tn_j_3rd_3s < nanmean(regm_tn_j_3rd) - i*nanstd(regm_tn_j_3rd)) =NaN;
regm_tn_j_3rd_ns(i,:)= regm_tn_j_3rd_3s;

regm_tp_j_1st_3s = regm_tp_j_1st;
regm_tp_j_1st_3s(regm_tp_j_1st_3s > nanmean(regm_tp_j_1st) + i*nanstd(regm_tp_j_1st)) =NaN;
regm_tp_j_1st_3s(regm_tp_j_1st_3s < nanmean(regm_tp_j_1st) - i*nanstd(regm_tp_j_1st)) =NaN;
regm_tp_j_1st_ns(i,:)= regm_tp_j_1st_3s;

regm_tp_j_2nd_3s = regm_tp_j_2nd;
regm_tp_j_2nd_3s(regm_tp_j_2nd_3s > nanmean(regm_tp_j_2nd) + i*nanstd(regm_tp_j_2nd)) =NaN;
regm_tp_j_2nd_3s(regm_tp_j_2nd_3s < nanmean(regm_tp_j_2nd) - i*nanstd(regm_tp_j_2nd)) =NaN;
regm_tp_j_2nd_ns(i,:)= regm_tp_j_2nd_3s;

regm_tp_j_3rd_3s = regm_tp_j_3rd;
regm_tp_j_3rd_3s(regm_tp_j_3rd_3s > nanmean(regm_tp_j_3rd) + i*nanstd(regm_tp_j_3rd)) =NaN;
regm_tp_j_3rd_3s(regm_tp_j_3rd_3s < nanmean(regm_tp_j_3rd) - i*nanstd(regm_tp_j_3rd)) =NaN;
regm_tp_j_3rd_ns(i,:)= regm_tp_j_3rd_3s;

% others 7points (except wallne & jungheung) without regime shift
% 1.광양산업단지폐수종말처리장
clearvars tempd
tempd=tn_g_2set;
tempd(tempd > nanmean(tn_g_2set) + i*nanstd(tn_g_2set)) =NaN;
tempd(tempd < nanmean(tn_g_2set) - i*nanstd(tn_g_2set)) =NaN;
tn_g_ns(i,:)= tempd;

clearvars tempd
tempd=tp_g_2set;
tempd(tempd > nanmean(tp_g_2set) + i*nanstd(tp_g_2set)) =NaN;
tempd(tempd < nanmean(tp_g_2set) - i*nanstd(tp_g_2set)) =NaN;
tp_g_ns(i,:)= tempd;
% 2.광양중앙하수종말처리장
clearvars tempd
tempd=tn_gc_2set;
tempd(tempd > nanmean(tn_gc_2set) + i*nanstd(tn_gc_2set)) =NaN;
tempd(tempd < nanmean(tn_gc_2set) - i*nanstd(tn_gc_2set)) =NaN;
tn_gc_ns(i,:)= tempd;

clearvars tempd
tempd=tp_gc_2set;
tempd(tempd > nanmean(tp_gc_2set) + i*nanstd(tp_gc_2set)) =NaN;
tempd(tempd < nanmean(tp_gc_2set) - i*nanstd(tp_gc_2set)) =NaN;
tp_gc_ns(i,:)= tempd;

% 3.광양하수종말처리장
clearvars tempd
tempd=tn_gh_2set;
tempd(tempd > nanmean(tn_gh_2set) + i*nanstd(tn_gh_2set)) =NaN;
tempd(tempd < nanmean(tn_gh_2set) - i*nanstd(tn_gh_2set)) =NaN;
tn_gh_ns(i,:)= tempd;

clearvars tempd
tempd=tp_gh_2set;
tempd(tempd > nanmean(tp_gh_2set) + i*nanstd(tp_gh_2set)) =NaN;
tempd(tempd < nanmean(tp_gh_2set) - i*nanstd(tp_gh_2set)) =NaN;
tp_gh_ns(i,:)= tempd;

% 4.광영하수종말처리장
clearvars tempd
tempd=tn_gh2_2set;
tempd(tempd > nanmean(tn_gh2_2set) + i*nanstd(tn_gh2_2set)) =NaN;
tempd(tempd < nanmean(tn_gh2_2set) - i*nanstd(tn_gh2_2set)) =NaN;
tn_gh2_ns(i,:)= tempd;

clearvars tempd
tempd=tp_gh2_2set;
tempd(tempd > nanmean(tp_gh2_2set) + i*nanstd(tp_gh2_2set)) =NaN;
tempd(tempd < nanmean(tp_gh2_2set) - i*nanstd(tp_gh2_2set)) =NaN;
tp_gh2_ns(i,:)= tempd;

% 6.여수율촌산업단지폐수종말처리장
clearvars tempd
tempd=tn_y_2set;
tempd(tempd > nanmean(tn_y_2set) + i*nanstd(tn_y_2set)) =NaN;
tempd(tempd < nanmean(tn_y_2set) - i*nanstd(tn_y_2set)) =NaN;
tn_y_ns(i,:)= tempd;

clearvars tempd
tempd=tp_y_2set;
tempd(tempd > nanmean(tp_y_2set) + i*nanstd(tp_y_2set)) =NaN;
tempd(tempd < nanmean(tp_y_2set) - i*nanstd(tp_y_2set)) =NaN;
tp_y_ns(i,:)= tempd;

% 8.여수하수종말처리장
clearvars tempd
tempd=tn_ye_2set;
tempd(tempd > nanmean(tn_ye_2set) + i*nanstd(tn_ye_2set)) =NaN;
tempd(tempd < nanmean(tn_ye_2set) - i*nanstd(tn_ye_2set)) =NaN;
tn_ye_ns(i,:)= tempd;

clearvars tempd
tempd=tp_ye_2set;
tempd(tempd > nanmean(tp_ye_2set) + i*nanstd(tp_ye_2set)) =NaN;
tempd(tempd < nanmean(tp_ye_2set) - i*nanstd(tp_ye_2set)) =NaN;
tp_ye_ns(i,:)= tempd;

% 9.진월하수종말처리장
clearvars tempd
tempd=tn_jw_2set;
tempd(tempd > nanmean(tn_jw_2set) + i*nanstd(tn_jw_2set)) =NaN;
tempd(tempd < nanmean(tn_jw_2set) - i*nanstd(tn_jw_2set)) =NaN;
tn_jw_ns(i,:)= tempd;

clearvars tempd
tempd=tp_jw_2set;
tempd(tempd > nanmean(tp_jw_2set) + i*nanstd(tp_jw_2set)) =NaN;
tempd(tempd < nanmean(tp_jw_2set) - i*nanstd(tp_jw_2set)) =NaN;
tp_jw_ns(i,:)= tempd;
end
% % 5.여수월내산업단지폐수종말처리장
% % 7.여수중흥산업단지폐수종말처리장

clearvars clim_tp_* clim_tn_*
for i = 1:12 
    for j = 1:3 %sigma
    clim_tp_w_1st_ns(j,i)=nanmean(regm_tp_w_1st_ns(j,i:12:end));
    clim_tp_w_2nd_ns(j,i)=nanmean(regm_tp_w_2nd_ns(j,i:12:end));
    clim_tp_w_3rd_ns(j,i)=nanmean(regm_tp_w_3rd_ns(j,i:12:end));

    clim_tn_w_1st_ns(j,i)=nanmean(regm_tn_w_1st_ns(j,i:12:end));
    clim_tn_w_2nd_ns(j,i)=nanmean(regm_tn_w_2nd_ns(j,i:12:end));
    clim_tn_w_3rd_ns(j,i)=nanmean(regm_tn_w_3rd_ns(j,i:12:end));

    clim_tp_j_1st_ns(j,i)=nanmean(regm_tp_j_1st_ns(j,i:12:end));
    clim_tp_j_2nd_ns(j,i)=nanmean(regm_tp_j_2nd_ns(j,i:12:end));
    clim_tp_j_3rd_ns(j,i)=nanmean(regm_tp_j_3rd_ns(j,i:12:end));


    clim_tn_j_1st_ns(j,i)=nanmean(regm_tn_j_1st_ns(j,i:12:end));
    clim_tn_j_2nd_ns(j,i)=nanmean(regm_tn_j_2nd_ns(j,i:12:end));
    clim_tn_j_3rd_ns(j,i)=nanmean(regm_tn_j_3rd_ns(j,i:12:end));

% others 7points (except wallne & jungheung) without regime shift
        % 1.광양산업단지폐수종말처리장
        clim_tn_g_ns(j,i)=nanmean(tn_g_ns(j,i:12:end));
        clim_tp_g_ns(j,i)=nanmean(tp_g_ns(j,i:12:end));
        clim_tn_g_raw(j,i)=nanmean(tn_g_2set(1,i:12:end));
        clim_tp_g_raw(j,i)=nanmean(tp_g_2set(1,i:12:end));
        % 2.광양중앙하수종말처리장
        clim_tn_gc_ns(j,i)=nanmean(tn_gc_ns(j,i:12:end));
        clim_tp_gc_ns(j,i)=nanmean(tp_gc_ns(j,i:12:end));
        clim_tn_gc_raw(j,i)=nanmean(tn_gc_2set(1,i:12:end));
        clim_tp_gc_raw(j,i)=nanmean(tp_gc_2set(1,i:12:end));

        % 3.광양하수종말처리장
        clim_tn_gh_ns(j,i)=nanmean(tn_gh_ns(j,i:12:end));
        clim_tp_gh_ns(j,i)=nanmean(tp_gh_ns(j,i:12:end));
        clim_tn_gh_raw(j,i)=nanmean(tn_gh_2set(1,i:12:end));
        clim_tp_gh_raw(j,i)=nanmean(tp_gh_2set(1,i:12:end));
        
        % 4.광영하수종말처리장
        clim_tn_gh2_ns(j,i)=nanmean(tn_gh2_ns(j,i:12:end));
        clim_tp_gh2_ns(j,i)=nanmean(tp_gh2_ns(j,i:12:end));
        clim_tn_gh2_raw(j,i)=nanmean(tn_gh2_2set(1,i:12:end));
        clim_tp_gh2_raw(j,i)=nanmean(tp_gh2_2set(1,i:12:end));
        
        % 6.여수율촌산업단지폐수종말처리장
        clim_tn_y_ns(j,i)=nanmean(tn_y_ns(j,i:12:end));
        clim_tp_y_ns(j,i)=nanmean(tp_y_ns(j,i:12:end));
        clim_tn_y_raw(j,i)=nanmean(tn_y_2set(1,i:12:end));
        clim_tp_y_raw(j,i)=nanmean(tp_y_2set(1,i:12:end));
        
        % 8.여수하수종말처리장
        clim_tn_ye_ns(j,i)=nanmean(tn_ye_ns(j,i:12:end));
        clim_tp_ye_ns(j,i)=nanmean(tp_ye_ns(j,i:12:end));
        clim_tn_ye_raw(j,i)=nanmean(tn_ye_2set(1,i:12:end));
        clim_tp_ye_raw(j,i)=nanmean(tp_ye_2set(1,i:12:end));
        
        % 9.진월하수종말처리장
        clim_tn_jw_ns(j,i)=nanmean(tn_jw_ns(j,i:12:end));
        clim_tp_jw_ns(j,i)=nanmean(tp_jw_ns(j,i:12:end));     
        clim_tn_jw_raw(j,i)=nanmean(tn_jw_2set(1,i:12:end));
        clim_tp_jw_raw(j,i)=nanmean(tp_jw_2set(1,i:12:end));    
    end
end

for i = 1:12 
clim_tp_w_1st(i)=nanmean(regm_tp_w_1st(i:12:end));
clim_tp_w_2nd(i)=nanmean(regm_tp_w_2nd(i:12:end));
clim_tp_w_3rd(i)=nanmean(regm_tp_w_3rd(i:12:end));

clim_tn_w_1st(i)=nanmean(regm_tn_w_1st(i:12:end));
clim_tn_w_2nd(i)=nanmean(regm_tn_w_2nd(i:12:end));
clim_tn_w_3rd(i)=nanmean(regm_tn_w_3rd(i:12:end));

clim_tp_j_1st(i)=nanmean(regm_tp_j_1st(i:12:end));
clim_tp_j_2nd(i)=nanmean(regm_tp_j_2nd(i:12:end));
clim_tp_j_3rd(i)=nanmean(regm_tp_j_3rd(i:12:end));

clim_tn_j_1st(i)=nanmean(regm_tn_j_1st(i:12:end));
clim_tn_j_2nd(i)=nanmean(regm_tn_j_2nd(i:12:end));
clim_tn_j_3rd(i)=nanmean(regm_tn_j_3rd(i:12:end));

clim_dis_g(i)=nanmean(dis_g_2set(i:12:end));
clim_dis_gc(i)=nanmean(dis_gc_2set(i:12:end));
clim_dis_gh(i)=nanmean(dis_gh_2set(i:12:end));
clim_dis_gh2(i)=nanmean(dis_gh2_2set(i:12:end));
clim_dis(i)=nanmean(dis_2set(i:12:end));
clim_dis_y(i)=nanmean(dis_y_2set(i:12:end));
clim_dis_j(i)=nanmean(dis_j_2set(i:12:end));
clim_dis_ye(i)=nanmean(dis_ye_2set(i:12:end));
clim_dis_jw(i)=nanmean(dis_jw_2set(i:12:end));
end


%% make daily data

% make eom_d
k=0
for i = 1980:1980
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end
t_tick_pre=sum(eom_d,2);

for i = 1:size(eom_d,1)
    for j = 1:size(eom_d,2)
        eom_d_each(i,j) = sum(eom_d(i,1:j));
    end
end

% make monthly mean
clearvars mod_m_*
for i =1:12
    if i ==1 
        %dis
        clim_dis_g_d(1:eom_d_each(1,i))=clim_dis_g(i);
        clim_dis_gc_d(1:eom_d_each(1,i))=clim_dis_gc(i);
        clim_dis_gh_d(1:eom_d_each(1,i))=clim_dis_gh(i);
        clim_dis_gh2_d(1:eom_d_each(1,i))=clim_dis_gh2(i);
        clim_dis_d(1:eom_d_each(1,i))=clim_dis(i);
        clim_dis_y_d(1:eom_d_each(1,i))=clim_dis_y(i);
        clim_dis_j_d(1:eom_d_each(1,i))=clim_dis_j(i);
        clim_dis_ye_d(1:eom_d_each(1,i))=clim_dis_ye(i);
        clim_dis_jw_d(1:eom_d_each(1,i))=clim_dis_jw(i);
        %TN
        clim_tn_g_ns_d(1:eom_d_each(1,i))=clim_tn_g_ns(3,i);
        clim_tn_gc_ns_d(1:eom_d_each(1,i))=clim_tn_gc_ns(3,i);
        clim_tn_gh_ns_d(1:eom_d_each(1,i))=clim_tn_gh_ns(3,i);
        clim_tn_gh2_ns_d(1:eom_d_each(1,i))=clim_tn_gh2_ns(3,i);
        clim_tn_w_1st_ns_d(1:eom_d_each(1,i))=clim_tn_w_1st_ns(3,i);
        clim_tn_w_2nd_ns_d(1:eom_d_each(1,i))=clim_tn_w_2nd_ns(3,i);
        clim_tn_w_3rd_ns_d(1:eom_d_each(1,i))=clim_tn_w_3rd_ns(3,i);
        clim_tn_y_ns_d(1:eom_d_each(1,i))=clim_tn_y_ns(3,i);
        clim_tn_j_1st_ns_d(1:eom_d_each(1,i))=clim_tn_j_1st_ns(3,i);
        clim_tn_j_2nd_ns_d(1:eom_d_each(1,i))=clim_tn_j_2nd_ns(3,i);
        clim_tn_j_3rd_ns_d(1:eom_d_each(1,i))=clim_tn_j_3rd_ns(3,i);
        clim_tn_ye_ns_d(1:eom_d_each(1,i))=clim_tn_ye_ns(3,i);
        clim_tn_jw_ns_d(1:eom_d_each(1,i))=clim_tn_jw_ns(3,i);
        %TP
        clim_tp_g_ns_d(1:eom_d_each(1,i))=clim_tp_g_ns(3,i);
        clim_tp_gc_ns_d(1:eom_d_each(1,i))=clim_tp_gc_ns(3,i);
        clim_tp_gh_ns_d(1:eom_d_each(1,i))=clim_tp_gh_ns(3,i);
        clim_tp_gh2_ns_d(1:eom_d_each(1,i))=clim_tp_gh2_ns(3,i);
        clim_tp_w_1st_ns_d(1:eom_d_each(1,i))=clim_tp_w_1st_ns(3,i);
        clim_tp_w_2nd_ns_d(1:eom_d_each(1,i))=clim_tp_w_2nd_ns(3,i);
        clim_tp_w_3rd_ns_d(1:eom_d_each(1,i))=clim_tp_w_3rd_ns(3,i);
        clim_tp_y_ns_d(1:eom_d_each(1,i))=clim_tp_y_ns(3,i);
        clim_tp_j_1st_ns_d(1:eom_d_each(1,i))=clim_tp_j_1st_ns(3,i);
        clim_tp_j_2nd_ns_d(1:eom_d_each(1,i))=clim_tp_j_2nd_ns(3,i);
        clim_tp_j_3rd_ns_d(1:eom_d_each(1,i))=clim_tp_j_3rd_ns(3,i);
        clim_tp_ye_ns_d(1:eom_d_each(1,i))=clim_tp_ye_ns(3,i);
        clim_tp_jw_ns_d(1:eom_d_each(1,i))=clim_tp_jw_ns(3,i);
    else
        %dis
        clim_dis_g_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_dis_g(i);
        clim_dis_gc_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_dis_gc(i);
        clim_dis_gh_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_dis_gh(i);
        clim_dis_gh2_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_dis_gh2(i);
        clim_dis_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_dis(i);
        clim_dis_y_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_dis_y(i);
        clim_dis_j_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_dis_j(i);
        clim_dis_ye_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_dis_ye(i);
        clim_dis_jw_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_dis_jw(i);
        %TN
        clim_tn_g_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_g_ns(3,i);
        clim_tn_gc_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_gc_ns(3,i);
        clim_tn_gh_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_gh_ns(3,i);
        clim_tn_gh2_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_gh2_ns(3,i);
        clim_tn_w_1st_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_w_1st_ns(3,i);
        clim_tn_w_2nd_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_w_2nd_ns(3,i);
        clim_tn_w_3rd_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_w_3rd_ns(3,i);
        clim_tn_y_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_y_ns(3,i);
        clim_tn_j_1st_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_j_1st_ns(3,i);
        clim_tn_j_2nd_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_j_2nd_ns(3,i);
        clim_tn_j_3rd_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_j_3rd_ns(3,i);
        clim_tn_ye_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_ye_ns(3,i);
        clim_tn_jw_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tn_jw_ns(3,i);
        %TP
        clim_tp_g_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_g_ns(3,i);
        clim_tp_gc_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_gc_ns(3,i);
        clim_tp_gh_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_gh_ns(3,i);
        clim_tp_gh2_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_gh2_ns(3,i);
        clim_tp_w_1st_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_w_1st_ns(3,i);
        clim_tp_w_2nd_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_w_2nd_ns(3,i);
        clim_tp_w_3rd_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_w_3rd_ns(3,i);
        clim_tp_y_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_y_ns(3,i);
        clim_tp_j_1st_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_j_1st_ns(3,i);
        clim_tp_j_2nd_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_j_2nd_ns(3,i);
        clim_tp_j_3rd_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_j_3rd_ns(3,i);
        clim_tp_ye_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_ye_ns(3,i);
        clim_tp_jw_ns_d(eom_d_each(1,i-1)+1:eom_d_each(1,i))=clim_tp_jw_ns(3,i);
    end
end

return

figure; hold on
yyaxis left
plot(1:12,clim_dis); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_w_1st); ylabel('TP(mg/L)'); 
plot(1:12,clim_tp_w_2nd);
plot(1:12,clim_tp_w_3rd);
% plot(1:12,clim_tp_w_1st_ns);
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);

figure; hold on
yyaxis left
plot(1:12,clim_dis); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_w_1st); ylabel('TN(mg/L)'); 
plot(1:12,clim_tn_w_2nd);
plot(1:12,clim_tn_w_3rd);
% plot(1:12,clim_tp_w_1st_ns);
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);


figure; hold on
% yyaxis left
% plot(1:12,clim_dis); ylabel('m^3/day')
% yyaxis right
plot(1:12,clim_tn_w_1st); ylabel('TN(mg/L)'); 
% plot(1:12,clim_tn_w_2nd);
% plot(1:12,clim_tn_w_3rd);
plot(1:12,clim_tn_w_1st_ns(1,:)); plot(1:12,clim_tn_w_1st_ns(2,:)); plot(1:12,clim_tn_w_1st_ns(3,:));
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);


figure; hold on
% yyaxis left
% plot(1:12,clim_dis); ylabel('m^3/day')
% yyaxis right
plot(1:12,clim_tp_j_1st); ylabel('TP(mg/L)'); 
% plot(1:12,clim_tp_j_2nd);
plot(1:12,clim_tp_j_1st_ns(1,:)); plot(1:12,clim_tp_j_1st_ns(2,:)); plot(1:12,clim_tp_j_1st_ns(3,:));
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
legend


figure; hold on
yyaxis left
plot(1:12,clim_dis_j); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_j_1st); ylabel('TN (mg/L)'); 
plot(1:12,clim_tn_j_2nd);
plot(1:12,clim_tn_j_3rd);
% plot(1:12,clim_tn_j_1st_ns(3,:));
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);

figure; hold on
yyaxis left
plot(1:12,clim_dis_j); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_j_1st); ylabel('TP (mg/L)'); 
plot(1:12,clim_tp_j_2nd);
plot(1:12,clim_tp_j_3rd);
% plot(1:12,clim_tp_j_1st_ns(3,:));
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);

% 1.광양산업단지폐수종말처리장
figure; hold on
yyaxis left
plot(1:12,clim_dis_g); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_g_raw);
plot(1:12,clim_tn_g_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tn_g_ns(2,:),'k');
plot(1:12,clim_tn_g_ns(3,:),'g');
% plot(1:12,clim_tn_j_1st_ns(3,:));
title('1.광양산업단지폐수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 17])

figure; hold on
yyaxis left
plot(1:12,clim_dis_g); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_g_raw);
plot(1:12,clim_tp_g_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tp_g_ns(2,:),'k');
plot(1:12,clim_tp_g_ns(3,:),'g');
title('1.광양산업단지폐수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 1.2])

% 2.광양중앙하수종말처리장
figure; hold on
yyaxis left
plot(1:12,clim_dis_gc); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_gc_raw);
plot(1:12,clim_tn_gc_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tn_gc_ns(2,:),'k');
plot(1:12,clim_tn_gc_ns(3,:),'g');
% plot(1:12,clim_tn_j_1st_ns(3,:));
title('2.광양중앙하수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 17])

figure; hold on
yyaxis left
plot(1:12,clim_dis_gc); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_gc_raw);
plot(1:12,clim_tp_gc_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tp_gc_ns(2,:),'k');
plot(1:12,clim_tp_gc_ns(3,:),'g');
title('2.광양중앙하수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 1.2])

% 3.광양하수종말처리장
figure; hold on
yyaxis left
plot(1:12,clim_dis_gh); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_gh_raw);
plot(1:12,clim_tn_gh_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tn_gh_ns(2,:),'k');
plot(1:12,clim_tn_gh_ns(3,:),'g');
% plot(1:12,clim_tn_j_1st_ns(3,:));
title('3.광양하수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 17])

figure; hold on
yyaxis left
plot(1:12,clim_dis_gh); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_gh_raw);
plot(1:12,clim_tp_gh_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tp_gh_ns(2,:),'k');
plot(1:12,clim_tp_gh_ns(3,:),'g');
title('3.광양하수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 1.2])

% 4.광영하수종말처리장
figure; hold on
yyaxis left
plot(1:12,clim_dis_gh2); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_gh2_raw);
plot(1:12,clim_tn_gh2_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tn_gh2_ns(2,:),'k');
plot(1:12,clim_tn_gh2_ns(3,:),'g');
% plot(1:12,clim_tn_j_1st_ns(3,:));
title('4.광영하수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 17])

figure; hold on
yyaxis left
plot(1:12,clim_dis_gh2); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_gh2_raw);
plot(1:12,clim_tp_gh2_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tp_gh2_ns(2,:),'k');
plot(1:12,clim_tp_gh2_ns(3,:),'g');
title('4.광영하수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 1.2])

% 6.여수율촌산업단지폐수종말처리장
figure; hold on
yyaxis left
plot(1:12,clim_dis_y); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_y_raw);
plot(1:12,clim_tn_y_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tn_y_ns(2,:),'k');
plot(1:12,clim_tn_y_ns(3,:),'g');
% plot(1:12,clim_tn_j_1st_ns(3,:));
title('6.여수율촌산업단지폐수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 17])

figure; hold on
yyaxis left
plot(1:12,clim_dis_y); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_y_raw);
plot(1:12,clim_tp_y_ns(1,:),'m'); ylabel('TP (mg/L)'); 
plot(1:12,clim_tp_y_ns(2,:),'k');
plot(1:12,clim_tp_y_ns(3,:),'g');
title('6.여수율촌산업단지폐수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
% ylim([0 1.2])
ylim([0 inf])

% 8.여수하수종말처리장
figure; hold on
yyaxis left
plot(1:12,clim_dis_ye); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_ye_raw);
plot(1:12,clim_tn_ye_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tn_ye_ns(2,:),'k');
plot(1:12,clim_tn_ye_ns(3,:),'g');
% plot(1:12,clim_tn_j_1st_ns(3,:));
title('8.여수하수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 17])

figure; hold on
yyaxis left
plot(1:12,clim_dis_ye); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_ye_raw);
plot(1:12,clim_tp_ye_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tp_ye_ns(2,:),'k');
plot(1:12,clim_tp_ye_ns(3,:),'g');
title('8.여수하수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 1.2])

% 9.진월하수종말처리장
figure; hold on
yyaxis left
plot(1:12,clim_dis_jw); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_jw_raw);
plot(1:12,clim_tn_jw_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tn_jw_ns(2,:),'k');
plot(1:12,clim_tn_jw_ns(3,:),'g');
% plot(1:12,clim_tn_j_1st_ns(3,:));
title('9.진월하수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 17])

figure; hold on
yyaxis left
plot(1:12,clim_dis_jw); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_jw_raw);
plot(1:12,clim_tp_jw_ns(1,:),'m'); ylabel('TN (mg/L)'); 
plot(1:12,clim_tp_jw_ns(2,:),'k');
plot(1:12,clim_tp_jw_ns(3,:),'g');
title('9.진월하수종말처리장');
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);
ylim([0 1.2])

% % 5.여수월내산업단지폐수종말처리장
% % 7.여수중흥산업단지폐수종말처리장


figure; hold on; grid on;
plot(clim_tp_w_1st); plot(clim_tp_w_2nd,'r');  plot(clim_tp_w_3rd,'g');
xlim([1 12]); xticks(1:12);
ylabel('TP (mg/L)')
title('여천월내 regime climate')
legend('1st','2nd','3rd');
xlabel('month')
set(gca,'fontsize',12);
ylim([0 6]) 

figure; hold on; grid on;
plot(clim_tn_w_1st); plot(clim_tn_w_2nd,'r');  plot(clim_tn_w_3rd,'g');
xlim([1 12]); xticks(1:12);
ylabel('TN (mg/L)')
title('여천월내 regime climate')
legend('1st','2nd','3rd');
xlabel('month')
set(gca,'fontsize',12);
ylim([0 200]) 

figure; hold on; grid on;
plot(clim_tp_j_1st); plot(clim_tp_j_2nd,'r'); % plot(clim_tp_j_3rd,'g');
xlim([1 12]); xticks(1:12);
ylabel('TP (mg/L)')
title('여천중흥 regime climate')
legend('1st','2nd','3rd');
xlabel('month')
set(gca,'fontsize',12);
ylim([0 6]) 


figure; hold on; grid on;
plot(clim_tn_j_1st); plot(clim_tn_j_2nd,'r');  plot(clim_tn_j_3rd,'g');
xlim([1 12]); xticks(1:12);
ylabel('TN (mg/L)')
title('여천중흥 regime climate')
legend('1st','2nd','3rd');
xlabel('month')
set(gca,'fontsize',12);
ylim([0 200]) 

figure; hold on; 
plot(clim_dis_g,'b'); plot(clim_dis_gc,'g'); plot(clim_dis_gh,'k');
plot(clim_dis_gh2,'color',[148/255 0 211/255]); plot(clim_dis,'r'); plot(clim_dis_y,'c');
 plot(clim_dis_j,'m');  plot(clim_dis_ye,'color',[255/255 192/255 203/255]);  plot(clim_dis_jw,'color',[255/255 165/255 0]);
xlabel('month')
xticks(1:12);
xlim([1 12]); grid on;
ylabel('m^3/day')
legend(name_tag);
set(gca,'fontsize',12)

%% nut. * dis = nut. flux
figure; hold on; 
plot(clim_dis .* clim_tn_w_1st ./10^6,'r');  plot(clim_dis_j .* clim_tn_j_1st ./10^6,'m'); %m^3/day * g/m^3 = g/day
plot(clim_dis .* clim_tn_w_2nd ./10^6,'b--');  plot(clim_dis_j .* clim_tn_j_2nd ./10^6,'g--');
plot(clim_dis .* clim_tn_w_3rd ./10^6,'k-.');  plot(clim_dis_j .* clim_tn_j_3rd ./10^6,'y-.');
xlabel('month')
xticks(1:12);
xlim([1 12]); grid on;
ylim([0 10])
ylabel('ton/day')
set(gca,'fontsize',12)

figure; hold on; 
plot(clim_dis .* clim_tp_w_1st ./10^6,'r');  plot(clim_dis_j .* clim_tp_j_1st ./10^6,'m'); %m^3/day * g/m^3 = g/day
plot(clim_dis .* clim_tp_w_2nd ./10^6,'b--');  plot(clim_dis_j .* clim_tp_j_2nd ./10^6,'g--');
plot(clim_dis .* clim_tp_w_3rd ./10^6,'k-.');  %plot(clim_dis_j .* clim_tp_j_3rd ./10^6,'y-.');
xlabel('month')
xticks(1:12);
xlim([1 12]); grid on;
ylim([0 .6])
ylabel('ton/day')
set(gca,'fontsize',12)

mean(clim_tn_w_1st(1:6))

figure; hold on
yyaxis left
plot(1:12,clim_dis); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_w_1st); ylabel('mg/L'); 
plot(1:12,clim_tn_w_2nd);
plot(1:12,clim_tn_w_3rd);
% plot(1:12,smoothdata(clim_tn_w_1st,'loess',6));
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);

figure; hold on
yyaxis left
plot(1:12,clim_dis_j); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tn_j_1st); ylabel('mg/L'); 
plot(1:12,clim_tn_j_2nd);
plot(1:12,clim_tn_j_3rd);
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);

%%
figure; hold on
yyaxis left
plot(1:12,clim_dis); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_w_1st); ylabel('mg/L'); 
plot(1:12,clim_tp_w_2nd);
plot(1:12,clim_tp_w_3rd);
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);

figure; hold on
yyaxis left
plot(1:12,clim_dis_j); ylabel('m^3/day')
yyaxis right
plot(1:12,clim_tp_j_1st); ylabel('mg/L'); 
plot(1:12,clim_tp_j_2nd);
% plot(1:12,clim_tp_j_3rd);
set(gca,'fontsize',12);
grid on; xlim([1 12]); xlabel('month'); xticks(1:12);

save('sewer_monthly_climate_fix_to06to15.mat','clim*');



% load sewer_monthly_climate.mat

%% all st, plot
figure; hold on; 
plot(clim_tn_g_ns(3,:),'b'); plot(clim_tn_gc_ns(3,:),'g'); plot(clim_tn_gh_ns(3,:),'k');
plot(clim_tn_gh2_ns(3,:),'color',[148/255 0 211/255]); plot(clim_tn_w_1st_ns(3,:),'r'); plot(clim_tn_y_ns(3,:),'c');
 plot(clim_tn_j_1st_ns(3,:),'m');  plot(clim_tn_ye_ns(3,:),'color',[255/255 192/255 203/255]);  plot(clim_tn_jw_ns(3,:),'color',[255/255 165/255 0]);
% xticklabels(1997:2020);
xticks(1:12);
xlim([1 12]); grid on;
xtickangle(45);
ylabel('TN (mg/L)')
legend(name_tag);
set(gca,'fontsize',12);
% ylim([0 200]) % y axis match with raw

figure; hold on; 
plot(clim_tn_g_ns(3,:),'b'); plot(clim_tn_gc_ns(3,:),'g'); plot(clim_tn_gh_ns(3,:),'k');
plot(clim_tn_gh2_ns(3,:),'color',[148/255 0 211/255]); %plot(clim_tn_w_1st_ns(3,:),'r'); 
plot(clim_tn_y_ns(3,:),'c');
% plot(clim_tn_j_1st_ns(3,:),'m');  
plot(clim_tn_ye_ns(3,:),'color',[255/255 192/255 203/255]);  plot(clim_tn_jw_ns(3,:),'color',[255/255 165/255 0]);
ylim([0 15]); yticks(1:15); yticklabels(1:15);
xticks(1:12);
xlim([1 12]); grid on;
xlabel('month')
ylabel('TN (mg/L)')
% legend(name_tag);
set(gca,'fontsize',12);



figure; hold on; 
plot(clim_tp_g_ns(3,:),'b'); plot(clim_tp_gc_ns(3,:),'g'); plot(clim_tp_gh_ns(3,:),'k');
plot(clim_tp_gh2_ns(3,:),'color',[148/255 0 211/255]); plot(clim_tp_w_1st_ns(3,:),'r'); plot(clim_tp_y_ns(3,:),'c');
 plot(clim_tp_j_1st_ns(3,:),'m');  plot(clim_tp_ye_ns(3,:),'color',[255/255 192/255 203/255]);  plot(clim_tp_jw_ns(3,:),'color',[255/255 165/255 0]);
% xticklabels(1997:2020);
xticks(1:12);
xlim([1 12]); grid on;
xlabel('month');
ylabel('TP (mg/L)')
legend(name_tag);
set(gca,'fontsize',12);
%  ylim([0 45]) % y axis match with raw

figure; hold on; 
plot(clim_tp_g_ns(3,:),'b'); plot(clim_tp_gc_ns(3,:),'g'); plot(clim_tp_gh_ns(3,:),'k');
plot(clim_tp_gh2_ns(3,:),'color',[148/255 0 211/255]); %plot(clim_tp_w_1st_ns(3,:),'r'); 
% plot(clim_tp_y_ns(3,:),'c');
% plot(clim_tp_j_1st_ns(3,:),'m'); 
plot(clim_tp_ye_ns(3,:),'color',[255/255 192/255 203/255]);  plot(clim_tp_jw_ns(3,:),'color',[255/255 165/255 0]);
% xticklabels(1997:2020);
xticks(1:12);
xlim([1 12]); grid on;
xlabel('month');
ylabel('TP (mg/L)')
% legend(name_tag);
ylim([0 1.1])
set(gca,'fontsize',12);

figure; hold on; 
plot(clim_dis_g,'b'); plot(clim_dis_gc,'g'); plot(clim_dis_gh,'k');
plot(clim_dis_gh2,'color',[148/255 0 211/255]); plot(clim_dis,'r'); plot(clim_dis_y,'c');
 plot(clim_dis_j,'m');  plot(clim_dis_ye,'color',[255/255 192/255 203/255]);  plot(clim_dis_jw,'color',[255/255 165/255 0]);
xticks(1:12);
xlim([1 12]); grid on;
ylim([0 1500]);
xlabel('month');
ylabel('m^3/day')
legend(name_tag);
set(gca,'fontsize',12);

%% raw (no sigma extraction climate)
figure; hold on; 
plot(clim_tn_g_raw(3,:),'b'); plot(clim_tn_gc_raw(3,:),'g'); plot(clim_tn_gh_raw(3,:),'k');
plot(clim_tn_gh2_raw(3,:),'color',[148/255 0 211/255]); plot(clim_tn_w_1st,'r'); 
plot(clim_tn_y_raw(3,:),'c');
plot(clim_tn_j_1st,'m');  
plot(clim_tn_ye_raw(3,:),'color',[255/255 192/255 203/255]);  plot(clim_tn_jw_raw(3,:),'color',[255/255 165/255 0]);
ylim([0 15]); yticks(1:15); yticklabels(1:15);
xticks(1:12);
xlim([1 12]); grid on;
xlabel('month')
ylabel('TN (mg/L)')
% legend(name_tag);
set(gca,'fontsize',12);


figure; hold on; 
plot(clim_tp_g_raw(3,:),'b'); plot(clim_tp_gc_raw(3,:),'g'); plot(clim_tp_gh_raw(3,:),'k');
plot(clim_tp_gh2_raw(3,:),'color',[148/255 0 211/255]); plot(clim_tp_w_1st,'r'); %plot(clim_tp_y_raw(3,:),'c');
 plot(clim_tp_j_1st,'m');  plot(clim_tp_ye_raw(3,:),'color',[255/255 192/255 203/255]);  plot(clim_tp_jw_raw(3,:),'color',[255/255 165/255 0]);
% xticklabels(1997:2020);
ylim([0 1.1]);
xticks(1:12);
xlim([1 12]); grid on;
xlabel('month');
ylabel('TP (mg/L)')
legend(name_tag);
set(gca,'fontsize',12);




% 여수-율촌산업단지 time (TP)
figure; hold on; 
plot(tp_y_2set); %plot(po4_sur./1000, 'g+');
plot(tp_y_ns(3,:),'g');% plot(tn_j_2set,'r');
xticklabels(1997:2020);
xticks(1:12:288);
xlim([169 288]); grid on;
xtickangle(45);
ylabel('TP (mg/L)')
legend('여수-율촌산업단지','율촌산업단지-3sig');
set(gca,'fontsize',12);
ylim([0 150])




