close all; clear; clc;

% 광양항, 광양4, 광양3, 광양2, 광양1, 광양5, 여수2, 여수3, 여수1

% 2011~ 이후 lon, lat 및 해당 지점 수심 기입되기 시작

%% 1997~2012
% 광양항 = 314
% 광양 1:5 = 136:140
% 여수 1:3 = 145:147
%% 2013~2019
%2013~
% 2월 (not 2월_ppb)

% 2013
% 광양항 = 354
% 광양 1:5 =159:163
% 여수 1:3 = 171:173

% 2014~2021
% 광양항 = 384
% 광양 1:5 = 174:178
% 여수 1:3 = 186:188

%% make 1997 to 2019 'yymm' form
k=0
for i = 1997:2021
    for j = 1:12
        k=k+1;
        ref_date{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end 
    mon_d=NaN(16,length(ref_date)); % day on month
    time=NaN(16,length(ref_date)); % time on day
    temp_sur=NaN(16,length(ref_date));
    temp_bot=NaN(16,length(ref_date));
    salt_sur=NaN(16,length(ref_date));
    salt_bot=NaN(16,length(ref_date));
    no3_sur=NaN(16,length(ref_date));
    no3_bot=NaN(16,length(ref_date));
    nh4_sur=NaN(16,length(ref_date));
    nh4_bot=NaN(16,length(ref_date));
    DIN_sur=NaN(16,length(ref_date));
    DIN_bot=NaN(16,length(ref_date));
    chl_sur=NaN(16,length(ref_date));
    chl_bot=NaN(16,length(ref_date));
    po4_sur=NaN(16,length(ref_date));
    po4_bot=NaN(16,length(ref_date));
    do_sur=NaN(16,length(ref_date));
    do_bot=NaN(16,length(ref_date)); 
    ss_sur=NaN(16,length(ref_date));
    ss_bot=NaN(16,length(ref_date)); 
    si_sur=NaN(16,length(ref_date));
    si_bot=NaN(16,length(ref_date)); 
    secchi_sur=NaN(16,length(ref_date));

set_path = 'D:\장기생태\Dynamic\KOEM\koem_yoonja_kang\측정망-메이스(문화)\';

for i = 1997:2021
        clearvars raw2 txt2 raw5 txt5 raw8 txt8 raw11 txt11 obs_tdx sheet_tail
    if i <= 2010
        col_temp = [7, 8];
        col_salt = [9, 10];
        col_no3 = [21, 22];
        col_nh4 = [17, 18];
        col_chl = [35, 36];
        col_po4 = [27, 28];
        col_do = [13, 14];
        col_DIN = [23, 24];
        col_ss=[33, 34];
        col_si=[31. 32];
        col_secchi=[37];
        row_pick = [314, 136:144, 368:370,  145:147]' - 2;   %368:370 = NaN 
        sheet_tail = '_ppb';
    elseif i > 2010 & i <= 2012
        col_time = [6, 7]; %day, time
        col_temp = [13, 14];
        col_salt = [15, 16];
        col_no3 = [27, 28];
        col_nh4 = [23, 24];
        col_chl = [41, 42];
        col_po4 = [33, 34];
        col_do = [19, 20];
        col_DIN = [29, 30];
        col_ss=[39, 40];
        col_si=[37. 38];
        col_secchi=[43];
        row_pick = [314, 136:144, 368:370, 145:147]'- 2; %368:370 = NaN 
        sheet_tail = '_ppb';
    elseif i == 2013
        col_time = [6, 7]; %day, time
        col_temp = [13, 14];
        col_salt = [15, 16];
        col_no3 = [27, 28];
        col_nh4 = [23, 24];
        col_chl = [41, 42];
        col_po4 = [33, 34];
        col_do = [19, 20];
        col_DIN = [29, 30];
        col_ss=[39, 40];
        col_si=[37. 38];
        col_secchi=[43];
        row_pick = [354, 159:170, 171:173]'- 2;
        sheet_tail = [];
    elseif i >= 2014
        col_time = [6, 7]; %day, time
        col_temp = [13, 14];
        col_salt = [15, 16];
        col_no3 = [27, 28];
        col_nh4 = [23, 24];
        col_chl = [41, 42];
        col_po4 = [33, 34];
        col_do = [19, 20];
        col_DIN = [29, 30];
        col_ss=[39, 40];
        col_si=[37. 38];
        col_secchi=[43];
        row_pick = [384, 174:185, 186:188]'- 2;
        sheet_tail = [];
    end
        

    [raw2 txt2]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['2월',sheet_tail]);
    [raw5 txt5]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['5월',sheet_tail]);
    [raw8 txt8]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['8월',sheet_tail]);
    [raw11 txt11]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['11월',sheet_tail]);
    
    obs_tdx = 2+((i-1997)*12):3:12+((i-1997)*12);
    
    % 2mon
    if i > 2010
    mon_d(:,obs_tdx(1))=raw2(row_pick,col_time(1)); % day on month
    time(:,obs_tdx(1))=raw2(row_pick,col_time(2))*24; % time on day
    end
    temp_sur(:,obs_tdx(1))=raw2(row_pick,col_temp(1));
    temp_bot(:,obs_tdx(1))=raw2(row_pick,col_temp(2));
    salt_sur(:,obs_tdx(1))=raw2(row_pick,col_salt(1));
    salt_bot(:,obs_tdx(1))=raw2(row_pick,col_salt(2));
    no3_sur(:,obs_tdx(1))=raw2(row_pick,col_no3(1));
    no3_bot(:,obs_tdx(1))=raw2(row_pick,col_no3(2));
    nh4_sur(:,obs_tdx(1))=raw2(row_pick,col_nh4(1));
    nh4_bot(:,obs_tdx(1))=raw2(row_pick,col_nh4(2));
    chl_sur(:,obs_tdx(1))=raw2(row_pick,col_chl(1));
    chl_bot(:,obs_tdx(1))=raw2(row_pick,col_chl(2));
    po4_sur(:,obs_tdx(1))=raw2(row_pick,col_po4(1));
    po4_bot(:,obs_tdx(1))=raw2(row_pick,col_po4(2));
    do_sur(:,obs_tdx(1))=raw2(row_pick,col_do(1));
    do_bot(:,obs_tdx(1))=raw2(row_pick,col_do(2));
    DIN_sur(:,obs_tdx(1))=raw2(row_pick,col_DIN(1));
    DIN_bot(:,obs_tdx(1))=raw2(row_pick,col_DIN(2));
    ss_sur(:,obs_tdx(1))=raw2(row_pick,col_ss(1));
    ss_bot(:,obs_tdx(1))=raw2(row_pick,col_ss(2));   
    si_sur(:,obs_tdx(1))=raw2(row_pick,col_si(1));
    si_bot(:,obs_tdx(1))=raw2(row_pick,col_si(2));  
    secchi_sur(:,obs_tdx(1))=raw2(row_pick,col_secchi(1));
    
    % 5mon
    if i > 2010
    mon_d(:,obs_tdx(2))=raw5(row_pick,col_time(1)); % day on month
    time(:,obs_tdx(2))=raw5(row_pick,col_time(2))*24; % time on day   
    end
    temp_sur(:,obs_tdx(2))=raw5(row_pick,col_temp(1));
    temp_bot(:,obs_tdx(2))=raw5(row_pick,col_temp(2));
    salt_sur(:,obs_tdx(2))=raw5(row_pick,col_salt(1));
    salt_bot(:,obs_tdx(2))=raw5(row_pick,col_salt(2));
    no3_sur(:,obs_tdx(2))=raw5(row_pick,col_no3(1));
    no3_bot(:,obs_tdx(2))=raw5(row_pick,col_no3(2));
    nh4_sur(:,obs_tdx(2))=raw5(row_pick,col_nh4(1));
    nh4_bot(:,obs_tdx(2))=raw5(row_pick,col_nh4(2));
    chl_sur(:,obs_tdx(2))=raw5(row_pick,col_chl(1));
    chl_bot(:,obs_tdx(2))=raw5(row_pick,col_chl(2));
    po4_sur(:,obs_tdx(2))=raw5(row_pick,col_po4(1));
    po4_bot(:,obs_tdx(2))=raw5(row_pick,col_po4(2));
    do_sur(:,obs_tdx(2))=raw5(row_pick,col_do(1));
    do_bot(:,obs_tdx(2))=raw5(row_pick,col_do(2)); 
    DIN_sur(:,obs_tdx(2))=raw5(row_pick,col_DIN(1));
    DIN_bot(:,obs_tdx(2))=raw5(row_pick,col_DIN(2));
    ss_sur(:,obs_tdx(2))=raw5(row_pick,col_ss(1));
    ss_bot(:,obs_tdx(2))=raw5(row_pick,col_ss(2));   
    si_sur(:,obs_tdx(2))=raw5(row_pick,col_si(1));
    si_bot(:,obs_tdx(2))=raw5(row_pick,col_si(2));  
    secchi_sur(:,obs_tdx(2))=raw5(row_pick,col_secchi(1));
    
    % 8mon
    if i > 2010
    mon_d(:,obs_tdx(3))=raw8(row_pick,col_time(1)); % day on month
    time(:,obs_tdx(3))=raw8(row_pick,col_time(2))*24; % time on day
    end
    temp_sur(:,obs_tdx(3))=raw8(row_pick,col_temp(1));
    temp_bot(:,obs_tdx(3))=raw8(row_pick,col_temp(2));
    salt_sur(:,obs_tdx(3))=raw8(row_pick,col_salt(1));
    salt_bot(:,obs_tdx(3))=raw8(row_pick,col_salt(2));
    no3_sur(:,obs_tdx(3))=raw8(row_pick,col_no3(1));
    no3_bot(:,obs_tdx(3))=raw8(row_pick,col_no3(2));
    nh4_sur(:,obs_tdx(3))=raw8(row_pick,col_nh4(1));
    nh4_bot(:,obs_tdx(3))=raw8(row_pick,col_nh4(2));
    chl_sur(:,obs_tdx(3))=raw8(row_pick,col_chl(1));
    chl_bot(:,obs_tdx(3))=raw8(row_pick,col_chl(2));
    po4_sur(:,obs_tdx(3))=raw8(row_pick,col_po4(1));
    po4_bot(:,obs_tdx(3))=raw8(row_pick,col_po4(2));
    do_sur(:,obs_tdx(3))=raw8(row_pick,col_do(1));
    do_bot(:,obs_tdx(3))=raw8(row_pick,col_do(2)); 
    DIN_sur(:,obs_tdx(3))=raw8(row_pick,col_DIN(1));
    DIN_bot(:,obs_tdx(3))=raw8(row_pick,col_DIN(2));
    ss_sur(:,obs_tdx(3))=raw8(row_pick,col_ss(1));
    ss_bot(:,obs_tdx(3))=raw8(row_pick,col_ss(2));   
    si_sur(:,obs_tdx(3))=raw8(row_pick,col_si(1));
    si_bot(:,obs_tdx(3))=raw8(row_pick,col_si(2));  
    secchi_sur(:,obs_tdx(3))=raw8(row_pick,col_secchi(1));
    
    % 11mon
    if i > 2010
    mon_d(:,obs_tdx(4))=raw11(row_pick,col_time(1)); % day on month
    time(:,obs_tdx(4))=raw11(row_pick,col_time(2))*24; % time on day
    end
    temp_sur(:,obs_tdx(4))=raw11(row_pick,col_temp(1));
    temp_bot(:,obs_tdx(4))=raw11(row_pick,col_temp(2));
    salt_sur(:,obs_tdx(4))=raw11(row_pick,col_salt(1));
    salt_bot(:,obs_tdx(4))=raw11(row_pick,col_salt(2));
    no3_sur(:,obs_tdx(4))=raw11(row_pick,col_no3(1));
    no3_bot(:,obs_tdx(4))=raw11(row_pick,col_no3(2));
    nh4_sur(:,obs_tdx(4))=raw11(row_pick,col_nh4(1));
    nh4_bot(:,obs_tdx(4))=raw11(row_pick,col_nh4(2));
%     if size(raw11,2) < 36
%         chl_sur(:,obs_tdx(4))=NaN(9,1);
%         chl_bot(:,obs_tdx(4))=NaN(9,1);
%     else
        chl_sur(:,obs_tdx(4))=raw11(row_pick,col_chl(1));
        chl_bot(:,obs_tdx(4))=raw11(row_pick,col_chl(2));
%     end
    po4_sur(:,obs_tdx(4))=raw11(row_pick,col_po4(1));
    po4_bot(:,obs_tdx(4))=raw11(row_pick,col_po4(2));
    do_sur(:,obs_tdx(4))=raw11(row_pick,col_do(1));
    do_bot(:,obs_tdx(4))=raw11(row_pick,col_do(2)); 
    DIN_sur(:,obs_tdx(4))=raw11(row_pick,col_DIN(1));
    DIN_bot(:,obs_tdx(4))=raw11(row_pick,col_DIN(2)); 
    ss_sur(:,obs_tdx(4))=raw11(row_pick,col_ss(1));
    ss_bot(:,obs_tdx(4))=raw11(row_pick,col_ss(2));   
    si_sur(:,obs_tdx(4))=raw11(row_pick,col_si(1));
    si_bot(:,obs_tdx(4))=raw11(row_pick,col_si(2));  
    secchi_sur(:,obs_tdx(4))=raw11(row_pick,col_secchi(1));
    
    
end

save('yoonjakangs_koem_data_monthly_v3_16points_2021.mat');

return
figure; hold on;
for i = 1:9
plot(2:3:length(ref_date), nh4_sur(:,2:3:end));
end

figure; 
% plot(2:3:length(ref_date), po4_sur(:,2:3:end) ./ 30.973762);
plot(2:3:length(ref_date), po4_sur(:,2:3:end) ./94.971482);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 5]);
% legend('광양1', '광양2', '광양3', '광양4', '광양5');
xtickangle(45);



%     [raw20 txt20]=xlsread([set_path,num2str(2010),'_ppb.xlsx'],'2월_ppb');
%     [raw21 txt21]=xlsread([set_path,num2str(2011),'_ppb.xlsx'],'2월_ppb');
%     [raw1 txt1]=xlsread([set_path,num2str(2019),'_ppb.xlsx'],'2월');