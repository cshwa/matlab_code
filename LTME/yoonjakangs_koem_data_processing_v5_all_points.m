close all; clear; clc;

set_path = 'D:\장기생태\Dynamic\KOEM\koem_yoonja_kang\측정망-메이스(문화)\';

sheet_tail = [];

[raw2 txt2]=xlsread([set_path,'1997_ppb.xlsx'],['2월_ppb']);

[raw0 txt0]=xlsread([set_path,'2021_ppb.xlsx'],['2월',sheet_tail]);

% st. name = col 4 + col 5
clearvars st_ref
for i = 1:length(raw0(2:end,5))
    if length(char(txt0(i+3,4))) == 0  % col 4 empty than fill as latest name
%         st_ref_pre{i,1} = [char(txt2(latest_row,4)), num2str(raw2(i+1,5))];
        if isnan(raw0(i+1,5)) == 0
            st_ref{i,1} = [char(txt0(latest_row,4)), num2str(raw0(i+1,5))];
        elseif isnan(raw0(i+1,5)) == 1
            st_ref{i,1} = [char(txt0(latest_row,4)), char(txt0(i+3,5))];
        end
    elseif length(char(txt0(i+3,4))) ~= 0
        latest_row = i+3;
        if isnan(raw0(i+1,5)) == 0
            st_ref{i,1} = [char(txt0(i+3,4)), num2str(raw0(i+1,5))];
        elseif isnan(raw0(i+1,5)) == 1
            st_ref{i,1} = [char(txt0(i+3,4)), char(txt0(i+3,5))];
        end
    end
end


lat=txt0(4:end,9);
lon=txt0(4:end,10);
k=0;
for i = 1:length(lat)
    if length(lat{i}) > 0        
         lat_tempo = char(lat{i});
         lon_tempo = char(lon{i});
         
         lat_num(i) = str2num(lat_tempo(1,1:2)) + str2num(lat_tempo(1,5:6))/60 + str2num(lat_tempo(1,9:10))/3600;
         lon_num(i) = str2num(lon_tempo(1,1:3)) + str2num(lon_tempo(1,6:7))/60 + str2num(lon_tempo(1,10:11))/3600;
    elseif length(lat{i}) == 0
         k = k+1;
         nandx(k) = i;
         lat_num(i) = NaN; 
         lon_num(i) = NaN;
    end
end

% temp_sur(nandx+1,:)

for i = 1:length(lat)
    if length(lat{i}) > 0        
         lat_m{i} = char(lat{i});
         lon_m{i} = char(lon{i});
    elseif length(lat{i}) == 0
         k = k+1;
         nandx(k) = i;
         lat_m{i} = NaN; 
         lon_m{i} = NaN;
    end
end



latlim = [33.0 39.5]; 
lonlim = [124 130];
% latlim = [34.0 35.3]; 
% lonlim = [126.35 129.00];

numberOfAttempts = 5;
attempt = 0;
info = [];
serverURL = 'http://basemap.nationalmap.gov/ArcGIS/services/USGSImageryOnly/MapServer/WMSServer?';
while(isempty(info))
    try
        info = wmsinfo(serverURL);
        orthoLayer = info.Layer(1);
    catch e 
        
        attempt = attempt + 1;
        if attempt > numberOfAttempts
            throw(e);
        else
            fprintf('Attempting to connect to server:\n"%s"\n', serverURL)
        end        
    end
end

imageLength = 1024;
[A,R] = wmsread(orthoLayer,'Latlim',latlim, ...
                           'Lonlim',lonlim, ...
                           'ImageHeight',imageLength, ...
                           'ImageWidth',imageLength);

% axesm('utm', ...
%       'Zone',utmzone(latlim, lonlim), ...
%       'MapLatlimit', latlim, ...
%       'MapLonlimit', lonlim, ...
%       'Geoid', wgs84Ellipsoid)

axesm('MapProjection','mercator',...
      'MapLatlimit', latlim, ...
      'MapLonlimit', lonlim, ...
      'Geoid', wgs84Ellipsoid)

% figure; hold on;
% geoshow(A,R)
% % framem
% % gridm
% ylim([33 39.5])
% title({'KOEM station'})
% set(gca, 'fontsize',13);
% saveas(gcf,strcat(save_fig_path,'monthly_mean_yjtak_',num2str(file_num),'m_bias_cen_obs.png'),'png');
% close; 

fontSizeTick = 13
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
geoshow(A,R); ylim([33 38.5])
xlim([125.5 130]);
title({'KOEM station'})
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% for i = 1:465
% text(lon_num(i),lat_num(i), num2str(i),'fontsize',8,'fontweight','bold','color','w') ;
% end
axis off
scatter(lon_num,lat_num,'r.')
saveas(gcf,strcat('KOEM_all_st.png'),'png');
close; 

lat = [34.2 35.1];
lon = [126.5 128.50];                       
fontSizeTick = 13
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;                       
for i = 1:465
text(lon_num(i),lat_num(i), num2str(i),'fontsize',8,'fontweight','bold','color','k') ;
end
% axis off
scatter(lon_num,lat_num,'r.')
ylim(lat)
xlim(lon)
plot_google_map('MapType','terrain','Scale',2,'Resize',2,'ShowLabels',0)   % overlay google map
ylim(lat)
xlim(lon) 
% [roadmap, satellite, terrain, hybrid]


% %% test %%
% clearvars st_pick len_raw_i
% len_raw_j=length(raw2(2:end,5)) - length(find(raw2(:,1) == 1)); %% cut extra row.
% for j = 1:len_raw_j
%     if length(txt2(:,5)) >= j+3 % i-index not exceed col 4 size
%         if length(char(txt2(j+3,4))) == 0 % col 4 empty than fill as latest name
%     %         st_pick_pre{i,1} = [char(txt2(latest_row,4)), num2str(raw2(i+1,5))];
%             if isnan(raw2(j+1,5)) == 0 
%                 st_pick{j,1} = [char(txt2(latest_row,4)), num2str(raw2(j+1,5))];
%             elseif isnan(raw2(j+1,5)) == 1 % st. num empty than H? on the txt.
%                 st_pick{j,1} = [char(txt2(latest_row,4)), char(txt2(j+3,5))];
%             end
%         elseif length(char(txt2(j+3,4))) ~= 0 % col 4 empty not empty
%             latest_row = j+3;
%             if isnan(raw2(j+1,5)) == 0 
%                 st_pick{j,1} = [char(txt2(j+3,4)), num2str(raw2(j+1,5))];
%             elseif isnan(raw2(j+1,5)) == 1 % st. num empty than H? on the txt.
%                 st_pick{j,1} = [char(txt2(j+3,4)), char(txt2(j+3,5))];
%             end
%         end
%     elseif length(txt2(:,5)) < j+3 % i-index exceed col 4 size
%             if isnan(raw2(j+1,5)) == 0 
%                 st_pick{j,1} = [char(txt2(latest_row,4)), num2str(raw2(j+1,5))];
%             elseif isnan(raw2(i+1,5)) == 1 % st. num empty than H? on the txt.
%                 st_pick{j,1} = [char(txt2(latest_row,4)), char(txt2(j+3,5))];
%             end
%     end       
% end
% 
% 
% for k = 1:length(st_ref)
% 
%     tempo_check=strcmp(st_ref{k},st_pick);
%     
%     check_idx = find(tempo_check == 1);
%     
%     if length(check_idx) == 1
%         row_pick(k) = check_idx;
%     elseif length(check_idx) == 0
%         row_pick(k) = NaN;
%     end
% end
% 
% temp_sur=NaN(length(st_ref),length(ref_date));
% i = 1997; col_temp = [7, 8];
% obs_tdx = 2+((i-1997)*12):3:12+((i-1997)*12);
% 
% for l = 1:length(row_pick)
%     if isnan(row_pick(l)) == 0
%         temp_sur(l,obs_tdx(1))=raw2(row_pick(l),col_temp(1));
%     else isnan(row_pick(l)) == 1
%         temp_sur(l,obs_tdx(1))=NaN;
%     end
% end

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
    mon_d=NaN(length(st_ref),length(ref_date)); % day on month
    time=NaN(length(st_ref),length(ref_date)); % time on day
    temp_sur=NaN(length(st_ref),length(ref_date));
    temp_bot=NaN(length(st_ref),length(ref_date));
    salt_sur=NaN(length(st_ref),length(ref_date));
    salt_bot=NaN(length(st_ref),length(ref_date));
    no3_sur=NaN(length(st_ref),length(ref_date));
    no3_bot=NaN(length(st_ref),length(ref_date));
    nh4_sur=NaN(length(st_ref),length(ref_date));
    nh4_bot=NaN(length(st_ref),length(ref_date));
    DIN_sur=NaN(length(st_ref),length(ref_date));
    DIN_bot=NaN(length(st_ref),length(ref_date));
    ph_sur=NaN(length(st_ref),length(ref_date));
    ph_bot=NaN(length(st_ref),length(ref_date));
    DIP_sur=NaN(length(st_ref),length(ref_date));
    DIP_bot=NaN(length(st_ref),length(ref_date));
    cod_sur=NaN(length(st_ref),length(ref_date));
    cod_bot=NaN(length(st_ref),length(ref_date));
    TN_sur=NaN(length(st_ref),length(ref_date));
    TN_bot=NaN(length(st_ref),length(ref_date));
    TP_sur=NaN(length(st_ref),length(ref_date));
    TP_bot=NaN(length(st_ref),length(ref_date)); 
    chl_sur=NaN(length(st_ref),length(ref_date));
    chl_bot=NaN(length(st_ref),length(ref_date));
    po4_sur=NaN(length(st_ref),length(ref_date));
    po4_bot=NaN(length(st_ref),length(ref_date));
    do_sur=NaN(length(st_ref),length(ref_date));
    do_bot=NaN(length(st_ref),length(ref_date)); 
    ss_sur=NaN(length(st_ref),length(ref_date));
    ss_bot=NaN(length(st_ref),length(ref_date)); 
    si_sur=NaN(length(st_ref),length(ref_date));
    si_bot=NaN(length(st_ref),length(ref_date)); 
    secchi_sur=NaN(length(st_ref),length(ref_date));


l2 = 0; l5 = 0; l8 = 0; l11 = 0;
for i = 1997:2021
        clearvars raw2 txt2 raw5 txt5 raw8 txt8 raw11 txt11 obs_tdx sheet_tail
    if i <= 2010
        col_temp = [7, 8];
        col_salt = [9, 10];
        col_ph = [11, 12];
        col_do = [13, 14];
        col_cod = [15, 16];
        col_no3 = [21, 22];
        col_nh4 = [17, 18];
        col_chl = [35, 36];
        col_po4 = [27, 28];
        col_do = [13, 14];
        col_no2 = [19, 20];
        col_DIN = [23, 24];
        col_TN = [25 26];
        col_dip = [27 28];
        col_TP = [29 30];        
        col_ss=[33, 34];
        col_si=[31, 32];
        col_secchi=[37];
        sheet_tail = '_ppb';
    elseif i > 2010 & i <= 2012
        col_time = [6, 7]; %day, time
        col_temp = [13, 14];
        col_salt = [15, 16];
        col_ph = [17, 18];
        col_no3 = [27, 28];
        col_nh4 = [23, 24];
        col_chl = [41, 42];
        col_po4 = [33, 34];
        col_do = [19, 20];
        col_cod = [21, 22];
        col_DIN = [29, 30];
        col_no2 = [25, 26];
        col_TN = [31,32];
        col_dip = [33,34];
        col_TP = [35,36];
        col_ss=[39, 40];
        col_si=[37, 38];
        col_secchi=[43];
        sheet_tail = '_ppb';
    elseif i >= 2013
        col_time = [6, 7]; %day, time
        col_temp = [13, 14];
        col_salt = [15, 16];
        col_ph = [17, 18];
        col_no3 = [27, 28];
        col_nh4 = [23, 24];
        col_chl = [41, 42];
        col_po4 = [33, 34];
        col_do = [19, 20];
        col_cod = [21, 22];
        col_DIN = [29, 30];
        col_no2 = [25, 26];
        col_TN = [31,32];
        col_dip = [33,34];
        col_TP = [35,36];
        col_ss=[39, 40];
        col_si=[37, 38];
        col_secchi=[43]; 
        sheet_tail = [];
    end
        
    [raw2 txt2]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['2월',sheet_tail]);
    [raw5 txt5]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['5월',sheet_tail]);
    [raw8 txt8]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['8월',sheet_tail]);
    [raw11 txt11]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['11월',sheet_tail]);
    
    
    %% load each year-month name list
    clearvars st_pick len_raw_j
    len_raw_j=length(raw2(2:end,5)) - length(find(raw2(:,1) == 1)); %% cut extra row.
    for j = 1:len_raw_j
        if length(txt2(:,5)) >= j+3 % i-index not exceed col 4 size
            if length(char(txt2(j+3,4))) == 0 % col 4 empty than fill as latest name
        %         st_pick_pre{i,1} = [char(txt2(latest_row,4)), num2str(raw2(i+1,5))];
                if isnan(raw2(j+1,5)) == 0 
                    st_pick{j,1} = [char(txt2(latest_row,4)), num2str(raw2(j+1,5))];
                elseif isnan(raw2(j+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick{j,1} = [char(txt2(latest_row,4)), char(txt2(j+3,5))];
                end
            elseif length(char(txt2(j+3,4))) ~= 0 % col 4 empty not empty
                latest_row = j+3;
                if isnan(raw2(j+1,5)) == 0 
                    st_pick{j,1} = [char(txt2(j+3,4)), num2str(raw2(j+1,5))];
                elseif isnan(raw2(j+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick{j,1} = [char(txt2(j+3,4)), char(txt2(j+3,5))];
                end
            end
        elseif length(txt2(:,5)) < j+3 % i-index exceed col 4 size
                if isnan(raw2(j+1,5)) == 0 
                    st_pick{j,1} = [char(txt2(latest_row,4)), num2str(raw2(j+1,5))];
                elseif isnan(raw2(i+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick{j,1} = [char(txt2(latest_row,4)), char(txt2(j+3,5))];
                end
        end       
    end
    
    
        clearvars st_pick_5m len_raw_j
    len_raw_j=length(raw5(2:end,5)) - length(find(raw5(:,1) == 1)); %% cut extra row.
    for j = 1:len_raw_j
        if length(txt5(:,5)) >= j+3 % i-index not exceed col 4 size
            if length(char(txt5(j+3,4))) == 0 % col 4 empty than fill as latest name
        %         st_pick_pre{i,1} = [char(txt5(latest_row,4)), num2str(raw5(i+1,5))];
                if isnan(raw5(j+1,5)) == 0 
                    st_pick_5m{j,1} = [char(txt5(latest_row,4)), num2str(raw5(j+1,5))];
                elseif isnan(raw5(j+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick_5m{j,1} = [char(txt5(latest_row,4)), char(txt5(j+3,5))];
                end
            elseif length(char(txt5(j+3,4))) ~= 0 % col 4 empty not empty
                latest_row = j+3;
                if isnan(raw5(j+1,5)) == 0 
                    st_pick_5m{j,1} = [char(txt5(j+3,4)), num2str(raw5(j+1,5))];
                elseif isnan(raw5(j+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick_5m{j,1} = [char(txt5(j+3,4)), char(txt5(j+3,5))];
                end
            end
        elseif length(txt5(:,5)) < j+3 % i-index exceed col 4 size
                if isnan(raw5(j+1,5)) == 0 
                    st_pick_5m{j,1} = [char(txt5(latest_row,4)), num2str(raw5(j+1,5))];
                elseif isnan(raw5(i+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick_5m{j,1} = [char(txt5(latest_row,4)), char(txt5(j+3,5))];
                end
        end       
    end
    
    
        clearvars st_pick_8m len_raw_j
    len_raw_j=length(raw8(2:end,5)) - length(find(raw8(:,1) == 1)); %% cut extra row.
    for j = 1:len_raw_j
        if length(txt8(:,5)) >= j+3 % i-index not exceed col 4 size
            if length(char(txt8(j+3,4))) == 0 % col 4 empty than fill as latest name
        %         st_pick_8m_pre{i,1} = [char(txt8(latest_row,4)), num2str(raw8(i+1,5))];
                if isnan(raw8(j+1,5)) == 0 
                    st_pick_8m{j,1} = [char(txt8(latest_row,4)), num2str(raw8(j+1,5))];
                elseif isnan(raw8(j+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick_8m{j,1} = [char(txt8(latest_row,4)), char(txt8(j+3,5))];
                end
            elseif length(char(txt8(j+3,4))) ~= 0 % col 4 empty not empty
                latest_row = j+3;
                if isnan(raw8(j+1,5)) == 0 
                    st_pick_8m{j,1} = [char(txt8(j+3,4)), num2str(raw8(j+1,5))];
                elseif isnan(raw8(j+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick_8m{j,1} = [char(txt8(j+3,4)), char(txt8(j+3,5))];
                end
            end
        elseif length(txt8(:,5)) < j+3 % i-index exceed col 4 size
                if isnan(raw8(j+1,5)) == 0 
                    st_pick_8m{j,1} = [char(txt8(latest_row,4)), num2str(raw8(j+1,5))];
                elseif isnan(raw8(i+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick_8m{j,1} = [char(txt8(latest_row,4)), char(txt8(j+3,5))];
                end
        end       
    end
    
    
        clearvars st_pick_11m len_raw_j
    len_raw_j=length(raw11(2:end,5)) - length(find(raw11(:,1) == 1)); %% cut extra row.
    for j = 1:len_raw_j
        if length(txt11(:,5)) >= j+3 % i-index not exceed col 4 size
            if length(char(txt11(j+3,4))) == 0 % col 4 empty than fill as latest name
        %         st_pick_11m_pre{i,1} = [char(txt11(latest_row,4)), num2str(raw11(i+1,5))];
                if isnan(raw11(j+1,5)) == 0 
                    st_pick_11m{j,1} = [char(txt11(latest_row,4)), num2str(raw11(j+1,5))];
                elseif isnan(raw11(j+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick_11m{j,1} = [char(txt11(latest_row,4)), char(txt11(j+3,5))];
                end
            elseif length(char(txt11(j+3,4))) ~= 0 % col 4 empty not empty
                latest_row = j+3;
                if isnan(raw11(j+1,5)) == 0 
                    st_pick_11m{j,1} = [char(txt11(j+3,4)), num2str(raw11(j+1,5))];
                elseif isnan(raw11(j+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick_11m{j,1} = [char(txt11(j+3,4)), char(txt11(j+3,5))];
                end
            end
        elseif length(txt11(:,5)) < j+3 % i-index exceed col 4 size
                if isnan(raw11(j+1,5)) == 0 
                    st_pick_11m{j,1} = [char(txt11(latest_row,4)), num2str(raw11(j+1,5))];
                elseif isnan(raw11(i+1,5)) == 1 % st. num empty than H? on the txt.
                    st_pick_11m{j,1} = [char(txt11(latest_row,4)), char(txt11(j+3,5))];
                end
        end       
    end
    
    %% make row_pick
    clearvars tempo_check check_idx row_pick*
    
    for k = 1:length(st_ref)

    tempo_check=strcmp(st_ref{k},st_pick);
    
    check_idx = find(tempo_check == 1);
    
        if length(check_idx) >= 1               
            if length(check_idx) > 1 %if there are multiple matches
                disp(st_ref{k});    
                if check_idx(end) < k  %if matrix size smaller than ref
                row_pick(k) = check_idx(end);    
                elseif check_idx(end) >= k 
                row_pick(k) = check_idx(find(k == check_idx));
                end
           else
                row_pick(k) = check_idx;
            end
        elseif length(check_idx) == 0
            row_pick(k) = NaN;
        end
    end
    
    clearvars tempo_check check_idx 
        for k = 1:length(st_ref)

    tempo_check=strcmp(st_ref{k},st_pick_5m);
    
    check_idx = find(tempo_check == 1);
    
        if length(check_idx) >= 1
            if length(check_idx) > 1
                disp(st_ref{k});    
                if check_idx(end) < k 
                row_pick_5m(k) = check_idx(end);    
                elseif check_idx(end) >= k 
                row_pick_5m(k) = check_idx(find(k == check_idx));
                end
           else
                row_pick_5m(k) = check_idx;
            end
        elseif length(check_idx) == 0
            row_pick_5m(k) = NaN;
        end
        end
    
        clearvars tempo_check check_idx
        for k = 1:length(st_ref)

    tempo_check=strcmp(st_ref{k},st_pick_8m);
    
    check_idx = find(tempo_check == 1);
    
        if length(check_idx) >= 1
            if length(check_idx) > 1
                disp(st_ref{k});    
                if check_idx(end) < k 
                row_pick_8m(k) = check_idx(end);    
                elseif check_idx(end) >= k 
                row_pick_8m(k) = check_idx(find(k == check_idx));
                end
            else
                row_pick_8m(k) = check_idx;
            end
        elseif length(check_idx) == 0
            row_pick_8m(k) = NaN;
        end
        end
    
            clearvars tempo_check check_idx
        for k = 1:length(st_ref)

    tempo_check=strcmp(st_ref{k},st_pick_11m);
    
    check_idx = find(tempo_check == 1);
    
        if length(check_idx) >= 1
            if length(check_idx) > 1
                disp(st_ref{k});    
                if check_idx(end) < k 
                row_pick_11m(k) = check_idx(end);    
                elseif check_idx(end) >= k 
                row_pick_11m(k) = check_idx(find(k == check_idx));
                end
            else
                row_pick_11m(k) = check_idx;
            end
        elseif length(check_idx) == 0
            row_pick_11m(k) = NaN;
        end
        end
    
        
    obs_tdx = 2+((i-1997)*12):3:12+((i-1997)*12);
    
    % 2mon
    for l = 1:length(row_pick)
    
    if isnan(row_pick(l)) == 0
        if i > 2010
        mon_d(l,obs_tdx(1))=raw2(row_pick(l),col_time(1)); % day on month
        time(l,obs_tdx(1))=raw2(row_pick(l),col_time(2))*24; % time on day
        end
%         temp_sur(i,obs_tdx(1))=raw2(row_pick(l),col_temp(1));
        temp_sur(l,obs_tdx(1))=raw2(row_pick(l),col_temp(1));
        temp_bot(l,obs_tdx(1))=raw2(row_pick(l),col_temp(2));
        salt_sur(l,obs_tdx(1))=raw2(row_pick(l),col_salt(1));
        salt_bot(l,obs_tdx(1))=raw2(row_pick(l),col_salt(2));
        no3_sur(l,obs_tdx(1))=raw2(row_pick(l),col_no3(1));
        no3_bot(l,obs_tdx(1))=raw2(row_pick(l),col_no3(2));
        nh4_sur(l,obs_tdx(1))=raw2(row_pick(l),col_nh4(1));
        nh4_bot(l,obs_tdx(1))=raw2(row_pick(l),col_nh4(2));
        chl_sur(l,obs_tdx(1))=raw2(row_pick(l),col_chl(1));
        chl_bot(l,obs_tdx(1))=raw2(row_pick(l),col_chl(2));
        po4_sur(l,obs_tdx(1))=raw2(row_pick(l),col_po4(1));
        po4_bot(l,obs_tdx(1))=raw2(row_pick(l),col_po4(2));
        do_sur(l,obs_tdx(1))=raw2(row_pick(l),col_do(1));
        do_bot(l,obs_tdx(1))=raw2(row_pick(l),col_do(2));
        DIN_sur(l,obs_tdx(1))=raw2(row_pick(l),col_DIN(1));
        DIN_bot(l,obs_tdx(1))=raw2(row_pick(l),col_DIN(2));
        ss_sur(l,obs_tdx(1))=raw2(row_pick(l),col_ss(1));
        ss_bot(l,obs_tdx(1))=raw2(row_pick(l),col_ss(2));   
        si_sur(l,obs_tdx(1))=raw2(row_pick(l),col_si(1));
        si_bot(l,obs_tdx(1))=raw2(row_pick(l),col_si(2));      
        ph_sur(l,obs_tdx(1))=raw2(row_pick(l),col_ph(1));
        ph_bot(l,obs_tdx(1))=raw2(row_pick(l),col_ph(2));        
        cod_sur(l,obs_tdx(1))=raw2(row_pick(l),col_cod(1));
        cod_bot(l,obs_tdx(1))=raw2(row_pick(l),col_cod(2));        
        no2_sur(l,obs_tdx(1))=raw2(row_pick(l),col_no2(1));
        no2_bot(l,obs_tdx(1))=raw2(row_pick(l),col_no2(2));
        TN_sur(l,obs_tdx(1))=raw2(row_pick(l),col_TN(1));
        TN_bot(l,obs_tdx(1))=raw2(row_pick(l),col_TN(2));       
        TP_sur(l,obs_tdx(1))=raw2(row_pick(l),col_TP(1));
        TP_bot(l,obs_tdx(1))=raw2(row_pick(l),col_TP(2));      
        DIP_sur(l,obs_tdx(1))=raw2(row_pick(l),col_dip(1));
        DIP_bot(l,obs_tdx(1))=raw2(row_pick(l),col_dip(2));     
        secchi_sur(l,obs_tdx(1))=raw2(row_pick(l),col_secchi(1));
    else isnan(row_pick(l)) == 1
        temp_sur(l,obs_tdx(1))=NaN;
        temp_bot(l,obs_tdx(1))=NaN;
        salt_sur(l,obs_tdx(1))=NaN;
        salt_bot(l,obs_tdx(1))=NaN;
        no3_sur(l,obs_tdx(1))=NaN;
        no3_bot(l,obs_tdx(1))=NaN;
        nh4_sur(l,obs_tdx(1))=NaN;
        nh4_bot(l,obs_tdx(1))=NaN;
        chl_sur(l,obs_tdx(1))=NaN;
        chl_bot(l,obs_tdx(1))=NaN;
        po4_sur(l,obs_tdx(1))=NaN;
        po4_bot(l,obs_tdx(1))=NaN;
        do_sur(l,obs_tdx(1))=NaN;
        do_bot(l,obs_tdx(1))=NaN;
        DIN_sur(l,obs_tdx(1))=NaN;
        DIN_bot(l,obs_tdx(1))=NaN;
        ss_sur(l,obs_tdx(1))=NaN;
        ss_bot(l,obs_tdx(1))=NaN;   
        si_sur(l,obs_tdx(1))=NaN;
        si_bot(l,obs_tdx(1))=NaN;        
        ph_sur(l,obs_tdx(1))=NaN;
        ph_bot(l,obs_tdx(1))=NaN;        
        cod_sur(l,obs_tdx(1))=NaN;
        cod_bot(l,obs_tdx(1))=NaN;        
        no2_sur(l,obs_tdx(1))=NaN;
        no2_bot(l,obs_tdx(1))=NaN;
        TN_sur(l,obs_tdx(1))=NaN;
        TN_bot(l,obs_tdx(1))=NaN;       
        TP_sur(l,obs_tdx(1))=NaN;
        TP_bot(l,obs_tdx(1))=NaN;      
        DIP_sur(l,obs_tdx(1))=NaN;
        DIP_bot(l,obs_tdx(1))=NaN;     
        secchi_sur(l,obs_tdx(1))=NaN;
    end
    end

    
    % 5mon
       for l = 1:length(row_pick_5m)
      
    if isnan(row_pick_5m(l)) == 0
        if i > 2010
        mon_d(l,obs_tdx(2))=raw5(row_pick_5m(l),col_time(1)); % day on month
        time(l,obs_tdx(2))=raw5(row_pick_5m(l),col_time(2))*24; % time on day
        end  
%         temp_sur(i,obs_tdx(2))=raw5(row_pick_5m(l),col_temp(1));
        temp_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_temp(1));
        temp_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_temp(2));
        salt_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_salt(1));
        salt_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_salt(2));
        no3_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_no3(1));
        no3_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_no3(2));
        nh4_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_nh4(1));
        nh4_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_nh4(2));
        chl_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_chl(1));
        chl_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_chl(2));
        po4_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_po4(1));
        po4_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_po4(2));
        do_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_do(1));
        do_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_do(2));
        DIN_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_DIN(1));
        DIN_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_DIN(2));
        ss_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_ss(1));
        ss_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_ss(2));   
        si_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_si(1));
        si_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_si(2));      
        ph_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_ph(1));
        ph_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_ph(2));        
        cod_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_cod(1));
        cod_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_cod(2));        
        no2_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_no2(1));
        no2_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_no2(2));
        TN_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_TN(1));
        TN_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_TN(2));       
        TP_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_TP(1));
        TP_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_TP(2));      
        DIP_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_dip(1));
        DIP_bot(l,obs_tdx(2))=raw5(row_pick_5m(l),col_dip(2));     
        secchi_sur(l,obs_tdx(2))=raw5(row_pick_5m(l),col_secchi(1));
    else isnan(row_pick_5m(l)) == 1
        temp_sur(l,obs_tdx(2))=NaN;
        temp_bot(l,obs_tdx(2))=NaN;
        salt_sur(l,obs_tdx(2))=NaN;
        salt_bot(l,obs_tdx(2))=NaN;
        no3_sur(l,obs_tdx(2))=NaN;
        no3_bot(l,obs_tdx(2))=NaN;
        nh4_sur(l,obs_tdx(2))=NaN;
        nh4_bot(l,obs_tdx(2))=NaN;
        chl_sur(l,obs_tdx(2))=NaN;
        chl_bot(l,obs_tdx(2))=NaN;
        po4_sur(l,obs_tdx(2))=NaN;
        po4_bot(l,obs_tdx(2))=NaN;
        do_sur(l,obs_tdx(2))=NaN;
        do_bot(l,obs_tdx(2))=NaN;
        DIN_sur(l,obs_tdx(2))=NaN;
        DIN_bot(l,obs_tdx(2))=NaN;
        ss_sur(l,obs_tdx(2))=NaN;
        ss_bot(l,obs_tdx(2))=NaN;   
        si_sur(l,obs_tdx(2))=NaN;
        si_bot(l,obs_tdx(2))=NaN;        
        ph_sur(l,obs_tdx(2))=NaN;
        ph_bot(l,obs_tdx(2))=NaN;        
        cod_sur(l,obs_tdx(2))=NaN;
        cod_bot(l,obs_tdx(2))=NaN;        
        no2_sur(l,obs_tdx(2))=NaN;
        no2_bot(l,obs_tdx(2))=NaN;
        TN_sur(l,obs_tdx(2))=NaN;
        TN_bot(l,obs_tdx(2))=NaN;       
        TP_sur(l,obs_tdx(2))=NaN;
        TP_bot(l,obs_tdx(2))=NaN;      
        DIP_sur(l,obs_tdx(2))=NaN;
        DIP_bot(l,obs_tdx(2))=NaN;     
        secchi_sur(l,obs_tdx(2))=NaN;
    end
    end
    
    % 8mon

       for l = 1:length(row_pick_8m)
     
    if isnan(row_pick_8m(l)) == 0
         if i > 2010
        mon_d(l,obs_tdx(3))=raw8(row_pick_8m(l),col_time(1)); % day on month
        time(l,obs_tdx(3))=raw8(row_pick_8m(l),col_time(2))*24; % time on day
        end  
%         temp_sur(i,obs_tdx(3))=raw8(row_pick_8m(l),col_temp(1));
        temp_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_temp(1));
        temp_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_temp(2));
        salt_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_salt(1));
        salt_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_salt(2));
        no3_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_no3(1));
        no3_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_no3(2));
        nh4_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_nh4(1));
        nh4_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_nh4(2));
        chl_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_chl(1));
        chl_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_chl(2));
        po4_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_po4(1));
        po4_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_po4(2));
        do_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_do(1));
        do_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_do(2));
        DIN_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_DIN(1));
        DIN_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_DIN(2));
        ss_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_ss(1));
        ss_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_ss(2));   
        si_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_si(1));
        si_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_si(2));      
        ph_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_ph(1));
        ph_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_ph(2));        
        cod_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_cod(1));
        cod_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_cod(2));        
        no2_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_no2(1));
        no2_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_no2(2));
        TN_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_TN(1));
        TN_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_TN(2));       
        TP_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_TP(1));
        TP_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_TP(2));      
        DIP_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_dip(1));
        DIP_bot(l,obs_tdx(3))=raw8(row_pick_8m(l),col_dip(2));     
        secchi_sur(l,obs_tdx(3))=raw8(row_pick_8m(l),col_secchi(1));
    else isnan(row_pick_8m(l)) == 1
        temp_sur(l,obs_tdx(3))=NaN;
        temp_bot(l,obs_tdx(3))=NaN;
        salt_sur(l,obs_tdx(3))=NaN;
        salt_bot(l,obs_tdx(3))=NaN;
        no3_sur(l,obs_tdx(3))=NaN;
        no3_bot(l,obs_tdx(3))=NaN;
        nh4_sur(l,obs_tdx(3))=NaN;
        nh4_bot(l,obs_tdx(3))=NaN;
        chl_sur(l,obs_tdx(3))=NaN;
        chl_bot(l,obs_tdx(3))=NaN;
        po4_sur(l,obs_tdx(3))=NaN;
        po4_bot(l,obs_tdx(3))=NaN;
        do_sur(l,obs_tdx(3))=NaN;
        do_bot(l,obs_tdx(3))=NaN;
        DIN_sur(l,obs_tdx(3))=NaN;
        DIN_bot(l,obs_tdx(3))=NaN;
        ss_sur(l,obs_tdx(3))=NaN;
        ss_bot(l,obs_tdx(3))=NaN;   
        si_sur(l,obs_tdx(3))=NaN;
        si_bot(l,obs_tdx(3))=NaN;        
        ph_sur(l,obs_tdx(3))=NaN;
        ph_bot(l,obs_tdx(3))=NaN;        
        cod_sur(l,obs_tdx(3))=NaN;
        cod_bot(l,obs_tdx(3))=NaN;        
        no2_sur(l,obs_tdx(3))=NaN;
        no2_bot(l,obs_tdx(3))=NaN;
        TN_sur(l,obs_tdx(3))=NaN;
        TN_bot(l,obs_tdx(3))=NaN;       
        TP_sur(l,obs_tdx(3))=NaN;
        TP_bot(l,obs_tdx(3))=NaN;      
        DIP_sur(l,obs_tdx(3))=NaN;
        DIP_bot(l,obs_tdx(3))=NaN;     
        secchi_sur(l,obs_tdx(3))=NaN;
    end
    end
    
    % 11mon
       for l = 1:length(row_pick_11m)
      
    if isnan(row_pick_11m(l)) == 0
        if i > 2010
        mon_d(l,obs_tdx(4))=raw11(row_pick_11m(l),col_time(1)); % day on month
        time(l,obs_tdx(4))=raw11(row_pick_11m(l),col_time(2))*24; % time on day
        end  
%         temp_sur(i,obs_tdx(4))=raw11(row_pick_11m(l),col_temp(1));
        temp_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_temp(1));
        temp_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_temp(2));
        salt_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_salt(1));
        salt_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_salt(2));
        no3_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_no3(1));
        no3_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_no3(2));
        nh4_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_nh4(1));
        nh4_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_nh4(2));
        chl_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_chl(1));
        chl_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_chl(2));
        po4_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_po4(1));
        po4_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_po4(2));
        do_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_do(1));
        do_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_do(2));
        DIN_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_DIN(1));
        DIN_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_DIN(2));
        ss_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_ss(1));
        ss_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_ss(2));   
        si_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_si(1));
        si_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_si(2));      
        ph_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_ph(1));
        ph_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_ph(2));        
        cod_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_cod(1));
        cod_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_cod(2));        
        no2_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_no2(1));
        no2_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_no2(2));
        TN_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_TN(1));
        TN_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_TN(2));       
        TP_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_TP(1));
        TP_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_TP(2));      
        DIP_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_dip(1));
        DIP_bot(l,obs_tdx(4))=raw11(row_pick_11m(l),col_dip(2));     
        secchi_sur(l,obs_tdx(4))=raw11(row_pick_11m(l),col_secchi(1));
    else isnan(row_pick_11m(l)) == 1
        temp_sur(l,obs_tdx(4))=NaN;
        temp_bot(l,obs_tdx(4))=NaN;
        salt_sur(l,obs_tdx(4))=NaN;
        salt_bot(l,obs_tdx(4))=NaN;
        no3_sur(l,obs_tdx(4))=NaN;
        no3_bot(l,obs_tdx(4))=NaN;
        nh4_sur(l,obs_tdx(4))=NaN;
        nh4_bot(l,obs_tdx(4))=NaN;
        chl_sur(l,obs_tdx(4))=NaN;
        chl_bot(l,obs_tdx(4))=NaN;
        po4_sur(l,obs_tdx(4))=NaN;
        po4_bot(l,obs_tdx(4))=NaN;
        do_sur(l,obs_tdx(4))=NaN;
        do_bot(l,obs_tdx(4))=NaN;
        DIN_sur(l,obs_tdx(4))=NaN;
        DIN_bot(l,obs_tdx(4))=NaN;
        ss_sur(l,obs_tdx(4))=NaN;
        ss_bot(l,obs_tdx(4))=NaN;   
        si_sur(l,obs_tdx(4))=NaN;
        si_bot(l,obs_tdx(4))=NaN;        
        ph_sur(l,obs_tdx(4))=NaN;
        ph_bot(l,obs_tdx(4))=NaN;        
        cod_sur(l,obs_tdx(4))=NaN;
        cod_bot(l,obs_tdx(4))=NaN;        
        no2_sur(l,obs_tdx(4))=NaN;
        no2_bot(l,obs_tdx(4))=NaN;
        TN_sur(l,obs_tdx(4))=NaN;
        TN_bot(l,obs_tdx(4))=NaN;       
        TP_sur(l,obs_tdx(4))=NaN;
        TP_bot(l,obs_tdx(4))=NaN;      
        DIP_sur(l,obs_tdx(4))=NaN;
        DIP_bot(l,obs_tdx(4))=NaN;     
        secchi_sur(l,obs_tdx(4))=NaN;
    end
    end
    
    
end

save('yoonjakangs_koem_data_monthly_v5_all_points_2021.mat');


sig_extract =1;

if sig_extract == 1
%% extract over 3sig
sp_gy=size(do_bot,1); %% num of spatial point

for varlist = {'TN','TP','DIP','do','chl','ss','si','po4','no3','no2','nh4','DIN','temp','salt','secchi'}
        clearvars varname 
        varname = char(varlist);
        for j=1:sp_gy %% num of spatial point
            for sig=1:6 %% sigma
                clearvars data data_ref
                eval(['data_ref = ', varname,'_sur(j,:);']);
                eval(['data = ', varname,'_sur(j,:);']);
                data(data_ref > nanmean(data_ref) + sig*nanstd(data_ref)) =NaN;
                data(data_ref < nanmean(data_ref) - sig*nanstd(data_ref)) =NaN;
                eval([varname,'_s(j,:,sig) = data;']);  % extracted data OUT [sp_gy, T, sig]
                disp([varname,'_s'])
           
                if strcmp(varname,'secchi') == 0 % there are no secchi_bot
                    clearvars data data_ref
                    eval(['data_ref = ', varname,'_bot(j,:);']);
                    eval(['data = ', varname,'_bot(j,:);']);
                    data(data_ref > nanmean(data_ref) + sig*nanstd(data_ref)) =NaN;
                    data(data_ref < nanmean(data_ref) - sig*nanstd(data_ref)) =NaN;
                    eval([varname,'_b_s(j,:,sig) = data;']);  % extracted data OUT [sp_gy, T, sig]
                   disp([varname,'_b_s'])
                end
            end
        end
    end
end

for i=1:300/12
temp_s_yr(:,i,:)=mean(temp_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
temp_b_s_yr(:,i,:)=mean(temp_b_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
salt_s_yr(:,i,:)=mean(salt_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
salt_b_s_yr(:,i,:)=mean(salt_b_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');  
chl_s_yr(:,i,:)=mean(chl_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
chl_b_s_yr(:,i,:)=mean(chl_b_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
DIN_s_yr(:,i,:)=mean(DIN_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
TN_s_yr(:,i,:)=mean(TN_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
TP_s_yr(:,i,:)=mean(TP_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
no3_s_yr(:,i,:)=mean(no3_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
nh4_s_yr(:,i,:)=mean(nh4_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
DIP_s_yr(:,i,:)=mean(DIP_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
ss_s_yr(:,i,:)=mean(ss_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
si_s_yr(:,i,:)=mean(si_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
secchi_s_yr(:,i,:)=mean(secchi_s(:,(i-1)*12+1:i*12,:),[2],'omitnan');
end

save('plot_KOEM_1997to2021.mat');

return
close all; clear all; clc;
load plot_KOEM_1997to2021
sj=load('D:\장기생태\Dynamic\06_river\data\sj_1980to1996\plot_load_songjung_1989to2021.mat'); %make_3clim_river_06to15_v4_mixed_daily.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% river vs. KOEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TP load vs. chl
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_tp_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 20]); xlim([1 length(sj.yr_tp_load)]);
ylabel('TP load(g/s)');
plot(sj.yr_tp_load,'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('TP_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_tp_load)]);
saveas(gcf,strcat('TP_load_vs_bay_chl_fr1997.png'),'png');

%% monthly
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_tp_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 20]); xlim([1 length(sj.yr_tp_load)]);
ylabel('TP load(g/s)');
plot(sj.mon_tp_load(i:12:end),'k');
yyaxis right
plot(9:33,mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylim([0 15]);
ylabel('Chl.a(ug/L)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('TP_load_vs_bay_chl_',num2str(i),'mon.png'),'png');
xlim([9 inf]);
saveas(gcf,strcat('TP_load_vs_bay_chl_fr1997_',num2str(i),'mon.png'),'png');
end

%% DIP load vs. chl
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_po4_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 20]); xlim([1 length(sj.yr_tp_load)]);
ylabel('PO4 load(g/s)');
plot(sj.yr_po4_load,'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('PO4_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_tp_load)]);
saveas(gcf,strcat('PO4_load_vs_bay_chl_fr1997.png'),'png');

%% monthly
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_po4_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 20]); xlim([1 length(sj.yr_po4_load)]);
ylabel('PO4 load(g/s)');
plot(sj.mon_po4_load(i:12:end),'k');
yyaxis right
plot(9:33,mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylim([0 15]);
ylabel('Chl.a(ug/L)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('PO4_load_vs_bay_chl_',num2str(i),'mon.png'),'png');
xlim([9 inf]);
saveas(gcf,strcat('PO4_load_vs_bay_chl_fr1997_',num2str(i),'mon.png'),'png');
end

%% TN load vs. chl
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_tn_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 350]); 
xlim([1 length(sj.yr_tn_load)]);
ylabel('TN load(g/s)');
plot(sj.yr_tn_load,'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('TN_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_tp_load)]);
saveas(gcf,strcat('TN_load_vs_bay_chl_fr1997.png'),'png');

%% monthly
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_tn_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 1300]); xlim([1 length(sj.yr_tn_load)]);
ylabel('TN load(g/s)');
plot(sj.mon_tn_load(i:12:end),'k');
yyaxis right
plot(9:33,mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylim([0 15]);
ylabel('Chl.a(ug/L)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('TN_load_vs_bay_chl_',num2str(i),'mon.png'),'png');
xlim([9 inf]);
saveas(gcf,strcat('TN_load_vs_bay_chl_fr1997_',num2str(i),'mon.png'),'png');
end


%% DIN load vs. chl
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_no3_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 250]); 
xlim([1 length(sj.yr_no3_load)]);
ylabel('DIN load(g/s)');
plot(sj.yr_no3_load+sj.yr_nh4_load,'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('DIN_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_no3_load)]);
saveas(gcf,strcat('DIN_load_vs_bay_chl_fr1997.png'),'png');

%% monthly
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_no3_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 1000]); xlim([1 length(sj.yr_no3_load)]);
ylabel('DIN load(g/s)');
plot(sj.mon_no3_load(i:12:end)+sj.mon_nh4_load(i:12:end),'k');
yyaxis right
plot(9:33,mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylim([0 15]);
ylabel('Chl.a(ug/L)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('DIN_load_vs_bay_chl_',num2str(i),'mon.png'),'png');
xlim([9 inf]);
saveas(gcf,strcat('DIN_load_vs_bay_chl_fr1997_',num2str(i),'mon.png'),'png');
end

%% Chl.a load vs chl
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_chl_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 2500]); 
xlim([1 length(sj.yr_chl_load)]);
ylabel('Chl.a load(mg/s)');
plot(sj.yr_chl_load,'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('chl_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_chl_load)]);
saveas(gcf,strcat('chl_load_vs_bay_chl_fr1997.png'),'png');

%% monthly
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_chl_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 2000]); xlim([1 length(sj.yr_chl_load)]);
ylabel('Chl.a load(mg/s)');
plot(sj.mon_chl_load(i:12:end),'k');
yyaxis right
plot(9:33,mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylim([0 15]);
ylabel('Chl.a(ug/L)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('chl_load_vs_bay_chl_',num2str(i),'mon.png'),'png');
xlim([9 inf]);
saveas(gcf,strcat('chl_load_vs_bay_chl_fr1997_',num2str(i),'mon.png'),'png');
end


%% SS load vs. chl
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_ss_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 2500]); 
xlim([1 length(sj.yr_ss_load)]);
ylabel('SS load(g/s)');
plot(sj.yr_ss_load,'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('SS_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_ss_load)]);
saveas(gcf,strcat('SS_load_vs_bay_chl_fr1997.png'),'png');

%% monthly
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_ss_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; ylim([0 2000]); xlim([1 length(sj.yr_ss_load)]);
ylabel('SS load(g/s)');
plot(sj.mon_ss_load(i:12:end),'k');
yyaxis right
plot(9:33,mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylim([0 15]);
ylabel('Chl.a(ug/L)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('SS_load_vs_bay_chl_',num2str(i),'mon.png'),'png');
xlim([9 inf]);
saveas(gcf,strcat('SS_load_vs_bay_chl_fr1997_',num2str(i),'mon.png'),'png');
end

gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 10]); xlim([3 inf]);
ylabel('Chl.a(ug/L)');
plot(mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_chl.png'),'png');

%% monthly
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 20]); xlim([3 inf]);
ylabel('Chl.a(ug/L)');
plot(mean(chl_s(164:165,i:12:end,sig),[1],'omitnan'),'c');
plot(mean(chl_s(193:195,i:12:end,sig),[1],'omitnan'),'k');
plot(mean(chl_s([191,192,189,188],i:12:end,sig),[1],'omitnan'),'b');
plot(mean(chl_s([204,205,207,206],i:12:end,sig),[1],'omitnan'),'m');
plot(mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_chl_2007_',num2str(i),'mon.png'),'png');
end

%%%

% gy_p_13 = [171 172 173 174 175 176 177 178 179 180 181 182 381];
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 10]); xlim([3 inf]);
ylabel('Chl.a(ug/L)');
plot(mean(chl_s_yr(164:165,:,sig),[1],'omitnan'),'c');
plot(mean(chl_s_yr(193:195,:,sig),[1],'omitnan'),'k');
plot(mean(chl_s_yr([191,192,189,188],:,sig),[1],'omitnan'),'b');
plot(mean(chl_s_yr([204,205,207,206],:,sig),[1],'omitnan'),'m');
plot(mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_chl.png'),'png');

xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_chl_2007.png'),'png');

plot(mean(chl_s_yr([238,237,236,234,235],:,sig),[1],'omitnan'),'color',[0.5 0.5 0.5]);

%% monthly
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 20]); xlim([3 inf]);
ylabel('Chl.a(ug/L)');
plot(mean(chl_s(164:165,i:12:end,sig),[1],'omitnan'),'c');
plot(mean(chl_s(193:195,i:12:end,sig),[1],'omitnan'),'k');
plot(mean(chl_s([191,192,189,188],i:12:end,sig),[1],'omitnan'),'b');
plot(mean(chl_s([204,205,207,206],i:12:end,sig),[1],'omitnan'),'m');
plot(mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_chl_2007_',num2str(i),'mon.png'),'png');
end

%% Chl.a
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 7]); xlim([3 inf]);
ylabel('Chl.a(ug/L)');
plot(mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
plot(mean(chl_b_s_yr(gy_p_9,:,sig),1,'omitnan'),'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_GY_only_chl.png'),'png');

xlim([11 inf]);

plot(11:19,repmat(mean(chl_s_yr(gy_p_9,11:19,sig),[1 2],'omitnan'),1,length(11:19)),'k--');
plot(20:25,repmat(mean(chl_s_yr(gy_p_9,20:25,sig),[1 2],'omitnan'),1,length(20:25)),'k--');
text(21,4,['차이 = ', num2str(mean(chl_s_yr(gy_p_9,11:19,sig),[1 2],'omitnan') - mean(chl_s_yr(gy_p_9,20:25,sig),[1 2],'omitnan'),'%0.2f'), 'ug/L']);
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_GY_only_chl_2007.png'),'png');


mean(chl_s_yr(gy_p_9,11:19,sig),[1 2],'omitnan') - mean(chl_s_yr(gy_p_9,20:25,sig),[1 2],'omitnan')


plot(mean(chl_s_yr([238,237,236,234,235],:,sig),[1],'omitnan'),'color',[0.5 0.5 0.5]);

%% monthly
sig=6;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 20]); xlim([3 inf]);
ylabel('Chl.a(ug/L)');
plot(mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
plot(mean(chl_b_s(gy_p_9,i:12:end,sig),1,'omitnan'),'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_GY_chl_',num2str(i),'mon.png'),'png');
% xlim([11 inf]);
% saveas(gcf,strcat('KOEM_compare_GY_chl_2007_',num2str(i),'mon.png'),'png');
end



%% 2016~2021 - 2007~2015

chl_s_diff_yr = squeeze(mean(chl_s(:,229:300,:) - mean(chl_s(:,121:228,:),[2],'omitnan'),[2],'omitnan'));

fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
plot(chl_s_diff_yr(:,3),'b.'); 
yline(mean(chl_s_diff_yr(:,3),'omitnan'),'r--','linew',2);
yline(mean(0,'omitnan'),'color',[0.5 0.5 0.5],'linew',2);
xlabel('station num.')
grid on; ylim([-15 3]);
xlim([1 length(chl_s_diff_yr)])
text(20,-13,'2016~2021 - 2007~2015 : Chl.a (ug/L)','fontsize',fontSizeTick+2,'fontweight','bold')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
saveas(gcf,strcat('KOEM_all_st_chl_period_diff.png'),'png');
close; 

clearvars chl_s_diff_mon
k=0
for i = [2, 5, 8, 11]
    k=k+1;
chl_s_diff_mon(:,:,k) = squeeze(mean(chl_s(:,229+(i-1):12:300,:) - mean(chl_s(:,121+(i-1):12:228,:),[2],'omitnan'),[2],'omitnan'));
end

mon=[2, 5, 8, 11];
for i = 1:4
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
plot(chl_s_diff_mon(:,3,i),'b.'); 
yline(mean(chl_s_diff_mon(:,3,i),'omitnan'),'r--','linew',2);
yline(mean(0,'omitnan'),'color',[0.5 0.5 0.5],'linew',2);
xlabel('station num.')
grid on; ylim([-24 8]);
xlim([1 length(chl_s_diff_yr)])
text(20,-20,['2016~2021 - 2007~2015 ',num2str(mon(i)),'월 : Chl.a (ug/L)'],'fontsize',fontSizeTick+2,'fontweight','bold')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
saveas(gcf,strcat('KOEM_all_st_chl_period_diff_',num2str(mon(i)),'mon.png'),'png');
% close; 
end

%% 2007~2015 obs repeats on each st. has to be bigger than "length(229:300)/4"
clearvars inout_ind
thresh_hold_t = length(229:300)/4;
for i = 1:465
    for sig_i = 1:6
        tempo_data = squeeze(chl_s(i,121:228,sig_i));
        if length(find(isnan(tempo_data) == 0)) >= thresh_hold_t;
            inout_ind(i,sig_i) = 0;
        elseif length(find(isnan(tempo_data) == 0)) < thresh_hold_t;
            inout_ind(i,sig_i) = NaN;
        end    
    end
end


fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
plot(chl_s_diff_yr(:,3) + inout_ind(:,3),'b.'); 
yline(mean(chl_s_diff_yr(:,3) + inout_ind(:,3),'omitnan'),'r--','linew',2);
yline(mean(0,'omitnan'),'color',[0.5 0.5 0.5],'linew',2);
xlabel('station num.')
grid on; ylim([-15 3]);
xlim([1 length(chl_s_diff_yr)])
text(20,-13,'2016~2021 - 2007~2015 : Chl.a (ug/L)','fontsize',fontSizeTick+2,'fontweight','bold')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
saveas(gcf,strcat('KOEM_part_st_chl_period_diff.png'),'png');
close; 


mon=[2, 5, 8, 11];
for i = 1:4
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
plot(chl_s_diff_mon(:,3,i)  + inout_ind(:,3),'b.'); 
yline(mean(chl_s_diff_mon(:,3,i) + inout_ind(:,3),'omitnan'),'r--','linew',2);
yline(mean(0,'omitnan'),'color',[0.5 0.5 0.5],'linew',2);
xlabel('station num.')
grid on; ylim([-24 8]);
xlim([1 length(chl_s_diff_yr)])
text(20,-20,['2016~2021 - 2007~2015 ',num2str(mon(i)),'월 : Chl.a (ug/L)'],'fontsize',fontSizeTick+2,'fontweight','bold')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
saveas(gcf,strcat('KOEM_part_st_chl_period_diff_',num2str(mon(i)),'mon.png'),'png');
% close; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spatial plot %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load D:\장기생태\Dynamic\KOEM\NIFS\Figure\jet_mod  % % set colormap (jet_modified)
run('D:\장기생태\Dynamic\KOEM\NIFS\Figure\bwr_map')  % % set colormap (blue and white)
colorbar_fontsize = 13;
colorbar_title_fontsize = 13;
colormap_style = jet_mod;  % % default
% trendlev = [-15 3];
trendlev = [-3 3];
% trendlev = [-24 8];
fontSizeTick = 13
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
geoshow(A,R); ylim([33 38.5])
xlim([125.5 130]);
title({'KOEM station'})
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on'); 

% 64
var_lin=linspace(trendlev(1),trendlev(2),100);
for i = 1:length(chl_s_diff_yr(:,3))
    if (isnan(chl_s_diff_yr(i,3)))==1
    else
    if (chl_s_diff_yr(i,3)<=min(trendlev))
       colind=1;
    elseif (chl_s_diff_yr(i,3)>=max(trendlev))
       colind=var_lin(end);
    else
       tempo_diff=(chl_s_diff_yr(i,3)-var_lin)
       colind=find(min(abs(tempo_diff)) == abs(tempo_diff));     
    end
%     colormap jet
%     cspec=jet;
%      h1=plot(lon_num(i),lat_num(i),'marker','o','color','k', ...
%                   'markerfacecolor',cspec(colind,:)); 
            h1=plot(lon_num(i),lat_num(i),'marker','o','color','k', ...
                  'markerfacecolor',bwrmap(colind,:),'MarkerSize',4); 
    end
end
 % set colorbar 
    h = colorbar;
    colormap(bwrmap);
    set(h,'fontsize',colorbar_fontsize);
    title(h,'Chl.a (ug/L)','fontsize',colorbar_title_fontsize);
    caxis(trendlev);
    saveas(gcf,strcat('KOEM_all_st_diff_spatial.png'),'png');

   
% month

colorbar_fontsize = 13;
colorbar_title_fontsize = 13;
colormap_style = jet_mod;  % % default
% trendlev = [-15 3];
trendlev = [-3 3];
% trendlev = [-24 8];
fontSizeTick = 13

% 64
var_lin=linspace(trendlev(1),trendlev(2),100);
mon=[2, 5, 8, 11];
for j = 1:4
    figPos = [0 0 5 4];
    figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
    geoshow(A,R); ylim([33 38.5])
    xlim([125.5 130]);
    title({'KOEM station'})
    set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on'); 
    for i = 1:length(chl_s_diff_mon(:,3,j))
        if (isnan(chl_s_diff_mon(i,3,j)))==1
        else
        if (chl_s_diff_mon(i,3,j)<=min(trendlev))
           colind=1;
        elseif (chl_s_diff_mon(i,3,j)>=max(trendlev))
           colind=var_lin(end);
        else
           tempo_diff=(chl_s_diff_mon(i,3,j)-var_lin)
           colind=find(min(abs(tempo_diff)) == abs(tempo_diff));     
        end
    %     colormap jet
    %     cspec=jet;
    %      h1=plot(lon_num(i),lat_num(i),'marker','o','color','k', ...
    %                   'markerfacecolor',cspec(colind,:)); 
            h1=plot(lon_num(i),lat_num(i),'marker','o','color','k', ...
                  'markerfacecolor',bwrmap(colind,:),'MarkerSize',4); 
        end
    end
    
    % set colorbar 
%     h = colorbar;
%     colormap(bwrmap);
%     set(h,'fontsize',colorbar_fontsize);
%     title(h,'Chl.a (ug/L)','fontsize',colorbar_title_fontsize);
%     caxis(trendlev);
    text(127, 37,['KOEM ',num2str(mon(j)),'월'],'color','w','fontsize',18)
    saveas(gcf,strcat('KOEM_all_st_diff_spatial_',num2str(mon(j)),'mon.png'),'png');
    hold off;
end
 

%% DIN

gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 300]); xlim([11 inf]);
ylabel('DIN(mmol/m^3)');
plot(mean(DIN_s_yr(164:165,:,sig),[1],'omitnan'),'c');
plot(mean(DIN_s_yr(193:195,:,sig),[1],'omitnan'),'k');
plot(mean(DIN_s_yr([191,192,189,188],:,sig),[1],'omitnan'),'b');
plot(mean(DIN_s_yr([204,205,207,206],:,sig),[1],'omitnan'),'m');
plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_DIN.png'),'png');

gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 1000]); xlim([11 inf]);
ylabel('TN (mmol/m^3)');
plot(mean(TN_s_yr(164:165,:,sig),[1],'omitnan'),'c');
plot(mean(TN_s_yr(193:195,:,sig),[1],'omitnan'),'k');
plot(mean(TN_s_yr([191,192,189,188],:,sig),[1],'omitnan'),'b');
plot(mean(TN_s_yr([204,205,207,206],:,sig),[1],'omitnan'),'m');
plot(mean(TN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on'); 
saveas(gcf,strcat('KOEM_compare_bay_TN.png'),'png');
xlim([3 inf]);
saveas(gcf,strcat('KOEM_compare_bay_TN_1999.png'),'png');


gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 2000]); xlim([11 inf]);
ylabel('TN (mmol/m^3)');
for i = 1:465
plot(TN_s_yr(i,:,sig),'.','color',[0.5 0.5 0.5]);
end
% plot(mean(TN_s_yr(193:195,:,sig),[1],'omitnan'),'k');
% plot(mean(TN_s_yr([191,192,189,188],:,sig),[1],'omitnan'),'b');
% plot(mean(TN_s_yr([204,205,207,206],:,sig),[1],'omitnan'),'m');
plot(mean(TN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_all_st_TN.png'),'png');
xlim([3 inf]);
saveas(gcf,strcat('KOEM_compare_all_st_TN_1999.png'),'png');

sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 1000]); xlim([3 inf]);
ylabel('DIP (mmol/m^3)');
plot(mean(DIP_s(164:165,i:12:end,sig),[1],'omitnan'),'c');
plot(mean(DIP_s(193:195,i:12:end,sig),[1],'omitnan'),'k');
plot(mean(DIP_s([191,192,189,188],i:12:end,sig),[1],'omitnan'),'b');
plot(mean(DIP_s([204,205,207,206],i:12:end,sig),[1],'omitnan'),'m');
plot(mean(DIP_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_DIP_2007_',num2str(i),'mon.png'),'png');
end

sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 1000]); xlim([3 inf]);
ylabel('TN (mmol/m^3)');
plot(mean(TN_s(164:165,i:12:end,sig),[1],'omitnan'),'c');
plot(mean(TN_s(193:195,i:12:end,sig),[1],'omitnan'),'k');
plot(mean(TN_s([191,192,189,188],i:12:end,sig),[1],'omitnan'),'b');
plot(mean(TN_s([204,205,207,206],i:12:end,sig),[1],'omitnan'),'m');
plot(mean(TN_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_TN_2007_',num2str(i),'mon.png'),'png');
xlim([3 inf]);
saveas(gcf,strcat('KOEM_compare_bay_TN_1999_',num2str(i),'mon.png'),'png');
end


gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 100]); xlim([11 inf]);
ylabel('TP (mmol/m^3)');
plot(mean(TP_s_yr(164:165,:,sig),[1],'omitnan'),'c');
plot(mean(TP_s_yr(193:195,:,sig),[1],'omitnan'),'k');
plot(mean(TP_s_yr([191,192,189,188],:,sig),[1],'omitnan'),'b');
plot(mean(TP_s_yr([204,205,207,206],:,sig),[1],'omitnan'),'m');
plot(mean(TP_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_TP.png'),'png');



sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 100]); xlim([3 inf]);
ylabel('TP (mmol/m^3)');
plot(mean(TP_s(164:165,i:12:end,sig),[1],'omitnan'),'c');
plot(mean(TP_s(193:195,i:12:end,sig),[1],'omitnan'),'k');
plot(mean(TP_s([191,192,189,188],i:12:end,sig),[1],'omitnan'),'b');
plot(mean(TP_s([204,205,207,206],i:12:end,sig),[1],'omitnan'),'m');
plot(mean(TP_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_TP_2007_',num2str(i),'mon.png'),'png');
end

gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 1200]); xlim([11 inf]);
ylabel('si (mmol/m^3)');
plot(mean(si_s_yr(164:165,:,sig),[1],'omitnan'),'c');
plot(mean(si_s_yr(193:195,:,sig),[1],'omitnan'),'k');
plot(mean(si_s_yr([191,192,189,188],:,sig),[1],'omitnan'),'b');
plot(mean(si_s_yr([204,205,207,206],:,sig),[1],'omitnan'),'m');
plot(mean(si_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_si.png'),'png');



sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 1200]); xlim([3 inf]);
ylabel('si (mmol/m^3)');
plot(mean(si_s(164:165,i:12:end,sig),[1],'omitnan'),'c');
plot(mean(si_s(193:195,i:12:end,sig),[1],'omitnan'),'k');
plot(mean(si_s([191,192,189,188],i:12:end,sig),[1],'omitnan'),'b');
plot(mean(si_s([204,205,207,206],i:12:end,sig),[1],'omitnan'),'m');
plot(mean(si_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_si_2007_',num2str(i),'mon.png'),'png');
end


gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 5]); xlim([11 inf]);
ylabel('secchi(m)');
plot(mean(secchi_s_yr(164:165,:,sig),[1],'omitnan'),'c');
plot(mean(secchi_s_yr(193:195,:,sig),[1],'omitnan'),'k');
plot(mean(secchi_s_yr([191,192,189,188],:,sig),[1],'omitnan'),'b');
plot(mean(secchi_s_yr([204,205,207,206],:,sig),[1],'omitnan'),'m');
plot(mean(secchi_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_secchi.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('KOEM_compare_bay_secchi_all.png'),'png');

sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 5]); xlim([3 inf]);
ylabel('secchi(m)');
plot(mean(secchi_s(164:165,i:12:end,sig),[1],'omitnan'),'c');
plot(mean(secchi_s(193:195,i:12:end,sig),[1],'omitnan'),'k');
plot(mean(secchi_s([191,192,189,188],i:12:end,sig),[1],'omitnan'),'b');
plot(mean(secchi_s([204,205,207,206],i:12:end,sig),[1],'omitnan'),'m');
plot(mean(secchi_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_secchi_2007_',num2str(i),'mon.png'),'png');
end

%% TP
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 60]); xlim([11 inf]);
ylabel('TP(mmol/m^3)');
plot(mean(TP_s_yr(164:165,:,sig),[1],'omitnan'),'c');
plot(mean(TP_s_yr(193:195,:,sig),[1],'omitnan'),'k');
plot(mean(TP_s_yr([191,192,189,188],:,sig),[1],'omitnan'),'b');
plot(mean(TP_s_yr([204,205,207,206],:,sig),[1],'omitnan'),'m');
plot(mean(TP_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_TP.png'),'png');


sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 81]); xlim([3 inf]);
ylabel('TP(mmol/m^3)');
plot(mean(TP_s(164:165,i:12:end,sig),[1],'omitnan'),'c');
plot(mean(TP_s(193:195,i:12:end,sig),[1],'omitnan'),'k');
plot(mean(TP_s([191,192,189,188],i:12:end,sig),[1],'omitnan'),'b');
plot(mean(TP_s([204,205,207,206],i:12:end,sig),[1],'omitnan'),'m');
plot(mean(TP_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_TP_2007_',num2str(i),'mon.png'),'png');
end


%% DIP
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 50]); xlim([11 inf]);
ylabel('DIP(mmol/m^3)');
plot(mean(DIP_s_yr(164:165,:,sig),[1],'omitnan'),'c');
plot(mean(DIP_s_yr(193:195,:,sig),[1],'omitnan'),'k');
plot(mean(DIP_s_yr([191,192,189,188],:,sig),[1],'omitnan'),'b');
plot(mean(DIP_s_yr([204,205,207,206],:,sig),[1],'omitnan'),'m');
plot(mean(DIP_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_DIP.png'),'png');


sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 50]); xlim([3 inf]);
ylabel('DIP(mmol/m^3)');
plot(mean(DIP_s(164:165,i:12:end,sig),[1],'omitnan'),'c');
plot(mean(DIP_s(193:195,i:12:end,sig),[1],'omitnan'),'k');
plot(mean(DIP_s([191,192,189,188],i:12:end,sig),[1],'omitnan'),'b');
plot(mean(DIP_s([204,205,207,206],i:12:end,sig),[1],'omitnan'),'m');
plot(mean(DIP_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_DIP_2007_',num2str(i),'mon.png'),'png');
end

%% N:P ratio
MW_N = 14.006720;
MW_P = 30.973762;

gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 120]); xlim([11 inf]);
ylabel('DIN/DIP ratio');
plot((mean(DIN_s_yr(164:165,:,sig),[1],'omitnan')./MW_N)./(mean(DIP_s_yr(164:165,:,sig),[1],'omitnan')./MW_P),'c');
plot((mean(DIN_s_yr(193:195,:,sig),[1],'omitnan')./MW_N)./(mean(DIP_s_yr(193:195,:,sig),[1],'omitnan')./MW_P),'k');
plot((mean(DIN_s_yr([191,192,189,188],:,sig),[1],'omitnan')./MW_N)./(mean(DIP_s_yr([191,192,189,188],:,sig),[1],'omitnan')./MW_P),'b');
plot((mean(DIN_s_yr([204,205,207,206],:,sig),[1],'omitnan')./MW_N)./(mean(DIP_s_yr([204,205,207,206],:,sig),[1],'omitnan')./MW_P),'m');
plot((mean(DIN_s_yr(gy_p_9,:,sig),[1],'omitnan')./MW_N)./(mean(DIP_s_yr(gy_p_9,:,sig),1,'omitnan')./MW_P),'r');
yline(16,'--','color',[0.5 0.5 0.5],'linew',2);
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_NP.png'),'png');


sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 120]); xlim([3 inf]);
ylabel('DIN/DIP ratio');
plot((mean(DIN_s(164:165,i:12:end,sig),[1],'omitnan')./MW_N)./(mean(DIP_s(164:165,i:12:end,sig),[1],'omitnan')./MW_P),'c');
plot((mean(DIN_s(193:195,i:12:end,sig),[1],'omitnan')./MW_N)./(mean(DIP_s(193:195,i:12:end,sig),[1],'omitnan')./MW_P),'k');
plot((mean(DIN_s([191,192,189,188],i:12:end,sig),[1],'omitnan')./MW_N)./(mean(DIP_s([191,192,189,188],i:12:end,sig),[1],'omitnan')./MW_P),'b');
plot((mean(DIN_s([204,205,207,206],i:12:end,sig),[1],'omitnan')./MW_N)./(mean(DIP_s([204,205,207,206],i:12:end,sig),[1],'omitnan')./MW_P),'m');
plot((mean(DIN_s(gy_p_9,i:12:end,sig),[1],'omitnan')./MW_N)./(mean(DIP_s(gy_p_9,i:12:end,sig),1,'omitnan')./MW_P),'r');
yline(16,'--','color',[0.5 0.5 0.5],'linew',2);
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_NP_2007_',num2str(i),'mon.png'),'png');
end

sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 50]); xlim([11 inf]);
ylabel('TN/TP ratio');
plot(mean(TN_s_yr(164:165,:,sig),[1],'omitnan')./mean(TP_s_yr(164:165,:,sig),[1],'omitnan'),'c');
plot(mean(TN_s_yr(193:195,:,sig),[1],'omitnan')./mean(TP_s_yr(193:195,:,sig),[1],'omitnan'),'k');
plot(mean(TN_s_yr([191,192,189,188],:,sig),[1],'omitnan')./mean(TP_s_yr([191,192,189,188],:,sig),[1],'omitnan'),'b');
plot(mean(TN_s_yr([204,205,207,206],:,sig),[1],'omitnan')./mean(TP_s_yr([204,205,207,206],:,sig),[1],'omitnan'),'m');
plot(mean(TN_s_yr(gy_p_9,:,sig),[1],'omitnan')./mean(TP_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
yline(16,'--','color',[0.5 0.5 0.5],'linew',2);
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('KOEM_compare_bay_TNP.png'),'png');


sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 50]); xlim([3 inf]);
ylabel('TN/TP ratio');
plot(mean(TN_s(164:165,i:12:end,sig),[1],'omitnan')./mean(TP_s(164:165,i:12:end,sig),[1],'omitnan'),'c');
plot(mean(TN_s(193:195,i:12:end,sig),[1],'omitnan')./mean(TP_s(193:195,i:12:end,sig),[1],'omitnan'),'k');
plot(mean(TN_s([191,192,189,188],i:12:end,sig),[1],'omitnan')./mean(TP_s([191,192,189,188],i:12:end,sig),[1],'omitnan'),'b');
plot(mean(TN_s([204,205,207,206],i:12:end,sig),[1],'omitnan')./mean(TP_s([204,205,207,206],i:12:end,sig),[1],'omitnan'),'m');
plot(mean(TN_s(gy_p_9,i:12:end,sig),[1],'omitnan')./mean(TP_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
yline(16,'--','color',[0.5 0.5 0.5],'linew',2);
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([11 inf]);
saveas(gcf,strcat('KOEM_compare_bay_TNP_2007_',num2str(i),'mon.png'),'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GY only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% KD
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 1.5]); xlim([11 inf]);
ylabel('Kd (1/m)');
plot(1.7./mean(secchi_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('KOEM_compare_kd_all.png'),'png');



sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 1.5]); %xlim([3 inf]);
ylabel('Kd (1/m)');
plot(1.7./mean(secchi_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('KOEM_compare_bay_kd_all_',num2str(i),'mon.png'),'png');
end

%% temp
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([12 18]); xlim([11 inf]);
ylabel('temp (^oC)');
plot(mean(temp_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
plot(mean(temp_b_s_yr(gy_p_9,:,sig),1,'omitnan'),'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('KOEM_compare_temp_all.png'),'png');



sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 30]); %xlim([3 inf]);
ylabel('temp (^oC)');
plot(mean(temp_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
plot(mean(temp_b_s(gy_p_9,i:12:end,sig),1,'omitnan'),'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('KOEM_compare_bay_temp_all_',num2str(i),'mon.png'),'png');
end

%% salt
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([24 34]); xlim([11 inf]);
ylabel('salt (g/kg)');
plot(mean(salt_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
plot(mean(salt_b_s_yr(gy_p_9,:,sig),1,'omitnan'),'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('KOEM_compare_salt_all.png'),'png');



sig=3;
fontSizeTick = 11
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([24 34]); %xlim([3 inf]);
ylabel('salt (g/kg)');
plot(mean(salt_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
plot(mean(salt_b_s(gy_p_9,i:12:end,sig),1,'omitnan'),'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('KOEM_compare_bay_salt_all_',num2str(i),'mon.png'),'png');
end


%% air 
load('D:\장기생태\Dynamic\KOEM\ERA5_forc_daily_1980to2021_spmean.mat');

path_xls = 'D:\장기생태\Dynamic\KOEM\monthly_ASOS_&_AWS\';

[raw_aws_yeosu txt_aws_yeosu]=xlsread([path_xls,'OBS_AWS_여수산단_월평균.xlsx'],'OBS_AWS_여수산단_월평균','');
[raw_asos_yeosu txt_asos_yeosu]=xlsread([path_xls,'OBS_ASOS_여수_월평균.xlsx'],'OBS_ASOS_여수_월평균','');
% 4 col : mean temp, 11 col : wind speed, 15 col : monthly culmulative rain

yeosu_aws_temp=raw_aws_yeosu(:,4);
yeosu_aws_mintemp=raw_aws_yeosu(:,6);
yeosu_aws_speed=raw_aws_yeosu(:,11);
yeosu_aws_rain=raw_aws_yeosu(:,15);
yeosu_aws_maxspeed=raw_aws_yeosu(:,12);
yeosu_aws_maxspeed_dir=raw_aws_yeosu(:,13);

%% make 1997 to 2019 'yymm' form
k=0
for i = 1942:2021
    for j = 1:12
        k=k+1;
        ref_date{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end 

clearvars date_asos_yeosu_c
date_asos_yeosu=char(txt_asos_yeosu(2:end,3));
for i=1:length(date_asos_yeosu)
date_asos_yeosu_c{i,1}= date_asos_yeosu(i,1:7);
end

for i = 1:length(ref_date)
    check_size=[];
    check_size=length(find(strcmp(ref_date{i},date_asos_yeosu_c)==1));
    
    if check_size >= 1
       idx_asos(i)=find(strcmp(ref_date{i},date_asos_yeosu_c)==1);
    else
        idx_asos(i)=NaN;
    end
end

for i = 1:length(idx_asos)
    if isnan(idx_asos(i)) == 0
    yeosu_asos_temp(i)=raw_asos_yeosu(idx_asos(i),4);
    yeosu_asos_pair(i)=raw_asos_yeosu(idx_asos(i),12);
    yeosu_asos_speed(i)=raw_asos_yeosu(idx_asos(i),39);
    yeosu_asos_rain(i)=raw_asos_yeosu(idx_asos(i),26);
    yeosu_asos_evap_s(i)=raw_asos_yeosu(idx_asos(i),33);
    yeosu_asos_evap_l(i)=raw_asos_yeosu(idx_asos(i),35);
    yeosu_asos_srad(i)=raw_asos_yeosu(idx_asos(i),51);
    yeosu_asos_soil_t(i)=raw_asos_yeosu(idx_asos(i),60);
    elseif isnan(idx_asos(i)) == 1
    yeosu_asos_temp(i)=NaN;
    yeosu_asos_pair(i)=NaN;
    yeosu_asos_speed(i)=NaN;
    yeosu_asos_rain(i)=NaN;
    yeosu_asos_evap_s(i)=NaN;
    yeosu_asos_evap_l(i)=NaN;
    yeosu_asos_srad(i)=NaN;
    yeosu_asos_soil_t(i)=NaN;
    end
end

clearvars yeosu_aws_*_yr
for i=1:300/12
yeosu_aws_temp_yr(i)=mean(yeosu_aws_temp((i-1)*12+1:i*12),[1],'omitnan');
yeosu_aws_mintemp_yr(i)=mean(yeosu_aws_mintemp((i-1)*12+1:i*12),[1],'omitnan');
yeosu_aws_speed_yr(i)=mean(yeosu_aws_speed((i-1)*12+1:i*12),[1],'omitnan');
yeosu_aws_rain_yr(i)=sum(yeosu_aws_rain((i-1)*12+1:i*12),[1],'omitnan');
yeosu_aws_maxspeed_dir_yr(i)=mean(yeosu_aws_maxspeed_dir((i-1)*12+1:i*12),[1],'omitnan');
yeosu_aws_maxspeed_yr(i)=mean(yeosu_aws_maxspeed((i-1)*12+1:i*12),[1],'omitnan');
end


clearvars yeosu_asos_*_yr
for i=1:length(idx_asos)/12
yeosu_asos_temp_yr(i)=mean(yeosu_asos_temp((i-1)*12+1:i*12),[2],'omitnan');
yeosu_asos_pair_yr(i)=mean(yeosu_asos_pair((i-1)*12+1:i*12),[2],'omitnan');
yeosu_asos_speed_yr(i)=mean(yeosu_asos_speed((i-1)*12+1:i*12),[2],'omitnan');
yeosu_asos_rain_yr(i)=sum(yeosu_asos_rain((i-1)*12+1:i*12),[2],'omitnan');
yeosu_asos_evap_s_yr(i)=sum(yeosu_asos_evap_s((i-1)*12+1:i*12),[2],'omitnan');
yeosu_asos_evap_l_yr(i)=sum(yeosu_asos_evap_l((i-1)*12+1:i*12),[2],'omitnan');
yeosu_asos_srad_yr(i)=mean(yeosu_asos_srad((i-1)*12+1:i*12),[2],'omitnan');
yeosu_asos_soil_t_yr(i)=mean(yeosu_asos_soil_t((i-1)*12+1:i*12),[2],'omitnan');
end


% airT
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([10 18]); xlim([11 inf]);
ylabel('air temp (^oC)');
plot(yeosu_aws_temp_yr,'r');
plot(yeosu_aws_mintemp_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_temp_all.png'),'png');

fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([14 18]); xlim([11 inf]);
ylabel('air temp (^oC)');
plot(yeosu_aws_temp_yr,'r');
% plot(yeosu_aws_mintemp_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_temp_yr.png'),'png');

fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; ylim([9.5 16.5]); 
ylabel('air temp (^oC)');
plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp_yr,'r');
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_temp_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_temp_asos_vs_aws_yr_full.png'),'png');

xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; ylim([13.5 16.5]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_temp_asos_vs_aws_yr.png'),'png');


figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 30]); %xlim([3 inf]);
ylabel('temp (^oC)');
plot(yeosu_aws_temp(i:12:end),'r-')
% plot(mean(temp_b_s(gy_p_9,i:12:end,sig),1,'omitnan'),'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_temp_',num2str(i),'mon.png'),'png');
end


figPos = [0 0 5 4];
% for i = [2,5,8,11]
for i = [1:12]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; ylim([0 30]); %xlim([3 inf]);
ylabel('temp (^oC)');
plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp(i:12:end),'r-')
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_temp(i:12:end),'b-')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_temp_asos_vs_aws_',num2str(i),'mon_full.png'),'png');
xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; ylim([0 30]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_temp_asos_vs_aws_',num2str(i),'mon.png'),'png');
end

% pair
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; %ylim([9.5 16.5]); 
ylabel('air pressure (hPa)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp_yr,'r');
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_pair_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_pair_asos_yr_full.png'),'png');

xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; ylim([1015 1017]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_pair_asos_yr.png'),'png');


figPos = [0 0 5 4];
% for i = [2,5,8,11]
for i = [1:12]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; ylim([1005 1027]); %xlim([3 inf]);
ylabel('air pressure (hPa)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp(i:12:end),'r-')
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_pair(i:12:end),'b-')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_pair_asos_',num2str(i),'mon_full.png'),'png');
xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([0 30]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_pair_asos_',num2str(i),'mon.png'),'png');
end

% speed
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([2 20]); xlim([11 inf]);
ylabel('wind speend (m/s)');
plot(yeosu_aws_speed_yr,'r');
plot(yeosu_aws_maxspeed_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_speed_all.png'),'png');

fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; %ylim([2 20]); 
xlim([11 inf]);
ylabel('wind direction (^o)');
plot(yeosu_aws_maxspeed_dir_yr,'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_wind_dir.png'),'png');

fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 5]); xlim([11 inf]);
ylabel('wind speend (m/s)');
plot(yeosu_aws_speed_yr,'r');
% plot(yeosu_aws_mintemp_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_speed_yr.png'),'png');


fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on;  xlim([11 inf]);
ylabel('wind speend (m/s)');
plot(yeosu_aws_maxspeed_yr,'r');
% plot(yeosu_aws_mintemp_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_maxspeed_yr.png'),'png');


figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([10 30]); %xlim([3 inf]);
ylabel('wind speend (m/s)');
plot(yeosu_aws_maxspeed(i:12:end),'r-')
% plot(mean(temp_b_s(gy_p_9,i:12:end,sig),1,'omitnan'),'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_maxspeed_',num2str(i),'mon.png'),'png');
end


fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; %ylim([9.5 16.5]); 
ylabel('wind speed (m/s)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp_yr,'r');
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_speed_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_speed_asos_yr_full.png'),'png');

xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([1015 1017]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_speed_asos_yr.png'),'png');


figPos = [0 0 5 4];
% for i = [2,5,8,11]
for i = [1:12]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; ylim([2.5 6]); %xlim([3 inf]);
ylabel('wind speed (m/s)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp(i:12:end),'r-')
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_speed(i:12:end),'b-')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_speed_asos_',num2str(i),'mon_full.png'),'png');
xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([0 30]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_speed_asos_',num2str(i),'mon.png'),'png');
end

%% dir
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on;  ylim([45 360]);
ylabel('wind direction (^o)');
plot(yeosu_aws_maxspeed_dir_yr,'r');
% plot(yeosu_aws_mintemp_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_maxspeed_dir_yr.png'),'png');


figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([45 360]); %xlim([3 inf]);
ylabel('wind direction (^o)');
plot(yeosu_aws_maxspeed_dir(i:12:end),'r-')
% plot(mean(temp_b_s(gy_p_9,i:12:end,sig),1,'omitnan'),'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_maxspeed_dir_',num2str(i),'mon.png'),'png');
end

%% rain
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on;  ylim([0 2300]);
ylabel('rain (mm)');
plot(yeosu_aws_rain_yr,'r');
% plot(yeosu_aws_mintemp_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_rain_yr.png'),'png');


figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:300/12)
xticklabels(1997:2021)
xtickangle(45)
grid on; ylim([0 800]); %xlim([3 inf]);
ylabel('rain (mm)');
plot(yeosu_aws_rain(i:12:end),'r-')
% plot(mean(temp_b_s(gy_p_9,i:12:end,sig),1,'omitnan'),'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_rain_',num2str(i),'mon.png'),'png');
end


fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; %ylim(500 2500]); 
ylabel('rain (mm)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp_yr,'r');
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_rain_yr,'b','linew',2);
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_rain_asos_yr_full.png'),'png');

xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([1015 1017]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_rain_asos_yr.png'),'png');
yyaxis right 
plot(length(1942:1997):length(1942:2021),sj.yr_trans(9:33),'r','linew',2);
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
ylabel('discharge (m^3/s)');
saveas(gcf,strcat('yeosu_aws_rain_asos_vs_dis_yr.png'),'png');
plot(length(1942:1997):length(1942:2021),sj.yr_no3_load(9:33)+sj.yr_nh4_load(9:33),'g-','linew',2);
plot(length(1942:1997):length(1942:2021),sj.yr_tn_load(9:33),'k-','linew',2);
saveas(gcf,strcat('yeosu_aws_rain_asos_vs_n_load_yr.png'),'png');
legend('강우','방류량','DIN','TN','NumColumns',2);

% plot(length(1942:1989):length(1942:2021),sj.yr_trans(1:33),'r','linew',2);
% plot(length(1942:1989):length(1942:2021),sj.yr_no3_load(1:33)+sj.yr_nh4_load(1:33),'g-','linew',2);
% plot(length(1942:1989):length(1942:2021),sj.yr_tn_load(1:33),'k-','linew',2);
% xlim([length(1942:1989) inf]);
% xtickangle(90)
plot(length(1942:1997):length(1942:2021),sj.yr_trans(9:33),'r','linew',2);



figPos = [0 0 5 4];
% for i = [2,5,8,11]
for i = [1:12]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; ylim([0 800]); %xlim([3 inf]);
ylabel('rain (mm)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp(i:12:end),'r-')
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_rain(i:12:end),'b-','linew',2)
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_rain_asos_',num2str(i),'mon_full.png'),'png');
xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([0 30]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_rain_asos_',num2str(i),'mon.png'),'png');
ylim([-inf inf]);
yyaxis right 
plot(length(1942:1989):length(1942:2021),sj.monthly_trans(i:12:end),'r','linew',2);
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
ylabel('discharge (m^3/s)');
saveas(gcf,strcat('yeosu_aws_rain_asos_vs_dis_',num2str(i),'mon.png'),'png');
plot(length(1942:1989):length(1942:2021),sj.mon_no3_load(i:12:end)+sj.mon_nh4_load(i:12:end),'g-','linew',2);
plot(length(1942:1989):length(1942:2021),sj.mon_tn_load(i:12:end),'k-','linew',2);
saveas(gcf,strcat('yeosu_aws_rain_asos_vs_n_load_',num2str(i),'mon.png'),'png');
% legend('강우','방류량','DIN','TN','NumColumns',2);
end


% evaporation
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; %ylim(500 2500]); 
ylabel('evaporation (mm)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp_yr,'r');
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_evap_s_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_evap_s_asos_yr_full.png'),'png');

xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([1015 1017]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_evap_s_asos_yr.png'),'png');


figPos = [0 0 5 4];
% for i = [2,5,8,11]
for i = [1:12]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; ylim([0 300]); %xlim([3 inf]);
ylabel('evaporation (mm)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp(i:12:end),'r-')
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_evap_s(i:12:end),'b-')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_evap_s_asos_',num2str(i),'mon_full.png'),'png');
xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([0 30]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_evap_s_asos_',num2str(i),'mon.png'),'png');
end


fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; %ylim(500 2500]); 
ylabel('P-E (mm)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp_yr,'r');
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_rain_yr-yeosu_asos_evap_s_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_P-E_asos_yr_full.png'),'png');

xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([1015 1017]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_P-E_asos_yr.png'),'png');
% yyaxis right 
% plot(length(1942:1997):length(1942:2021),sj.yr_trans(9:33),'r','linew',2);

figPos = [0 0 5 4];
% for i = [2,5,8,11]
for i = [1:12]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; ylim([-600 600]); %xlim([3 inf]);
ylabel('P-E (mm)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp(i:12:end),'r-')
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_rain(i:12:end)-yeosu_asos_evap_s(i:12:end),'b-')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_P-E_asos_',num2str(i),'mon_full.png'),'png');
xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([0 30]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_P-E_asos_',num2str(i),'mon.png'),'png');
end

fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; %ylim(500 2500]); 
ylabel('evaporation (mm)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp_yr,'r');
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_evap_l_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_evap_l_asos_yr_full.png'),'png');

xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([1015 1017]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_evap_l_asos_yr.png'),'png');


figPos = [0 0 5 4];
% for i = [2,5,8,11]
for i = [1:12]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; ylim([0 300]); %xlim([3 inf]);
ylabel('evaporation (mm)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp(i:12:end),'r-')
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_evap_l(i:12:end),'b-')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_evap_l_asos_',num2str(i),'mon_full.png'),'png');
xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([0 30]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_evap_l_asos_',num2str(i),'mon.png'),'png');
end

% radiation
fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; %ylim(500 2500]); 
ylabel('radiation (MJ/m^2)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp_yr,'r');
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_soil_t_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_srad_asos_yr_full.png'),'png');

xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([1015 1017]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_srad_asos_yr.png'),'png');


figPos = [0 0 5 4];
% for i = [2,5,8,11]
for i = [1:12]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; ylim([0 650]); %xlim([3 inf]);
ylabel('radiation (MJ/m^2)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp(i:12:end),'r-')
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_srad(i:12:end),'b-')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_srad_asos_',num2str(i),'mon_full.png'),'png');
xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([0 30]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_srad_asos_',num2str(i),'mon.png'),'png');
end

% soil temperature

fontSizeTick = 11
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; %ylim(500 2500]); 
ylabel('soil temperature (^oC)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp_yr,'r');
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_soil_t_yr,'b');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_soil_t_asos_yr_full.png'),'png');

xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([1015 1017]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_soil_t_asos_yr.png'),'png');


figPos = [0 0 5 4];
% for i = [2,5,8,11]
for i = [1:12]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:4:length(yeosu_asos_temp_yr))
xticklabels(1942:4:2021)
xtickangle(45)
grid on; ylim([0 32]); %xlim([3 inf]);
ylabel('soil temperature (^oC)');
% plot(length(1942:1997):length(yeosu_asos_temp_yr),yeosu_aws_temp(i:12:end),'r-')
plot(1:length(yeosu_asos_temp_yr),yeosu_asos_soil_t(i:12:end),'b-')
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('KOEM_compare_bay_chl_',num2str(i),'mon.png'),'png');
xlim([1 inf]);
saveas(gcf,strcat('yeosu_aws_soil_t_asos_',num2str(i),'mon_full.png'),'png');
xticks(1:1:length(yeosu_asos_temp_yr))
xticklabels(1942:1:2021)
xtickangle(45)
grid on; %ylim([0 30]); 
xlim([length(1942:1997) inf]);
saveas(gcf,strcat('yeosu_aws_soil_t_asos_',num2str(i),'mon.png'),'png');
end

%% chl vs. Srad
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sp_sw_yr))
xticklabels(1980:2021)
xtickangle(90)
grid on;  xlim([1 length(sp_sw_yr)]);
% ylim([0 20]);
ylabel('radiation (W/m^2)');
plot(sp_sw_yr,'k');
yyaxis right
plot(18:length(sp_sw_yr),mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 7]);
xlim([18 length(sp_sw_yr)]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
saveas(gcf,strcat('SW_vs_bay_chl_yr_fr1997.png'),'png');


%% monthly
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sp_sw_yr))
xticklabels(1980:2021)
xtickangle(90)
grid on; 
% ylim([130 180]); 
xlim([1 length(sp_sw_yr)]);
ylabel('radiation (W/m^2)');
plot(sp_sw_mon(i,:),'k');
yyaxis right
plot(18:length(sp_sw_yr),mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylim([0 15]);
ylabel('Chl.a(ug/L)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
xlim([18 length(sp_sw_yr)]);
saveas(gcf,strcat('SW_vs_bay_chl_fr1997_',num2str(i),'mon.png'),'png');
end


%% NH4 vs. chl
gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=6;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(chl_s_yr))
xticklabels(1989:2021)
xtickangle(90)
grid on;% ylim([0 50]); 
xlim([1 length(sj.yr_tn_load)]);
ylabel('NH4-N (mmol/L)');
plot(9:33,mean(nh4_s_yr(gy_p_9,:,sig),1,'omitnan'),'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('TN_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_tp_load)]);
% saveas(gcf,strcat('TN_load_vs_bay_chl_fr1997.png'),'png');

gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=6;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(chl_s_yr))
xticklabels(1989:2021)
xtickangle(90)
grid on;% ylim([0 50]); 
xlim([1 length(sj.yr_tn_load)]);
ylabel('NH4/DIN (mmol/L)');
plot(9:33,mean(nh4_s_yr(gy_p_9,:,sig),1,'omitnan')./(mean(nh4_s_yr(gy_p_9,:,sig),1,'omitnan')+mean(no3_s_yr(gy_p_9,:,sig),1,'omitnan')),'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('TN_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_tp_load)]);


gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=6;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(chl_s_yr))
xticklabels(1989:2021)
xtickangle(90)
grid on;% ylim([0 50]); 
xlim([1 length(sj.yr_tn_load)]);
ylabel('NH4/TN (mmol/L)');
plot(9:33,mean(nh4_s_yr(gy_p_9,:,sig),1,'omitnan')./(mean(TN_s_yr(gy_p_9,:,sig),1,'omitnan')),'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('TN_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_tp_load)]);


gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(chl_s_yr))
xticklabels(1989:2021)
xtickangle(90)
grid on;% ylim([0 50]); 
xlim([1 length(sj.yr_tn_load)]);
ylabel('NH4/(TN-DIN) (mmol/L)');
plot(9:33,mean(nh4_s_yr(gy_p_9,:,sig),1,'omitnan')./(mean(TN_s_yr(gy_p_9,:,sig),1,'omitnan') - mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan')),'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('TN_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_tp_load)]);


sig=3
NH4_ov_DON=mean(nh4_s_yr(gy_p_9,:,sig),1,'omitnan')./(mean(TN_s_yr(gy_p_9,:,sig),1,'omitnan') - mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'));
chl_gy=mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan');
corrcoef(chl_gy(~isnan(NH4_ov_DON)),NH4_ov_DON(~isnan(NH4_ov_DON)))



gy_p_9 = [171 172 173 174 175 381 183 184 185 ];
sig=3;
fontSizeTick = 10
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(chl_s_yr))
xticklabels(1989:2021)
xtickangle(90)
grid on;% ylim([0 50]); 
xlim([1 length(sj.yr_tn_load)]);
ylabel('(TN-DIN)/NH4 (mmol/L)');
plot(9:33,(mean(TN_s_yr(gy_p_9,:,sig),1,'omitnan') - mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'))./mean(nh4_s_yr(gy_p_9,:,sig),1,'omitnan'),'k');
yyaxis right
plot(9:33,mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylabel('Chl.a(ug/L)');
ylim([0 15]);
% plot(mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'),'r');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('TN_load_vs_bay_chl.png'),'png');
xlim([9 length(sj.yr_tp_load)]);


sig=3
NH4_ov_DON=mean(nh4_s_yr(gy_p_9,:,sig),1,'omitnan')./(mean(TN_s_yr(gy_p_9,:,sig),1,'omitnan') - mean(DIN_s_yr(gy_p_9,:,sig),1,'omitnan'));
chl_gy=mean(chl_s_yr(gy_p_9,:,sig),1,'omitnan');
corrcoef(chl_gy(~isnan(NH4_ov_DON)),NH4_ov_DON(~isnan(NH4_ov_DON)))



%% monthly
sig=6;
fontSizeTick = 10
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_tn_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; %ylim([0 1300]); 
xlim([1 length(sj.yr_tn_load)]);
ylabel('(TN-DIN)/NH4 (mmol/L)');
plot(9:33,(mean(TN_s(gy_p_9,i:12:end,sig),1,'omitnan') - mean(DIN_s(gy_p_9,i:12:end,sig),1,'omitnan'))./mean(nh4_s(gy_p_9,i:12:end,sig),1,'omitnan'),'k');
yyaxis right
plot(9:33,mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylim([0 15]);
ylabel('Chl.a(ug/L)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('TN_load_vs_bay_chl_',num2str(i),'mon.png'),'png');
xlim([9 inf]);

NH4_ov_DON=(mean(TN_s(gy_p_9,i:12:end,sig),1,'omitnan') - mean(DIN_s(gy_p_9,i:12:end,sig),1,'omitnan'))./mean(nh4_s(gy_p_9,i:12:end,sig),1,'omitnan');
chl_gy=mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan');
corrcoef(chl_gy(~isnan(NH4_ov_DON)),NH4_ov_DON(~isnan(NH4_ov_DON)))

% saveas(gcf,strcat('TN_load_vs_bay_chl_fr1997_',num2str(i),'mon.png'),'png');
end




%% monthly
sig=6;
fontSizeTick = 10
figPos = [0 0 5 4];
for i = [2,5,8,11]
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
xticks(1:length(sj.yr_tn_load))
xticklabels(1989:2021)
xtickangle(90)
grid on; %ylim([0 1300]); 
xlim([1 length(sj.yr_tn_load)]);
ylabel('(TN-DIN)/NH4 (mmol/L)');
plot(9:33,(mean(TN_s(gy_p_9,i:12:end,sig),1,'omitnan') - mean(DIN_s(gy_p_9,i:12:end,sig),1,'omitnan'))./mean(nh4_s(gy_p_9,i:12:end,sig),1,'omitnan'),'k');
yyaxis right
plot(9:33,mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan'),'r');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
ylim([0 15]);
ylabel('Chl.a(ug/L)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
% saveas(gcf,strcat('TN_load_vs_bay_chl_',num2str(i),'mon.png'),'png');
xlim([9 inf]);

NH4_ov_DON=(mean(TN_s(gy_p_9,i:12:end,sig),1,'omitnan') - mean(DIN_s(gy_p_9,i:12:end,sig),1,'omitnan'))./mean(nh4_s(gy_p_9,i:12:end,sig),1,'omitnan');
chl_gy=mean(chl_s(gy_p_9,i:12:end,sig),1,'omitnan');
corrcoef(chl_gy(~isnan(NH4_ov_DON)),NH4_ov_DON(~isnan(NH4_ov_DON)))

% saveas(gcf,strcat('TN_load_vs_bay_chl_fr1997_',num2str(i),'mon.png'),'png');
end


plot(yeosu_asos_temp(2:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])


plot(yeosu_asos_temp(656+1:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])

plot(yeosu_asos_pair(656+1:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])

plot(yeosu_asos_srad(656+1:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])

plot(yeosu_asos_speed(656+1:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])

%% evap !!
plot(yeosu_asos_evap_s(656+0:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])

%% soil !!
plot(yeosu_asos_soil_t(656+1:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])



plot(yeosu_aws_mintemp(2:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])

plot(raw_aws_yeosu(2:12:end,7),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])



plot(yeosu_aws_speed(2:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])

plot(yeosu_aws_maxspeed(2:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])

plot(yeosu_aws_maxspeed_dir(2:12:end),'r-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])


plot(yeosu_aws_rain(2:12:end),'r-'); hold on;
plot(yeosu_aws_rain(1:12:end),'b-')
xticks(1:length(1997:2021))
xticklabels(1997:2021)
xtickangle(45)
grid on;
xlim([1 length(1997:2021)])




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