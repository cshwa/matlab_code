close all; clear; clc;

set_path = 'D:\Àå±â»ýÅÂ\Dynamic\KOEM\koem_yoonja_kang\ÃøÁ¤¸Á-¸ÞÀÌ½º(¹®È­)\';

sheet_tail = [];

[raw2 txt2]=xlsread([set_path,'1997_ppb.xlsx'],['2¿ù_ppb']);

[raw0 txt0]=xlsread([set_path,'2021_ppb.xlsx'],['2¿ù',sheet_tail]);

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
geoshow(A,R); ylim([33 39.5]);
title({'KOEM station'})
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');  
for i = 1:465
text(lon_num(i),lat_num(i), num2str(i),'fontsize',8,'fontweight','bold','color','w') ;
end
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

% ±¤¾çÇ×, ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5, ¿©¼ö2, ¿©¼ö3, ¿©¼ö1

% 2011~ ÀÌÈÄ lon, lat ¹× ÇØ´ç ÁöÁ¡ ¼ö½É ±âÀÔµÇ±â ½ÃÀÛ

%% 1997~2012
% ±¤¾çÇ× = 314
% ±¤¾ç 1:5 = 136:140
% ¿©¼ö 1:3 = 145:147
%% 2013~2019
%2013~
% 2¿ù (not 2¿ù_ppb)

% 2013
% ±¤¾çÇ× = 354
% ±¤¾ç 1:5 =159:163
% ¿©¼ö 1:3 = 171:173

% 2014~2021
% ±¤¾çÇ× = 384
% ±¤¾ç 1:5 = 174:178
% ¿©¼ö 1:3 = 186:188

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
        col_ph = [17, 18]
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
        col_ph = [17, 18]
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
        
    [raw2 txt2]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['2¿ù',sheet_tail]);
    [raw5 txt5]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['5¿ù',sheet_tail]);
    [raw8 txt8]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['8¿ù',sheet_tail]);
    [raw11 txt11]=xlsread([set_path,num2str(i),'_ppb.xlsx'],['11¿ù',sheet_tail]);
    
    
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
    
    %% make row_pick
    for k = 1:length(st_ref)

    tempo_check=strcmp(st_ref{k},st_pick);
    
    check_idx = find(tempo_check == 1);
    
        if length(check_idx) == 1
            row_pick(k) = check_idx;
        elseif length(check_idx) == 0
            row_pick(k) = NaN;
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
       for l = 1:length(row_pick)
      
    if isnan(row_pick(l)) == 0
        if i > 2010
        mon_d(l,obs_tdx(2))=raw5(row_pick(l),col_time(1)); % day on month
        time(l,obs_tdx(2))=raw5(row_pick(l),col_time(2))*24; % time on day
        end  
%         temp_sur(i,obs_tdx(2))=raw5(row_pick(l),col_temp(1));
        temp_sur(l,obs_tdx(2))=raw5(row_pick(l),col_temp(1));
        temp_bot(l,obs_tdx(2))=raw5(row_pick(l),col_temp(2));
        salt_sur(l,obs_tdx(2))=raw5(row_pick(l),col_salt(1));
        salt_bot(l,obs_tdx(2))=raw5(row_pick(l),col_salt(2));
        no3_sur(l,obs_tdx(2))=raw5(row_pick(l),col_no3(1));
        no3_bot(l,obs_tdx(2))=raw5(row_pick(l),col_no3(2));
        nh4_sur(l,obs_tdx(2))=raw5(row_pick(l),col_nh4(1));
        nh4_bot(l,obs_tdx(2))=raw5(row_pick(l),col_nh4(2));
        chl_sur(l,obs_tdx(2))=raw5(row_pick(l),col_chl(1));
        chl_bot(l,obs_tdx(2))=raw5(row_pick(l),col_chl(2));
        po4_sur(l,obs_tdx(2))=raw5(row_pick(l),col_po4(1));
        po4_bot(l,obs_tdx(2))=raw5(row_pick(l),col_po4(2));
        do_sur(l,obs_tdx(2))=raw5(row_pick(l),col_do(1));
        do_bot(l,obs_tdx(2))=raw5(row_pick(l),col_do(2));
        DIN_sur(l,obs_tdx(2))=raw5(row_pick(l),col_DIN(1));
        DIN_bot(l,obs_tdx(2))=raw5(row_pick(l),col_DIN(2));
        ss_sur(l,obs_tdx(2))=raw5(row_pick(l),col_ss(1));
        ss_bot(l,obs_tdx(2))=raw5(row_pick(l),col_ss(2));   
        si_sur(l,obs_tdx(2))=raw5(row_pick(l),col_si(1));
        si_bot(l,obs_tdx(2))=raw5(row_pick(l),col_si(2));      
        ph_sur(l,obs_tdx(2))=raw5(row_pick(l),col_ph(1));
        ph_bot(l,obs_tdx(2))=raw5(row_pick(l),col_ph(2));        
        cod_sur(l,obs_tdx(2))=raw5(row_pick(l),col_cod(1));
        cod_bot(l,obs_tdx(2))=raw5(row_pick(l),col_cod(2));        
        no2_sur(l,obs_tdx(2))=raw5(row_pick(l),col_no2(1));
        no2_bot(l,obs_tdx(2))=raw5(row_pick(l),col_no2(2));
        TN_sur(l,obs_tdx(2))=raw5(row_pick(l),col_TN(1));
        TN_bot(l,obs_tdx(2))=raw5(row_pick(l),col_TN(2));       
        TP_sur(l,obs_tdx(2))=raw5(row_pick(l),col_TP(1));
        TP_bot(l,obs_tdx(2))=raw5(row_pick(l),col_TP(2));      
        DIP_sur(l,obs_tdx(2))=raw5(row_pick(l),col_dip(1));
        DIP_bot(l,obs_tdx(2))=raw5(row_pick(l),col_dip(2));     
        secchi_sur(l,obs_tdx(2))=raw5(row_pick(l),col_secchi(1));
    else isnan(row_pick(l)) == 1
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

       for l = 1:length(row_pick)
     
    if isnan(row_pick(l)) == 0
         if i > 2010
        mon_d(l,obs_tdx(3))=raw8(row_pick(l),col_time(1)); % day on month
        time(l,obs_tdx(3))=raw8(row_pick(l),col_time(2))*24; % time on day
        end  
%         temp_sur(i,obs_tdx(3))=raw8(row_pick(l),col_temp(1));
        temp_sur(l,obs_tdx(3))=raw8(row_pick(l),col_temp(1));
        temp_bot(l,obs_tdx(3))=raw8(row_pick(l),col_temp(2));
        salt_sur(l,obs_tdx(3))=raw8(row_pick(l),col_salt(1));
        salt_bot(l,obs_tdx(3))=raw8(row_pick(l),col_salt(2));
        no3_sur(l,obs_tdx(3))=raw8(row_pick(l),col_no3(1));
        no3_bot(l,obs_tdx(3))=raw8(row_pick(l),col_no3(2));
        nh4_sur(l,obs_tdx(3))=raw8(row_pick(l),col_nh4(1));
        nh4_bot(l,obs_tdx(3))=raw8(row_pick(l),col_nh4(2));
        chl_sur(l,obs_tdx(3))=raw8(row_pick(l),col_chl(1));
        chl_bot(l,obs_tdx(3))=raw8(row_pick(l),col_chl(2));
        po4_sur(l,obs_tdx(3))=raw8(row_pick(l),col_po4(1));
        po4_bot(l,obs_tdx(3))=raw8(row_pick(l),col_po4(2));
        do_sur(l,obs_tdx(3))=raw8(row_pick(l),col_do(1));
        do_bot(l,obs_tdx(3))=raw8(row_pick(l),col_do(2));
        DIN_sur(l,obs_tdx(3))=raw8(row_pick(l),col_DIN(1));
        DIN_bot(l,obs_tdx(3))=raw8(row_pick(l),col_DIN(2));
        ss_sur(l,obs_tdx(3))=raw8(row_pick(l),col_ss(1));
        ss_bot(l,obs_tdx(3))=raw8(row_pick(l),col_ss(2));   
        si_sur(l,obs_tdx(3))=raw8(row_pick(l),col_si(1));
        si_bot(l,obs_tdx(3))=raw8(row_pick(l),col_si(2));      
        ph_sur(l,obs_tdx(3))=raw8(row_pick(l),col_ph(1));
        ph_bot(l,obs_tdx(3))=raw8(row_pick(l),col_ph(2));        
        cod_sur(l,obs_tdx(3))=raw8(row_pick(l),col_cod(1));
        cod_bot(l,obs_tdx(3))=raw8(row_pick(l),col_cod(2));        
        no2_sur(l,obs_tdx(3))=raw8(row_pick(l),col_no2(1));
        no2_bot(l,obs_tdx(3))=raw8(row_pick(l),col_no2(2));
        TN_sur(l,obs_tdx(3))=raw8(row_pick(l),col_TN(1));
        TN_bot(l,obs_tdx(3))=raw8(row_pick(l),col_TN(2));       
        TP_sur(l,obs_tdx(3))=raw8(row_pick(l),col_TP(1));
        TP_bot(l,obs_tdx(3))=raw8(row_pick(l),col_TP(2));      
        DIP_sur(l,obs_tdx(3))=raw8(row_pick(l),col_dip(1));
        DIP_bot(l,obs_tdx(3))=raw8(row_pick(l),col_dip(2));     
        secchi_sur(l,obs_tdx(3))=raw8(row_pick(l),col_secchi(1));
    else isnan(row_pick(l)) == 1
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
       for l = 1:length(row_pick)
      
    if isnan(row_pick(l)) == 0
        if i > 2010
        mon_d(l,obs_tdx(4))=raw11(row_pick(l),col_time(1)); % day on month
        time(l,obs_tdx(4))=raw11(row_pick(l),col_time(2))*24; % time on day
        end  
%         temp_sur(i,obs_tdx(4))=raw11(row_pick(l),col_temp(1));
        temp_sur(l,obs_tdx(4))=raw11(row_pick(l),col_temp(1));
        temp_bot(l,obs_tdx(4))=raw11(row_pick(l),col_temp(2));
        salt_sur(l,obs_tdx(4))=raw11(row_pick(l),col_salt(1));
        salt_bot(l,obs_tdx(4))=raw11(row_pick(l),col_salt(2));
        no3_sur(l,obs_tdx(4))=raw11(row_pick(l),col_no3(1));
        no3_bot(l,obs_tdx(4))=raw11(row_pick(l),col_no3(2));
        nh4_sur(l,obs_tdx(4))=raw11(row_pick(l),col_nh4(1));
        nh4_bot(l,obs_tdx(4))=raw11(row_pick(l),col_nh4(2));
        chl_sur(l,obs_tdx(4))=raw11(row_pick(l),col_chl(1));
        chl_bot(l,obs_tdx(4))=raw11(row_pick(l),col_chl(2));
        po4_sur(l,obs_tdx(4))=raw11(row_pick(l),col_po4(1));
        po4_bot(l,obs_tdx(4))=raw11(row_pick(l),col_po4(2));
        do_sur(l,obs_tdx(4))=raw11(row_pick(l),col_do(1));
        do_bot(l,obs_tdx(4))=raw11(row_pick(l),col_do(2));
        DIN_sur(l,obs_tdx(4))=raw11(row_pick(l),col_DIN(1));
        DIN_bot(l,obs_tdx(4))=raw11(row_pick(l),col_DIN(2));
        ss_sur(l,obs_tdx(4))=raw11(row_pick(l),col_ss(1));
        ss_bot(l,obs_tdx(4))=raw11(row_pick(l),col_ss(2));   
        si_sur(l,obs_tdx(4))=raw11(row_pick(l),col_si(1));
        si_bot(l,obs_tdx(4))=raw11(row_pick(l),col_si(2));      
        ph_sur(l,obs_tdx(4))=raw11(row_pick(l),col_ph(1));
        ph_bot(l,obs_tdx(4))=raw11(row_pick(l),col_ph(2));        
        cod_sur(l,obs_tdx(4))=raw11(row_pick(l),col_cod(1));
        cod_bot(l,obs_tdx(4))=raw11(row_pick(l),col_cod(2));        
        no2_sur(l,obs_tdx(4))=raw11(row_pick(l),col_no2(1));
        no2_bot(l,obs_tdx(4))=raw11(row_pick(l),col_no2(2));
        TN_sur(l,obs_tdx(4))=raw11(row_pick(l),col_TN(1));
        TN_bot(l,obs_tdx(4))=raw11(row_pick(l),col_TN(2));       
        TP_sur(l,obs_tdx(4))=raw11(row_pick(l),col_TP(1));
        TP_bot(l,obs_tdx(4))=raw11(row_pick(l),col_TP(2));      
        DIP_sur(l,obs_tdx(4))=raw11(row_pick(l),col_dip(1));
        DIP_bot(l,obs_tdx(4))=raw11(row_pick(l),col_dip(2));     
        secchi_sur(l,obs_tdx(4))=raw11(row_pick(l),col_secchi(1));
    else isnan(row_pick(l)) == 1
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

save('yoonjakangs_koem_data_monthly_v4_all_points_2021.mat');

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
% legend('±¤¾ç1', '±¤¾ç2', '±¤¾ç3', '±¤¾ç4', '±¤¾ç5');
xtickangle(45);



%     [raw20 txt20]=xlsread([set_path,num2str(2010),'_ppb.xlsx'],'2¿ù_ppb');
%     [raw21 txt21]=xlsread([set_path,num2str(2011),'_ppb.xlsx'],'2¿ù_ppb');
%     [raw1 txt1]=xlsread([set_path,num2str(2019),'_ppb.xlsx'],'2¿ù');