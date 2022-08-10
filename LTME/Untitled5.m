% filename = 'KO0002.txt';
% 
% data = textread('KO0002.E', '', 'delimiter', ',','emptyvalue', NaN);

close all; clear;clc;

% ref date (yy-mm)
ik=0;
for iy = 1964:2010
    for im = 1:12
        ik= ik +1;
        iy_s = num2str(iy);
        im_s = num2str(im,'%02d');
        ref_time_yymm{ik} = [ iy_s(3:4),'-',im_s ]
    end
end

% reference location

ref_lon = 115.00:2:164.00;
ref_lat = 15.00:2:52.00;

% [ref_lon_m ref_lat_m] = meshgrid(ref_lon_pre,ref_lat_pre);
[ref_lat_m ref_lon_m] = meshgrid(ref_lat,ref_lon);


%standard depth 
std_depth = [[1:10:100]'; [200:100:5500]';];

%std form of each variables
       temp_3d=NaN(size(ref_lon_m,1),size(ref_lon_m,2),length(std_depth),length(ref_time_yymm));
       salt_3d=NaN(size(ref_lon_m,1),size(ref_lon_m,2),length(std_depth),length(ref_time_yymm));
       tp_3d=NaN(size(ref_lon_m,1),size(ref_lon_m,2),length(std_depth),length(ref_time_yymm));
       no3_3d=NaN(size(ref_lon_m,1),size(ref_lon_m,2),length(std_depth),length(ref_time_yymm));
       no2_3d=NaN(size(ref_lon_m,1),size(ref_lon_m,2),length(std_depth),length(ref_time_yymm));
       nh3_3d=NaN(size(ref_lon_m,1),size(ref_lon_m,2),length(std_depth),length(ref_time_yymm));
       ph_3d=NaN(size(ref_lon_m,1),size(ref_lon_m,2),length(std_depth),length(ref_time_yymm));
       chl_3d=NaN(size(ref_lon_m,1),size(ref_lon_m,2),length(std_depth),length(ref_time_yymm));
       pha_3d=NaN(size(ref_lon_m,1),size(ref_lon_m,2),length(std_depth),length(ref_time_yymm));


cd C:\Users\user\Desktop\jma
dir_list = dir('*.E'); % search dir.
for f_dir =  1:length(dir_list)
    clearvars -except dir_list *_3d std_depth ref_*_m ref_time_yymm f_dir
    
    cd(dir_list(f_dir).name)

    list_file = dir('*.E');
      
tz=0;       
for f = 1:length(list_file)
% for    f = length(list_file):length(list_file)
% f=1;

% filename = 'KO0002.E';

clearvars -except f dir_list f_dir list_file *_3d std_depth ref_*_m ref_time_yymm tz 

filename = list_file(f).name;

xx = filename(1:2);
yy = filename(3:4);
mm = filename(5:6);


%name : xxyymm
% xx: kruise title
% yy: year
% mm: month

file_cont = textread(filename, '%s', 'delimiter', '\n', ...
                'whitespace', '');
            
            
            
% readmatrix('KO0002.txt')
% 
% readtable(filename)

% j=0;
% for i = 1:length(file_cont)
%     if length(strfind(file_cont{i},'AO-')) ~= 0
%         j=j+1;
%         start_head2(j) =i;
%     end
% end

    j=0
    for i = 1:length(file_cont)
        if length(strfind(file_cont{i},[yy mm '='])) ~= 0
            j=j+1;
            start_head2(j) =i;
        end
    end


% % head
for i=1:length(start_head2)
% STATION NO 1
    st_num{i} = file_cont{start_head2(i)}(1:7);
% LATITUDE 9     deg,min
    lat_raw{i} = file_cont{start_head2(i)}(9:13);
% LONGITUDE 17 
    lon_raw{i} = file_cont{start_head2(i)}(17:22);
% DATE/TIME 26   mon,day,hr ( start / end  )
    time_raw_start{i} = file_cont{start_head2(i)}(26:35);
    time_raw_end{i} = file_cont{start_head2(i)}(36:46);
% SUB STN NO 116     ÜÍð¾ö´ïÃÛã?  -> sub- station num
    sub_st_num{i} = file_cont{start_head2(i)}(116:120);
end

% % no need head anymore
for i=1:length(start_head2)
% STATION NO 1
file_cont{start_head2(i)}=[];
end

clearvars depth_raw
k=0;
% % body
for i=1:length(start_head2) %st_num
    if i ~= length(start_head2)
       for j = start_head2(i)+2:start_head2(i+1)-1 %dep for each st.
           k=k+1;
    %     start_head2(i)+1+j:start_head2(i+1)-1 
    
    % depth for observ. 
    if length(str2num(file_cont{j}(17:20))) ~= 0
        depth_raw{i}(k) = str2num(file_cont{j}(17:20));
    end
        
    if length(str2num(file_cont{j}(17:20))) ~= 0  %if there is no depth then no data
    % TEMP(OBS) 22 F6.3 CTD â©? (1990 Ò´?ð·?ÓøÙÍàü (ITS-90))
        if length(str2num(file_cont{j}(22:26))) == 0
            temp_raw{i}(k) = NaN;
        else
            temp_raw{i}(k) = str2num(file_cont{j}(22:26));
        end
        % SAL(OBS) 28 F6.3 CTD ?ÝÂ (1978 Ò´?éÄ?ÝÂ (PSS-78))
        if length(str2num(file_cont{j}(28:34))) == 0
            salt_raw{i}(k) = NaN;
        else
            salt_raw{i}(k) = str2num(file_cont{j}(28:34));
        end
        % DO 35 I3 éÁðíß«áÈÒØÓø£¨?mol/l)
        if length(str2num(file_cont{j}(35:38))) == 0
            do_raw{i}(k) = NaN;
        else
            do_raw{i}(k) = str2num(file_cont{j}(35:38));
        end
        % PO4-P 39 F4.2 «ê«óß«?£¨?mol/l)
        if length(str2num(file_cont{j}(39:42))) == 0
            po4_raw{i}(k) = NaN;
        else    
            po4_raw{i}(k) = str2num(file_cont{j}(39:42));
        end

        % T-P 44 F4.2 îï«ê«óß«? (?mol/l)
        if length(str2num(file_cont{j}(44:47))) == 0
            tp_raw{i}(k) = NaN;
        else    
            tp_raw{i}(k) = str2num(file_cont{j}(44:47));
        end
        % NO3-N 49 F4.1 õ¦ß«?£¨?mol/l£©
        if length(str2num(file_cont{j}(49:52))) == 0
            no3_raw{i}(k) = NaN;
        else
            no3_raw{i}(k) = str2num(file_cont{j}(49:52));
        end
        % NO2-N 54 F4.2 ?õ¦ß«?£¨?mol/l£©
        if length(str2num(file_cont{j}(54:57))) == 0
            no2_raw{i}(k) = NaN;
        else
            no2_raw{i}(k) = str2num(file_cont{j}(54:57));
        end
        % NH3-N 59 F4.2 «¢«ó«â«Ë«¢ (?mol/l)
        if length(str2num(file_cont{j}(59:62))) == 0
            nh3_raw{i}(k) = NaN;
        else
            nh3_raw{i}(k) = str2num(file_cont{j}(59:62));
        end
        % PH 64 F4.2 25 ¡ÉªËªªª±ªëâ©áÈ«¤«ª«óÒØÓøò¦?£¨NBS «¹«±?«ë£©
        if length(str2num(file_cont{j}(64:67))) == 0
            ph_raw{i}(k) = NaN;
        else
            ph_raw{i}(k) = str2num(file_cont{j}(64:67));
        end
        % CHL 69 F6.2 «¯«í«í«Õ«£«ë a (?g/l)
        if length(str2num(file_cont{j}(69:74))) == 0
            chl_raw{i}(k) = NaN;
        else
            chl_raw{i}(k) = str2num(file_cont{j}(69:74));
        end
        % PHA 76 F6.2 «Õ«£«ª«Õ«£«Á«ó (?g/l)  
        if length(str2num(file_cont{j}(76:81))) == 0
            pha_raw{i}(k) = NaN;
        else
            pha_raw{i}(k) = str2num(file_cont{j}(76:81));
        end
    end
       end
       
        %     start_head2(i+1)-2 % finished line
        elseif i == length(start_head2)
           for j = start_head2(i)+2:length(file_cont) %dep for each st.
               k=k+1;
        %     start_head2(i)+1+j:start_head2(i+1)-1
        % STATION NO 1
        
    if length(str2num(file_cont{j}(17:20))) ~= 0
        depth_raw{i}(k) = str2num(file_cont{j}(17:20));
%     elseif j==start_head2(i)+2 & length(str2num(file_cont{j}(17:20))) == 0
%         depth_raw{i} = str2num(file_cont{j}(17:20));
    end
            
        if length(str2num(file_cont{j}(17:20))) ~= 0 
    % TEMP(OBS) 22 F6.3 CTD â©? (1990 Ò´?ð·?ÓøÙÍàü (ITS-90))
            if length(str2num(file_cont{j}(22:26))) == 0
                temp_raw{i}(k) = NaN;
            else
                temp_raw{i}(k) = str2num(file_cont{j}(22:26));
            end
            % SAL(OBS) 28 F6.3 CTD ?ÝÂ (1978 Ò´?éÄ?ÝÂ (PSS-78))
            if length(str2num(file_cont{j}(28:34))) == 0
                salt_raw{i}(k) = NaN;
            else
                salt_raw{i}(k) = str2num(file_cont{j}(28:34));
            end
            % DO 35 I3 éÁðíß«áÈÒØÓø£¨?mol/l)
            if length(str2num(file_cont{j}(35:38))) == 0
                do_raw{i}(k) = NaN;
            else
                do_raw{i}(k) = str2num(file_cont{j}(35:38));
            end
            % PO4-P 39 F4.2 «ê«óß«?£¨?mol/l)
            if length(str2num(file_cont{j}(39:42))) == 0
                po4_raw{i}(k) = NaN;
            else    
                po4_raw{i}(k) = str2num(file_cont{j}(39:42));
            end

            % T-P 44 F4.2 îï«ê«óß«? (?mol/l)
            if length(str2num(file_cont{j}(44:47))) == 0
                tp_raw{i}(k) = NaN;
            else    
                tp_raw{i}(k) = str2num(file_cont{j}(44:47));
            end
            % NO3-N 49 F4.1 õ¦ß«?£¨?mol/l£©
            if length(str2num(file_cont{j}(49:52))) == 0
                no3_raw{i}(k) = NaN;
            else
                no3_raw{i}(k) = str2num(file_cont{j}(49:52));
            end
            % NO2-N 54 F4.2 ?õ¦ß«?£¨?mol/l£©
            if length(str2num(file_cont{j}(54:57))) == 0
                no2_raw{i}(k) = NaN;
            else
                no2_raw{i}(k) = str2num(file_cont{j}(54:57));
            end
            % NH3-N 59 F4.2 «¢«ó«â«Ë«¢ (?mol/l)
            if length(str2num(file_cont{j}(59:62))) == 0
                nh3_raw{i}(k) = NaN;
            else
                nh3_raw{i}(k) = str2num(file_cont{j}(59:62));
            end
            % PH 64 F4.2 25 ¡ÉªËªªª±ªëâ©áÈ«¤«ª«óÒØÓøò¦?£¨NBS «¹«±?«ë£©
            if length(str2num(file_cont{j}(64:67))) == 0
                ph_raw{i}(k) = NaN;
            else
                ph_raw{i}(k) = str2num(file_cont{j}(64:67));
            end
            % CHL 69 F6.2 «¯«í«í«Õ«£«ë a (?g/l)
            if length(str2num(file_cont{j}(69:74))) == 0
                chl_raw{i}(k) = NaN;
            else
                chl_raw{i}(k) = str2num(file_cont{j}(69:74));
            end
            % PHA 76 F6.2 «Õ«£«ª«Õ«£«Á«ó (?g/l)  
            if length(str2num(file_cont{j}(76:81))) == 0
                pha_raw{i}(k) = NaN;
            else
                pha_raw{i}(k) = str2num(file_cont{j}(76:81));
            end
       end  
           end
    end      
   k=0; %initialize
end

% obs time
time_comp = {[yy,'-',mm]};



% time compare

% make it to cell-array
% for i = 1:length(yymmdd_txt)
%     yymmdd_txt_c{i,1} = yymmdd_txt(i,:); % raw
% end

% pick matched date from water temp date
clearvars indx_raw
indx_raw = find([strcmp(ref_time_yymm, time_comp)] == 1);


% obs location
for i=1:length(lon_raw)
lon_raw_deg(i) = str2num(lon_raw{i}(1:3));
lat_raw_deg(i) = str2num(lat_raw{i}(1:2));
lon_raw_min(i) = str2num(lon_raw{i}(5:6));
lat_raw_min(i) = str2num(lat_raw{i}(4:5));
end

lon_raw_decimal= lon_raw_deg + (lon_raw_min ./ 60);
lat_raw_decimal= lat_raw_deg + (lat_raw_min ./ 60);


% ref_lon_pre = 164.00:-2:115.00;
% ref_lat_pre = 52.00:-2:15.00;
% 
% [ref_lat_m ref_lon_m] = meshgrid(ref_lat_pre,ref_lon_pre);
% 
% clearvars ref_lon_m ref_lat_m
% ref_lon = 164.00:-2:115.00;
% ref_lat = 52.00:-2:15.00;
% 
% [ref_lon_m ref_lat_m] = meshgrid(ref_lon_pre,ref_lat_pre);
% 
% clearvars ref_lon_m ref_lat_m
% ref_lon = 115.00:2:164.00;
% ref_lat = 52.00:-2:15.00;
% 
% [ref_lon_m ref_lat_m] = meshgrid(ref_lon_pre,ref_lat_pre);


% location compare

for i = 1:length(lon_raw_decimal)
clearvars temp_2d
    temp_2d = sqrt((ref_lon_m - lon_raw_decimal(i)).^2 + (ref_lat_m - lat_raw_decimal(i)).^2);
    temp_2d = temp_2d %.* mask;
    temp_near = find(nanmin(nanmin(temp_2d))==temp_2d);
    near_point(i) = temp_near(1);
    dist_2d(:,:,i) = temp_2d;
end

for i = 1:length(near_point)
    clearvars row col
    [row,col]=ind2sub(size(ref_lon_m),near_point(i)); %% 1d to 2d index
    near_point_2d(i,:) = [row, col]; 
end

% ref_lon_m(near_point_2d(1,1),near_point_2d(1,2))
% ref_lat_m(near_point_2d(1,1),near_point_2d(1,2))
% 
% ref_lon_m(near_point_2d(end,1),near_point_2d(end,2))
% ref_lat_m(near_point_2d(end,1),near_point_2d(end,2))


% make it standard depth form

% make form of grid
clearvars *_num
% for i=1:length(start_head2) %st_num
%     % depth for observ.
%         depth_raw_num{i}=str2num(depth_raw{i});
%     % TEMP(OBS) 22 F6.3 CTD â©? (1990 Ò´?ð·?ÓøÙÍàü (ITS-90))
%         temp_raw_num{i}=str2num(temp_raw{i});
%     % SAL(OBS) 28 F6.3 CTD ?ÝÂ (1978 Ò´?éÄ?ÝÂ (PSS-78))
%         salt_raw_num{i}=str2num(salt_raw{i});
%     % DO 35 I3 éÁðíß«áÈÒØÓø£¨?mol/l)
%         do_raw_num{i}=str2num(do_raw{i});
%     % PO4-P 39 F4.2 «ê«óß«?£¨?mol/l)
%         po4_raw_num{i}=str2num(po4_raw{i});
%     % T-P 44 F4.2 îï«ê«óß«? (?mol/l)
%         tp_raw_num{i}=str2num(tp_raw{i});
%     % NO3-N 49 F4.1 õ¦ß«?£¨?mol/l£©
%         no3_raw_num{i}=str2num(no3_raw{i});
%     % NO2-N 54 F4.2 ?õ¦ß«?£¨?mol/l£©
%         no2_raw_num{i}=str2num(no2_raw{i});
%     % NH3-N 59 F4.2 «¢«ó«â«Ë«¢ (?mol/l)
%         nh3_raw_num{i}=str2num(nh3_raw{i});
%     % PH 64 F4.2 25 ¡ÉªËªªª±ªëâ©áÈ«¤«ª«óÒØÓøò¦?£¨NBS «¹«±?«ë£©
%         ph_raw_num{i}=str2num(ph_raw{i});
%     % CHL 69 F6.2 «¯«í«í«Õ«£«ë a (?g/l)
%         chl_raw_num{i}=str2num(chl_raw{i});
%     % PHA 76 F6.2 «Õ«£«ª«Õ«£«Á«ó (?g/l)  
%         pha_raw_num{i}=str2num(pha_raw{i});
% end


       
logical_depth = exist('depth_raw');

if logical_depth == 1  %logical switch for no depth then no analyze
for i=1:length(depth_raw) %st_num
    clearvars tempo_*_raw *_1d
    tempo_dep_raw =depth_raw{i}; %each st. depth
    tempo_temp_raw =temp_raw{i}; %each st. 
    tempo_salt_raw =salt_raw{i}; %each st.
    tempo_do_raw =do_raw{i}; %each st. 
    tempo_po4_raw =po4_raw{i}; %each st. 
    tempo_tp_raw =tp_raw{i}; %each st. 
    tempo_no3_raw =no3_raw{i}; %each st. 
    tempo_no2_raw =no2_raw{i}; %each st. 
    tempo_nh3_raw =nh3_raw{i}; %each st. 
    tempo_ph_raw =ph_raw{i}; %each st. 
    tempo_chl_raw =chl_raw{i}; %each st. 
    tempo_pha_raw =pha_raw{i}; %each st. 

   if length(min(tempo_dep_raw):max(tempo_dep_raw)) < 10 %filter out abnormal data for obs depth
    tempo_dep_raw(:) = NaN;
    tempo_temp_raw(:) = NaN;
    tempo_salt_raw(:) = NaN;
    tempo_do_raw(:) = NaN;
    tempo_po4_raw(:) = NaN;
    tempo_tp_raw(:) = NaN;
    tempo_no3_raw(:) = NaN;
    tempo_no2_raw(:) = NaN;
    tempo_ph_raw(:) = NaN;
    tempo_chl_raw(:) = NaN;
    tempo_pha_raw(:) = NaN;
   end

   
   [dep_unique_raw, ia,ic] = unique(tempo_dep_raw); % remaining first redundant depth   
   if length(dep_unique_raw) ~= length(tempo_dep_raw)
   clearvars tempo_dep_raw
   tempo_dep_raw = dep_unique_raw;
   tempo_temp_raw = tempo_temp_raw(ia);
    tempo_salt_raw = tempo_salt_raw(ia);
    tempo_do_raw = tempo_do_raw(ia);
    tempo_po4_raw = tempo_po4_raw(ia);
    tempo_tp_raw = tempo_tp_raw(ia);
    tempo_nh3_raw = tempo_nh3_raw(ia);
    tempo_no3_raw = tempo_no3_raw(ia);
    tempo_no2_raw = tempo_no2_raw(ia);
    tempo_ph_raw = tempo_ph_raw(ia);
    tempo_chl_raw = tempo_chl_raw(ia);
    tempo_pha_raw = tempo_pha_raw(ia);
   end
   
       
       % ignore last missing value
       clearvars nan_* *_1d_in
       
       nan_temp = isnan(tempo_temp_raw);
       nan_salt = isnan(tempo_salt_raw);
       nan_tp = isnan(tempo_tp_raw);
       nan_no3 = isnan(tempo_no3_raw);
       nan_no2 = isnan(tempo_no2_raw);
       nan_nh3 = isnan(tempo_nh3_raw);
       nan_ph = isnan(tempo_ph_raw);
       nan_chl = isnan(tempo_chl_raw);
       nan_pha = isnan(tempo_pha_raw);
    
       
       temp_1d_in = tempo_temp_raw;
       salt_1d_in = tempo_salt_raw;
       tp_1d_in = tempo_tp_raw;
       no3_1d_in = tempo_no3_raw;
       no2_1d_in = tempo_no2_raw;
       nh3_1d_in = tempo_nh3_raw;
       ph_1d_in = tempo_ph_raw;
       chl_1d_in = tempo_chl_raw;
       pha_1d_in = tempo_pha_raw;
       
       
       if sum(isnan(tempo_temp_raw)) ~= length(tempo_temp_raw) & length(tempo_dep_raw(~nan_temp)) >= 3
           if sum(nan_temp) ~= 0
              temp_1d_in(nan_temp)=interp1(tempo_dep_raw(~nan_temp),tempo_temp_raw(~nan_temp),tempo_dep_raw(nan_temp));
           elseif sum(nan_temp) == 0
               temp_1d_in = tempo_temp_raw; 
           end
       elseif sum(isnan(tempo_temp_raw)) == length(tempo_temp_raw) | length(tempo_dep_raw(~nan_salt)) < 3
           temp_1d_in = NaN(1,length(std_depth));
       end
       
       
       if sum(isnan(tempo_salt_raw)) ~= length(tempo_salt_raw) & length(tempo_dep_raw(~nan_salt)) >= 3
           if sum(nan_salt) ~= 0
              salt_1d_in(nan_salt)=interp1(tempo_dep_raw(~nan_salt),tempo_salt_raw(~nan_salt),tempo_dep_raw(nan_salt));
           elseif sum(nan_salt) == 0
               salt_1d_in = tempo_salt_raw; 
           end
       elseif sum(isnan(tempo_salt_raw)) == length(tempo_salt_raw) | length(tempo_dep_raw(~nan_salt)) < 3
           salt_1d_in = NaN(1,length(std_depth));
       end
       
       
       if sum(isnan(tempo_tp_raw)) ~= length(tempo_tp_raw) & length(tempo_dep_raw(~nan_tp)) >= 3 
           if sum(nan_tp) ~= 0 
              tp_1d_in(nan_tp)=interp1(tempo_dep_raw(~nan_tp),tempo_tp_raw(~nan_tp),tempo_dep_raw(nan_tp));
           elseif sum(nan_tp) == 0 
               tp_1d_in = tempo_tp_raw; 
           end
       elseif sum(isnan(tempo_tp_raw)) == length(tempo_tp_raw) | length(tempo_dep_raw(~nan_tp)) < 3
           tp_1d_in = NaN(1,length(std_depth));
       end
       
       
       if sum(isnan(tempo_no3_raw)) ~= length(tempo_no3_raw) & length(tempo_dep_raw(~nan_no3)) >= 3
           if sum(nan_no3) ~= 0  
              no3_1d_in(nan_no3)=interp1(tempo_dep_raw(~nan_no3),tempo_no3_raw(~nan_no3),tempo_dep_raw(nan_no3));
           elseif sum(nan_no3) == 0 
               no3_1d_in = tempo_no3_raw; 
           end
       elseif sum(isnan(tempo_no3_raw)) == length(tempo_no3_raw) | length(tempo_dep_raw(~nan_no3)) < 3
           no3_1d_in = NaN(1,length(std_depth));
       end
       
       
       if sum(isnan(tempo_no2_raw)) ~= length(tempo_no2_raw) & length(tempo_dep_raw(~nan_no2)) >= 3
           if sum(nan_no2) ~= 0 
              no2_1d_in(nan_no2)=interp1(tempo_dep_raw(~nan_no2),tempo_no2_raw(~nan_no2),tempo_dep_raw(nan_no2));
           elseif sum(nan_no2) == 0 
               no2_1d_in = tempo_no2_raw; 
           end
       elseif sum(isnan(tempo_no2_raw)) == length(tempo_no2_raw) | length(tempo_dep_raw(~nan_no2)) < 3
           no2_1d_in = NaN(1,length(std_depth));
       end
       
       if sum(isnan(tempo_nh3_raw)) ~= length(tempo_nh3_raw) & length(tempo_dep_raw(~nan_nh3)) >= 3
           nan_nh3 = isnan(tempo_nh3_raw);
           if sum(nan_nh3) ~= 0 
              nh3_1d_in(nan_nh3)=interp1(tempo_dep_raw(~nan_nh3),tempo_nh3_raw(~nan_nh3),tempo_dep_raw(nan_nh3));
           elseif sum(nan_nh3) == 0 
               nh3_1d_in = tempo_nh3_raw; 
           end
       elseif sum(isnan(tempo_nh3_raw)) == length(tempo_nh3_raw) | length(tempo_dep_raw(~nan_nh3)) < 3
           nh3_1d_in = NaN(1,length(std_depth));
       end
       
       if sum(isnan(tempo_ph_raw)) ~= length(tempo_ph_raw) & length(tempo_dep_raw(~nan_ph)) >= 3
           nan_ph = isnan(tempo_ph_raw);
           if sum(nan_ph) ~= 0 
              ph_1d_in(nan_ph)=interp1(tempo_dep_raw(~nan_ph),tempo_ph_raw(~nan_ph),tempo_dep_raw(nan_ph));
           elseif sum(nan_ph) == 0 
               ph_1d_in = tempo_ph_raw; 
           end
       elseif sum(isnan(tempo_ph_raw)) == length(tempo_ph_raw) | length(tempo_dep_raw(~nan_ph)) < 3
           ph_1d_in = NaN(1,length(std_depth));
       end
       
       if sum(isnan(tempo_chl_raw)) ~= length(tempo_chl_raw) & length(tempo_dep_raw(~nan_chl)) >= 3
           nan_chl = isnan(tempo_chl_raw);
           if sum(nan_chl) ~= 0 
              chl_1d_in(nan_chl)=interp1(tempo_dep_raw(~nan_chl),tempo_chl_raw(~nan_chl),tempo_dep_raw(nan_chl));
           elseif sum(nan_chl) == 0 
               chl_1d_in = tempo_chl_raw; 
           end
       elseif sum(isnan(tempo_chl_raw)) == length(tempo_chl_raw) | length(tempo_dep_raw(~nan_chl)) < 3
           chl_1d_in = NaN(1,length(std_depth));
       end
       
       if sum(isnan(tempo_pha_raw)) ~= length(tempo_pha_raw) & length(tempo_dep_raw(~nan_pha)) >= 3
           nan_pha = isnan(tempo_pha_raw);
           if sum(nan_pha) ~= 0 
              pha_1d_in(nan_pha)=interp1(tempo_dep_raw(~nan_pha),tempo_pha_raw(~nan_pha),tempo_dep_raw(nan_pha));
           elseif sum(nan_pha) == 0 
               pha_1d_in = tempo_pha_raw; 
           end
       elseif sum(isnan(tempo_pha_raw)) == length(tempo_pha_raw) | length(tempo_dep_raw(~nan_pha)) < 3
           pha_1d_in = NaN(1,length(std_depth));
       end
       
       
       % interp to make it standard depth form.
       clearvars nan_* *_1d_f
       
       nan_temp = isnan(temp_1d_in)
       if sum(isnan(temp_1d_in)) == length(temp_1d_in) | sum(isnan(tp_1d_in)) == length(tp_1d_in)
          temp_1d_f = NaN(length(std_depth),1);
       elseif sum(isnan(temp_1d_in)) ~= length(temp_1d_in)
          temp_1d_f = interp1(tempo_dep_raw(~nan_temp),temp_1d_in(~nan_temp),std_depth);
       end
       
       nan_salt = isnan(salt_1d_in) 
       if sum(isnan(salt_1d_in)) == length(salt_1d_in) | sum(isnan(tp_1d_in)) == length(tp_1d_in)
          salt_1d_f = NaN(length(std_depth),1);
       elseif sum(isnan(salt_1d_in)) ~= length(salt_1d_in)
          salt_1d_f = interp1(tempo_dep_raw(~nan_salt),salt_1d_in(~nan_salt),std_depth);
       end
       
%        salt_1d_f =interp1(tempo_dep_raw(~nan_salt),salt_1d_in(~nan_salt),std_depth);
       
       nan_tp = isnan(tp_1d_in);
       if sum(isnan(tp_1d_in)) == length(tp_1d_in)
          tp_1d_f = NaN(length(std_depth),1);
       elseif sum(isnan(tp_1d_in)) ~= length(tp_1d_in)
          tp_1d_f = interp1(tempo_dep_raw(~nan_tp),tp_1d_in(~nan_tp),std_depth);
       end       
       
%        tp_1d_f =interp1(tempo_dep_raw(~nan_tp),tp_1d_in(~nan_tp),std_depth);
       
       nan_no3 = isnan(no3_1d_in);
       if sum(isnan(no3_1d_in)) == length(no3_1d_in)
          no3_1d_f = NaN(length(std_depth),1);
       elseif sum(isnan(no3_1d_in)) ~= length(no3_1d_in)
          no3_1d_f = interp1(tempo_dep_raw(~nan_no3),no3_1d_in(~nan_no3),std_depth);
       end
       
%        no3_1d_f =interp1(tempo_dep_raw(~nan_no3),no3_1d_in(~nan_no3),std_depth);
       
       nan_no2 = isnan(no2_1d_in);
       if sum(isnan(no2_1d_in)) == length(no2_1d_in)
          no2_1d_f = NaN(length(std_depth),1);
       elseif sum(isnan(no2_1d_in)) ~= length(no2_1d_in)
          no2_1d_f = interp1(tempo_dep_raw(~nan_no2),no2_1d_in(~nan_no2),std_depth);
       end
%        no2_1d_f =interp1(tempo_dep_raw(~nan_no2),no2_1d_in(~nan_no2),std_depth);
       
       nan_nh3 = isnan(nh3_1d_in);
       if sum(isnan(nh3_1d_in)) == length(nh3_1d_in)
          nh3_1d_f = NaN(length(std_depth),1);
       elseif sum(isnan(nh3_1d_in)) ~= length(nh3_1d_in)
          nh3_1d_f = interp1(tempo_dep_raw(~nan_nh3),nh3_1d_in(~nan_nh3),std_depth);
       end
%        nh3_1d_f =interp1(tempo_dep_raw(~nan_nh3),nh3_1d_in(~nan_nh3),std_depth);
       
       nan_ph = isnan(ph_1d_in);
       if sum(isnan(ph_1d_in)) == length(ph_1d_in)
          ph_1d_f = NaN(length(std_depth),1);
       elseif sum(isnan(ph_1d_in)) ~= length(ph_1d_in)
          ph_1d_f = interp1(tempo_dep_raw(~nan_ph),ph_1d_in(~nan_ph),std_depth);
       end
%        ph_1d_f =interp1(tempo_dep_raw(~nan_ph),ph_1d_in(~nan_ph),std_depth);
       
       nan_chl = isnan(chl_1d_in);
       if sum(isnan(chl_1d_in)) == length(chl_1d_in)
          chl_1d_f = NaN(length(std_depth),1);
       elseif sum(isnan(chl_1d_in)) ~= length(chl_1d_in)
          chl_1d_f = interp1(tempo_dep_raw(~nan_chl),chl_1d_in(~nan_chl),std_depth);
       end
%        chl_1d_f =interp1(tempo_dep_raw(~nan_chl),chl_1d_in(~nan_chl),std_depth);
       
       nan_pha = isnan(pha_1d_in);
       if sum(isnan(pha_1d_in)) == length(pha_1d_in)
          pha_1d_f = NaN(length(std_depth),1);
       elseif sum(isnan(pha_1d_in)) ~= length(pha_1d_in)
          pha_1d_f = interp1(tempo_dep_raw(~nan_pha),pha_1d_in(~nan_pha),std_depth);
       end
%        pha_1d_f =interp1(tempo_dep_raw(~nan_pha),pha_1d_in(~nan_pha),std_depth);
temp_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)= nanmean([squeeze(temp_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)) , temp_1d_f],2);
salt_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)= nanmean([squeeze(salt_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)) , salt_1d_f],2);
tp_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)= nanmean([squeeze(tp_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)) , tp_1d_f],2);
no3_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)= nanmean([squeeze(no3_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)) , no3_1d_f],2);
no2_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)= nanmean([squeeze(no2_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)) , no2_1d_f],2);
nh3_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)= nanmean([squeeze(nh3_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)) , nh3_1d_f],2);
ph_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)= nanmean([squeeze(ph_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)) , ph_1d_f],2);
chl_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)= nanmean([squeeze(chl_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)) , chl_1d_f],2);
pha_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)= nanmean([squeeze(pha_3d(near_point_2d(i,1),near_point_2d(i,2),:,indx_raw)) , pha_1d_f],2);
    
end
    
end

tz = tz + 1
time_3d(tz) = indx_raw; %when there is a data
end

kk =0
for i = 1:size(temp_3d,4)
    if length(find(isnan(squeeze(temp_3d(:,:,1,i)))==0)) ~= 0
        kk = kk +1;
        ni(kk) = i
    end
end

kk =0
for i = 1:size(nh3_3d,4)
    if length(find(isnan(squeeze(nh3_3d(:,:,1,i)))==0)) ~= 0
        kk = kk +1;
        ni_nh3(kk) = i
    end
end

kk =0
for i = 1:size(no3_3d,4)
    if length(find(isnan(squeeze(no3_3d(:,:,1,i)))==0)) ~= 0
        kk = kk +1;
        ni_no3(kk) = i
    end
end
cd C:\Users\user\Desktop\jma
end
save('obs_jma_mat')

% check which file blow-up (mainly depth blow it up)
list_file(f).name
tempo_dep_raw(~nan_no3)



% nanmean(tempo_dep_raw(~nan_no3)ic

figure
pcolor(ref_lon_m, ref_lat_m, squeeze(temp_3d(:,:,1,ni(1))));
colorbar;
ylim([min(min(ref_lat_m)) max(max(ref_lat_m))])
xlim([min(min(ref_lon_m)) max(max(ref_lon_m))])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain')  % overlay google map

figure
pcolor(ref_lon_m, ref_lat_m, squeeze(temp_3d(:,:,1,ni(10))));
colorbar;
ylim([min(min(ref_lat_m)) max(max(ref_lat_m))])
xlim([min(min(ref_lon_m)) max(max(ref_lon_m))])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','satellite')  % overlay google map


figure
pcolor(ref_lon_m, ref_lat_m, squeeze(no3_3d(:,:,1,ni_no3(1))));
colorbar;
ylim([min(min(ref_lat_m)) max(max(ref_lat_m))])
xlim([min(min(ref_lon_m)) max(max(ref_lon_m))])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','satellite')  % overlay google map



figure
pcolor(ref_lon_m, ref_lat_m, squeeze(no3_3d(:,:,1,ni_nh3(1))));
colorbar;


figure
pcolor(ref_lon_m, ref_lat_m, squeeze(salt_3d(:,:,1,indx_raw)));
colorbar;
        
% figure
% pcolor(ref_lon_m, ref_lat_m, squeeze(nh3_3d(:,:,2,indx_raw)));
% colorbar;   
% 
% pcolor(ref_lon_m, ref_lat_m, squeeze(nh3_3d(:,:,3,indx_raw)));
% colorbar;  
        
% for i = 1:length(lon_raw_decimal)
%     near_point_2d(1,1),near_point_2d(1,2)
% 
% 
% end





% file{2}(1:7)
% 
% length(strfind(file_cont{2},'AF-'))



% HEADER-2 (ö´ïÃï×ÜÃ)
% é©áÈ ËÒã· «Õ«£?«ë«É «Õ«£?«ë«ÉªÎ?Ù¥
% êÈöÇ «¿«¤«×
% STATION NO 1 
% LATITUDE 9     deg,min
% LONGITUDE 17 
% DATE/TIME 26   ?ö´ªÎËÒã·ÐàªÓðûÖõªÎêÅ-ìí-ãÁ£¨ìíÜâøöñÞãÁ£© start / end
% W-DEPTH 48 I4 ú­î¼ªÞªÇªÎâ©ä¢ (m)  -> surface to bottom
% W-COLOR 54 I2 â©ßä£¨«Õ«©?«ì«ë?«¦?«ìøöñÞäûªÎÛã?£© -> water color
% TRANS 57  ÷âÙ¥Óø÷ùªËªèªë÷âÙ¥Óø(m?êÈ)ªÈ«ï«¤«ä?ªÎÌËÊÇ£¨Óø) -> secchi depth
% SSF-NO 102    ??ª¹ªëøúöµâ©?«Ç?«¿ªÎö´ïÃÛã?  ->  sst data obs station
% ACM-NO 109    ??ª¹ªëøúöµú­×µ«Ç?«¿ªÎö´ïÃÛã?  ->  ss-current data obs station
% SUB STN NO 116     ÜÍð¾ö´ïÃÛã?  -> sub- station num

% DATA (?ö´«Ç?«¿)
% é©áÈ ËÒã· «Õ«£?«ë«É «Õ«£?«ë«ÉªÎ?Ù¥
% êÈöÇ «¿«¤«×
% STATION NO 1 A3,I4 ö´ïÃÛã?£¨ÊÀöµ?ö´«³?«É + Ö§?ª·ª¿ 4 ùùªÎ?í®)
% TIME 9 2I2 óõâ©ãÁÊ¾ (ìíÜâøöñÞãÁ)
% DEPTH(OBS) 17 I4 óõâ©öµªÎä¢Óø (m)
% TEMP(OBS) 22 F6.3 CTD â©? (1990 Ò´?ð·?ÓøÙÍàü (ITS-90))
% SAL(OBS) 28 F6.3 CTD ?ÝÂ (1978 Ò´?éÄ?ÝÂ (PSS-78))
% DO 35 I3 éÁðíß«áÈÒØÓø£¨?mol/l)
% PO4-P 39 F4.2 «ê«óß«?£¨?mol/l)
% T-P 44 F4.2 îï«ê«óß«? (?mol/l)
% NO3-N 49 F4.1 õ¦ß«?£¨?mol/l£©
% NO2-N 54 F4.2 ?õ¦ß«?£¨?mol/l£©
% NH3-N 59 F4.2 «¢«ó«â«Ë«¢ (?mol/l)
% PH 64 F4.2 25 ¡ÉªËªªª±ªëâ©áÈ«¤«ª«óÒØÓøò¦?£¨NBS «¹«±?«ë£©
% CHL 69 F6.2 «¯«í«í«Õ«£«ë a (?g/l)
% PHA 76 F6.2 «Õ«£«ª«Õ«£«Á«ó (?g/l)
% (ADD PARAM) 83 (õÚÊ¥é©áÈ)
% ¡¸PRESSURE¡¹?Õô£¨104Pa), ¡¸COD¡¹ûù?îÜß«áÈé©
% Ï´Õá (mg/l), ¡¸SILICATE¡¹«±«¤ß«? (?mol/l),
% ¡¸TOTAL-N¡¹îïòòáÈ (?mol/l), ¡¸ALKALINITY¡¹«¢«ë
% «««êÓø (mmol/l), ¡¸TIC¡¹îï÷©ß« (mmol/l)
% DEPTH(STD) 94 I4 øöñÞöµ£¨m£©
% TEMP(STD) 99 F6.3 CTD ?ö´ªËªèªëâ©? (ITS-90)
% SAL(STD) 105 F6.3 CTD ?ö´ªËªèªë?ÝÂ (PSS-78)
% D-ST 116 I4 «µ?«â«¹«Æ«ê«Ã«¯«¢«Î«Þ«ê?£¨10?8m
% 3/kg£©
% DELTA-D 121 F5.3 «¸«ª«Ý«Æ«ó«·«ã«ë«¢«Î«Þ«ê? (10m2/sec2)
% REC IND 126 A1 ¡®@¡¯(ÊÀö´ïÃªÎÑÀ?ªÎðûÖõªòãÆª¹) ªÞª¿ªÏ ¡®=¡¯
% REC IND 126 A1 @, =  end of obs

% file{5}(22:27)
% str2num(file{5}(28:34))