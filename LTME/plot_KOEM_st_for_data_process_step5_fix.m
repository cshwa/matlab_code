close all; clear; clc;   % -v3

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
for i = 1:3
name_tag_1{i} = ['여수신항 H' num2str(i)] 
end

name_tag_2{1} = ['광양항 H' num2str(1)] 

name_tag_3{1} = ['삼천포항 H' num2str(1)] 

% estuary
for i = 1:5
name_tag_4{i} = ['가막만 ' num2str(i,'%02d')] 
end

for i = 1:25
name_tag_5{i} = ['섬진강하구 ' num2str(i,'%02d')] 
end

for i = 1:2
name_tag_6{i} = ['진주만 ' num2str(i,'%02d')] 
end

%coastal 
for i = 1:28
name_tag_7{i} = ['대한해협연안 ' num2str(i,'%02d')] 
end

% combining the tag and outter point excluding
name_tag = name_tag_1'; 
name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
size_tag = length(name_tag);

for i = 1:length(name_tag_4)
name_tag{size_tag+i} = name_tag_4{i}; 
end
size_tag = length(name_tag);

for i = 1:length(name_tag_5)
name_tag{size_tag+i} = name_tag_5{i}; 
end
size_tag = length(name_tag);

for i = 1:length(name_tag_6)
name_tag{size_tag+i} = name_tag_6{i}; 
end
size_tag = length(name_tag);

% %  when skip the outer point
% % for i = 1:length(name_tag_7)-18  % 연안 3, 6, 13:28 has to be remove (18)  [it's located at out of domain] 
% %  new_i = [1:28]; new_i(3)=[]; new_i(6)=[]; new_i(1)=[];
% %  name_tag{size_tag+i} = name_tag_7{new_i(i)};
% % end

for i = 1:length(name_tag_7) % 연안 3, 6, 13:28 has to be remove (18)  [it's located at out of domain] 
 name_tag{size_tag+i} = name_tag_7{i};
end

%% pick the row on the excel which has same name with tag
[raw_p txt_p]=xlsread('해양환경측정망(항만측정망).xls','sheet','');
txt_matc_p = txt_p(3:end,1); % name list
txt_date_p = txt_p(3:end,3); % date list
temp_sur_p = txt_p(3:end,4); 
temp_bot_p = txt_p(3:end,5); 
salt_sur_p = txt_p(3:end,6); 
salt_bot_p = txt_p(3:end,7);
do_sur_p = txt_p(3:end,10); 
do_bot_p = txt_p(3:end,11);
nh4_sur_p = txt_p(3:end,14); 
nh4_bot_p = txt_p(3:end,15);
no3_sur_p = txt_p(3:end,18); 
no3_bot_p = txt_p(3:end,19);
chl_sur_p = txt_p(3:end,32); 
chl_bot_p = txt_p(3:end,33);

[raw_es txt_es]=xlsread('해양환경측정망(하천영향및반폐쇄성해역환경측정망).xls','sheet','');
txt_matc_es = txt_es(3:end,1);
txt_date_es = txt_es(3:end,3);
temp_sur_es = txt_es(3:end,4); 
temp_bot_es = txt_es(3:end,5); 
salt_sur_es = txt_es(3:end,6); 
salt_bot_es = txt_es(3:end,7); 
do_sur_es = txt_es(3:end,10); 
do_bot_es = txt_es(3:end,11);
nh4_sur_es = txt_es(3:end,14); 
nh4_bot_es = txt_es(3:end,15);
no3_sur_es = txt_es(3:end,18); 
no3_bot_es = txt_es(3:end,19);
chl_sur_es = txt_es(3:end,32); 
chl_bot_es = txt_es(3:end,33);

[raw_co txt_co]=xlsread('해양환경측정망(연안해역환경측정망).xls','sheet','');
txt_matc_co = txt_co(3:end,1);
txt_date_co = txt_co(3:end,3);
temp_sur_co = txt_co(3:end,4); 
temp_bot_co = txt_co(3:end,5); 
salt_sur_co = txt_co(3:end,6); 
salt_bot_co = txt_co(3:end,7);
do_sur_co = txt_co(3:end,10); 
do_bot_co = txt_co(3:end,11);
nh4_sur_co = txt_co(3:end,14); 
nh4_bot_co = txt_co(3:end,15);
no3_sur_co = txt_co(3:end,18); 
no3_bot_co = txt_co(3:end,19);
chl_sur_co = txt_co(3:end,32); 
chl_bot_co = txt_co(3:end,33);

merge_txt = [txt_matc_p; txt_matc_es; txt_matc_co;]; % name list
merge_date = [txt_date_p; txt_date_es; txt_date_co;]; % date list
merge_temp_sur = [temp_sur_p; temp_sur_es; temp_sur_co;];
merge_temp_bot = [temp_bot_p; temp_bot_es; temp_bot_co;];
merge_salt_sur = [salt_sur_p; salt_sur_es; salt_sur_co;];
merge_salt_bot = [salt_bot_p; salt_bot_es; salt_bot_co;];
merge_do_sur = [do_sur_p; do_sur_es; do_sur_co;];
merge_do_bot = [do_bot_p; do_bot_es; do_bot_co;];
merge_nh4_sur = [nh4_sur_p; nh4_sur_es; nh4_sur_co;];
merge_nh4_bot = [nh4_bot_p; nh4_bot_es; nh4_bot_co;];
merge_no3_sur = [no3_sur_p; no3_sur_es; no3_sur_co;];
merge_no3_bot = [no3_bot_p; no3_bot_es; no3_bot_co;];
merge_chl_sur = [chl_sur_p; chl_sur_es; chl_sur_co;];
merge_chl_bot = [chl_bot_p; chl_bot_es; chl_bot_co;];

merge_data_txt = [txt_p(3:end,4:end); txt_es(3:end,4:end); txt_co(3:end,4:end);];

%% pick matched name with tag
% %port
% for i = 1:length(name_tag)
%    if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
%        indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)     
%    end
% end

%merge
for i = 1:length(name_tag)
   if  sum(strcmp(name_tag{i}, merge_txt)) ~= 0
       indx{i} = find([strcmp(name_tag{i}, merge_txt)] == 1)     
   end
end

%% make date to be 'yymm' form
for i = 1:length(merge_date)
temp = char(merge_date{i});
if size(temp) ~= 7
    temp=temp(1,1:7);
end
merge_yymm{i,1} = temp;
end

%% make date to be 'mm' form
for i = 1:length(merge_date)
temp = char(merge_date{i});
temp = temp(1,6:7);
merge_mm{i,1} = temp;
end

%% make 1997 to 2018 'yymm' form
k=0
for i = 1997:2018
    for j = 1:12
        k=k+1;
        ref_date{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end

%% make 1997 to 2018 'yymm' form
for j = 1:12
 ref_date_mm{j,1} = [num2str(j,'%02d')];
end

% matched date 'yymm' form
for j = 1:length(indx) % st. axis
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_yymm(indx{j}))) ~= 0
           indx_date{j,i} = find([strcmp(ref_date{i}, merge_yymm(indx{j}))] == 1);     
       end
    end
end

% matched date 'mm' form
for j = 1:length(indx) % st. axis
    for i = 1:length(ref_date_mm) % date axis
       if  sum(strcmp(ref_date_mm{i}, merge_mm(indx{j}))) ~= 0
           indx_date_mm{j,i} = find([strcmp(ref_date_mm{i}, merge_mm(indx{j}))] == 1);     
       end
    end
end


% make climate
%temp
clearvars temp
for i = 1:length(indx)
    temp = merge_temp_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        temp_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_temp_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        temp_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

%salt
clearvars temp
for i = 1:length(indx)
    temp = merge_salt_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        salt_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_salt_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        salt_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

%do
clearvars temp
for i = 1:length(indx)
    temp = merge_do_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        do_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_do_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        do_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

%nh4
clearvars temp
for i = 1:length(indx)
    temp = merge_nh4_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        nh4_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_nh4_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        nh4_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

%no3
clearvars temp
for i = 1:length(indx)
    temp = merge_no3_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        no3_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_no3_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        no3_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

%chl
clearvars temp
for i = 1:length(indx)
    temp = merge_chl_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        chl_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_chl_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        chl_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end


lon=ncread('grid_sumjin_v1970_fix_3m.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m.nc','mask_rho');
load KOEM_st_info_(decimal_deg).mat
lat_koem = [lat_1; lat_2; lat_3;];
lon_koem = [lon_1; lon_2; lon_3;];

%% interp

%temp
clearvars temp
for i = 1:size(temp_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,temp_sur_clim(:,i),lon,lat);
    temp_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(temp_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,temp_bot_clim(:,i),lon,lat);
    temp_bot_clim_g(:,:,i) = temp;
    clearvars temp
end

%salt
clearvars temp
for i = 1:size(salt_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,salt_sur_clim(:,i),lon,lat);
    salt_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(salt_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,salt_bot_clim(:,i),lon,lat);
    salt_bot_clim_g(:,:,i) = temp;
    clearvars temp
end

%do
clearvars temp
for i = 1:size(do_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,do_sur_clim(:,i),lon,lat);
    do_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(do_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,do_bot_clim(:,i),lon,lat);
    do_bot_clim_g(:,:,i) = temp;
    clearvars temp
end

%nh4
clearvars temp
for i = 1:size(nh4_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,nh4_sur_clim(:,i),lon,lat);
    nh4_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(nh4_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,nh4_bot_clim(:,i),lon,lat);
    nh4_bot_clim_g(:,:,i) = temp;
    clearvars temp
end

%no3
clearvars temp
for i = 1:size(no3_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,no3_sur_clim(:,i),lon,lat);
    no3_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(no3_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,no3_bot_clim(:,i),lon,lat);
    no3_bot_clim_g(:,:,i) = temp;
    clearvars temp
end

%chl
clearvars temp
for i = 1:size(chl_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,chl_sur_clim(:,i),lon,lat);
    chl_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(chl_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,chl_bot_clim(:,i),lon,lat);
    chl_bot_clim_g(:,:,i) = temp;
    clearvars temp
end


%% check
% for i = 1: 11; figure; pcolor(lon,lat,temp_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,salt_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,do_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,nh4_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,no3_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,chl_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% 
% for i = 1: 11; figure; pcolor(lon,lat,temp_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,salt_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,do_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,nh4_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,no3_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,chl_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end

%% filling nan grid on the ROMS

%temp
clearvars temp
for i = 1:size(temp_sur_clim,2)
    temp=squeeze(temp_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        temp_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        temp_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(temp_bot_clim,2)
    temp=squeeze(temp_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        temp_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        temp_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

%salt
clearvars temp
for i = 1:size(salt_sur_clim,2)
    temp=squeeze(salt_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        salt_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        salt_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(salt_bot_clim,2)
    temp=squeeze(salt_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        salt_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        salt_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end


%do
clearvars temp
for i = 1:size(do_sur_clim,2)
    temp=squeeze(do_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        do_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        do_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(do_bot_clim,2)
    temp=squeeze(do_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        do_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        do_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end


%nh4
clearvars temp
for i = 1:size(nh4_sur_clim,2)
    temp=squeeze(nh4_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        nh4_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        nh4_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(nh4_bot_clim,2)
    temp=squeeze(nh4_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        nh4_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        nh4_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

%no3
clearvars temp
for i = 1:size(no3_sur_clim,2)
    temp=squeeze(no3_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        no3_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        no3_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(no3_bot_clim,2)
    temp=squeeze(no3_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        no3_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        no3_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

%chl
clearvars temp
for i = 1:size(chl_sur_clim,2)
    temp=squeeze(chl_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        chl_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        chl_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(chl_bot_clim,2)
    temp=squeeze(chl_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        chl_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        chl_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

save('koem_climate.mat', '*_clim_g_f');

clearvars temp
for i = 1:size(lon,1)
    for j = 1:size(lon,2)
            % temp
            temp_t_s=squeeze(temp_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_t_b=squeeze(temp_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % salt
            temp_s_s=squeeze(salt_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_s_b=squeeze(salt_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % do
            temp_d_s=squeeze(do_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_d_b=squeeze(do_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % nh4
            temp_nh_s=squeeze(nh4_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_nh_b=squeeze(nh4_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % no3
            temp_no_s=squeeze(no3_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_no_b=squeeze(no3_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % chl
            temp_ch_s=squeeze(chl_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_ch_b=squeeze(chl_bot_clim_g_f(i,j,:)); % pick time series of one point

     
     temp_t_s = [temp_t_s(end); NaN; temp_t_s; NaN; NaN; temp_t_s(2);];
     temp_t_b = [temp_t_b(end); NaN; temp_t_b; NaN; NaN; temp_t_b(2);];
     
     temp_s_s = [temp_s_s(end); NaN; temp_s_s; NaN; NaN; temp_s_s(2);];
     temp_s_b = [temp_s_b(end); NaN; temp_s_b; NaN; NaN; temp_s_b(2);];
     
     temp_d_s = [temp_d_s(end); NaN; temp_d_s; NaN; NaN; temp_d_s(2);];
     temp_d_b = [temp_d_b(end); NaN; temp_d_b; NaN; NaN; temp_d_b(2);];
     
     temp_nh_s = [temp_nh_s(end); NaN; temp_nh_s; NaN; NaN; temp_nh_s(2);];
     temp_nh_b = [temp_nh_b(end); NaN; temp_nh_b; NaN; NaN; temp_nh_b(2);];
     
     temp_no_s = [temp_no_s(end); NaN; temp_no_s; NaN; NaN; temp_no_s(2);];
     temp_no_b = [temp_no_b(end); NaN; temp_no_b; NaN; NaN; temp_no_b(2);];
     
     temp_ch_s = [temp_ch_s(end); NaN; temp_ch_s; NaN; NaN; temp_ch_s(2);];
     temp_ch_b = [temp_ch_b(end); NaN; temp_ch_b; NaN; NaN; temp_ch_b(2);];
     
     % filling NaNs on the time order
     t = 1:length(temp_t_s); %general time
     temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
     temp_t_b(isnan(temp_t_b)) = interp1(t(~isnan(temp_t_b)), temp_t_b(~isnan(temp_t_b)), t(isnan(temp_t_b)) ); 
     temp_t_s(1:2)=[]; temp_t_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_t_b(1:2)=[]; temp_t_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_s_s(isnan(temp_s_s)) = interp1(t(~isnan(temp_s_s)), temp_s_s(~isnan(temp_s_s)), t(isnan(temp_s_s)) ); 
     temp_s_b(isnan(temp_s_b)) = interp1(t(~isnan(temp_s_b)), temp_s_b(~isnan(temp_s_b)), t(isnan(temp_s_b)) );
     temp_s_s(1:2)=[]; temp_s_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_s_b(1:2)=[]; temp_s_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_d_s(isnan(temp_d_s)) = interp1(t(~isnan(temp_d_s)), temp_d_s(~isnan(temp_d_s)), t(isnan(temp_d_s)) ); 
     temp_d_b(isnan(temp_d_b)) = interp1(t(~isnan(temp_d_b)), temp_d_b(~isnan(temp_d_b)), t(isnan(temp_d_b)) );
     temp_d_s(1:2)=[]; temp_d_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_d_b(1:2)=[]; temp_d_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_nh_s(isnan(temp_nh_s)) = interp1(t(~isnan(temp_nh_s)), temp_nh_s(~isnan(temp_nh_s)), t(isnan(temp_nh_s)) ); 
     temp_nh_b(isnan(temp_nh_b)) = interp1(t(~isnan(temp_nh_b)), temp_nh_b(~isnan(temp_nh_b)), t(isnan(temp_nh_b)) );
     temp_nh_s(1:2)=[]; temp_nh_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_nh_b(1:2)=[]; temp_nh_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_no_s(isnan(temp_no_s)) = interp1(t(~isnan(temp_no_s)), temp_no_s(~isnan(temp_no_s)), t(isnan(temp_no_s)) ); 
     temp_no_b(isnan(temp_no_b)) = interp1(t(~isnan(temp_no_b)), temp_no_b(~isnan(temp_no_b)), t(isnan(temp_no_b)) );
     temp_no_s(1:2)=[]; temp_no_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_no_b(1:2)=[]; temp_no_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_ch_s(isnan(temp_ch_s)) = interp1(t(~isnan(temp_ch_s)), temp_ch_s(~isnan(temp_ch_s)), t(isnan(temp_ch_s)) ); 
     temp_ch_b(isnan(temp_ch_b)) = interp1(t(~isnan(temp_ch_b)), temp_ch_b(~isnan(temp_ch_b)), t(isnan(temp_ch_b)) );
     temp_ch_s(1:2)=[]; temp_ch_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_ch_b(1:2)=[]; temp_ch_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
      % interp on the layer order  

     prepretemp = [NaN(20,12)];
     prepretemp(20,1:12) = temp_t_s'; prepretemp(1,1:12) = temp_t_b';
     
     prepresalt = [NaN(20,12)];
     prepresalt(20,1:12) = temp_s_s'; prepresalt(1,1:12) = temp_s_b';
     
     prepredo = [NaN(20,12)];
     prepredo(20,1:12) = temp_d_s'; prepredo(1,1:12) = temp_d_b';
     
     preprenh4 = [NaN(20,12)];
     preprenh4(20,1:12) = temp_nh_s'; preprenh4(1,1:12) = temp_nh_b';
     
     prepreno3 = [NaN(20,12)];
     prepreno3(20,1:12) = temp_no_s'; prepreno3(1,1:12) = temp_no_b';
     
     preprech = [NaN(20,12)];
     preprech(20,1:12) = temp_ch_s'; preprech(1,1:12) = temp_ch_b';
     
          for k = 1:12 %time
              temp_t = squeeze(prepretemp(:,k));
              temp_s = squeeze(prepresalt(:,k));
              temp_do =  squeeze(prepredo(:,k));
              temp_nh4 =  squeeze(preprenh4(:,k));
              temp_no3 =  squeeze(prepreno3(:,k));
              temp_ch =  squeeze(preprech(:,k));
              
              g_t = 1:length(temp_t); %general layer num.
              temp_t(isnan(temp_t))=interp1(g_t(~isnan(temp_t)), temp_t(~isnan(temp_t)), g_t(isnan(temp_t)) );
              temp_s(isnan(temp_s))=interp1(g_t(~isnan(temp_s)), temp_s(~isnan(temp_s)), g_t(isnan(temp_s)) );
              temp_do(isnan(temp_do))=interp1(g_t(~isnan(temp_do)), temp_do(~isnan(temp_do)), g_t(isnan(temp_do)) );
              temp_nh4(isnan(temp_nh4))=interp1(g_t(~isnan(temp_nh4)), temp_nh4(~isnan(temp_nh4)), g_t(isnan(temp_nh4)) );
              temp_no3(isnan(temp_no3))=interp1(g_t(~isnan(temp_no3)), temp_no3(~isnan(temp_no3)), g_t(isnan(temp_no3)) );
              temp_ch(isnan(temp_ch))=interp1(g_t(~isnan(temp_ch)), temp_ch(~isnan(temp_ch)), g_t(isnan(temp_ch)) );
              
              pretemp(:,k) = temp_t;
              presalt(:,k) = temp_s;
              predo(:,k) = temp_do;
              preno4(:,k) = temp_nh4;
              preno3(:,k) = temp_no3;
              prech(:,k) = temp_ch;
              clearvars temp_t temp_s temp_do temp_nh4 temp_no3 temp_ch g_t
          end
               
          % completet
          temperature(i,j,:,:) = pretemp;
          salt(i,j,:,:) = presalt;
          do(i,j,:,:) = predo;
          nh4(i,j,:,:) = preno4;
          no3(i,j,:,:) = preno3;
          chla(i,j,:,:) = prech;
          
          clearvars temp_t* temp_s_* temp_d* temp_n* temp_ch*  pre* 
    end
end

return

save('KOEM_interped_climate.mat'); 

%% figure

%% temp
temp_range=[4:1:28];
for i = 1:12
figure%('position',[550 550 550 550],'PaperUnits','inches','PaperPosition',[0 0 7.5 3.5]);
% set(gca,'Position',[0.2 0.17 0.73 0.65]);
hold on;
% [cs,h]=contour(lon,lat,squeeze(temp_sur_clim_g(:,:,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
% clabel(cs,h,'fontsize',17,'labelspacing',700,'fontweight','bold');
pcolor(lon,lat,squeeze(temperature(:,:,20,i)).*(mask./mask));
% colormap summer(20);xlim([0 0.3])

H1=colorbar;
set(get(H1,'title'),'string','^oC','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([temp_range(1) temp_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_Temp_sur_' num2str(i) 'mth']; 

print('-dpng',filename);
end

temp_range=[4:1:28];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(temperature(:,:,1,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','^oC','fontsize',18,'fontweight','bold');
caxis([temp_range(1) temp_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_Temp_bot_' num2str(i) 'mth']; 

print('-dpng',filename);
end

%% salt
salt_range=[2:1:35];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(salt(:,:,20,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','PSU','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([salt_range(1) salt_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_salt_sur_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

salt_range=[2:1:35];
for i = 1:size(salt_bot_clim,2)
figure;
hold on;
pcolor(lon,lat,squeeze(salt(:,:,1,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','PSU','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([salt_range(1) salt_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_salt_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

%% DO
do_range=[4:1:14];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(do(:,:,20,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','DO-mg/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([do_range(1) do_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_do_sur_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

do_range=[4:1:14];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(do(:,:,1,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','DO-mg/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([do_range(1) do_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_do_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

%% NH4
% nh4_range=[0:1:135];
% nh4_range=[0:1:90];
nh4_range=[0:1:205];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(nh4(:,:,20,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','nh4-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([nh4_range(1) nh4_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_nh4_sur_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

nh4_range=[0:1:205];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(nh4(:,:,1,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','nh4-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([nh4_range(1) nh4_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_nh4_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

%% NO3
no3_range=[0:1:650];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(no3(:,:,20,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','no3-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([no3_range(1) no3_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_no3_sur_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

no3_range=[0:1:650];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(no3(:,:,1,i)).*(mask./mask));
% colormap summer(20);xlim([0 0.3])

H1=colorbar;
set(get(H1,'title'),'string','no3-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([no3_range(1) no3_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_no3_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 
end


%% chla
chl_range=[0:1:18];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(chla(:,:,20,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','chl-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([chl_range(1) chl_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_chl_sur_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

chl_range=[0:1:18];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(chla(:,:,1,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','chl-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([chl_range(1) chl_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_chl_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

close all

    % save('KOEM_name_tag.mat','name_tag');